%% CGM, this script loads the MCMC and ET results for post-processing of forecast stats

clear variables; close all; clc;
%#ok<*UNRCH>

%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%% settings
models            = {'STATEtrendgapSV', 'STATEtrendgapSVnoise2',  'STATEscaleSV', 'STATEconst'};
DATALABELS        = {'RGDP', 'PGDP', 'UNRATE'};

mcmc_results_folder = localresultsMCMC;
MCMCdraws           = 3e3;
ET_results_folder   = localresultsET;
SPF_data_folder     = fullfile('..', 'kensingtonDataMatfiles');
SPFprob             = load(fullfile(SPF_data_folder,'kensingtonPROB.mat')); % we only use this to check if correct ET and CGM files are loaded

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx90     = [2 8];
ndxQ1     = 4;
ndxQ2     = 5; % median
ndxQ3     = 6;

Kavg         = 4;

percentiles  = 1:99;
Npercentiles = length(percentiles);

doPlotSingle = false;

initwrap

%% allocate memory
T      = 215; % note: hardcoded 
RELESS = NaN(T, 2, length(models), length(DATALABELS));
ESS    = NaN(T, 2, length(models), length(DATALABELS));
omega_x= NaN(2, 2, T, length(models), length(DATALABELS));
KLIC   = NaN(T, 2, length(models), length(DATALABELS));
zero_weight_all = NaN(T, 2, length(models), length(DATALABELS));
large_weight_counter_all = NaN(T, 2, length(models), length(DATALABELS));
large_weight_sum_all = NaN(T, 2, length(models), length(DATALABELS));
small_weight_counter_all = NaN(T, 2, length(models), length(DATALABELS));

ESS_threshold           = 1000; % threshold for ESS, "arbitrary"
large_weight_threshold  = 0.01; % threshold for "large" weights, "arbitrary"

%% loop over ET modes
for doBinsOnly = [true, false]

    if doBinsOnly
        ETlabel = 'binsOnly';
        ndxET   = 1;
    else
        ETlabel = 'binsAndMeans';
        ndxET   = 2;
    end

    %% start loop over models
    for m = 1 : length(models)

        modellabel = models{m};

        %% start loop over variables
        for thisD = 1 : length(DATALABELS)

            datalabel = DATALABELS{thisD};

            % set up log file
            diaryfilename=fullfile(localstore, sprintf('inspectweights-%s-%s-ET%s-Ndraws%d.log', modellabel, datalabel, ETlabel, MCMCdraws));
            if exist(diaryfilename, 'file')
                delete(diaryfilename);
            end
            diary(diaryfilename);

            % select the SPF results corresponding to the variable
            SPFprob_act = SPFprob.(datalabel);

            % load raw CGM results

            CGMmcmc = matfile(fullfile(mcmc_results_folder, sprintf('CGMmcmc-%s%s-Ndraws%d.mat', datalabel, modellabel, MCMCdraws)));
            Ndraws  = CGMmcmc.MCMCdraws * CGMmcmc.Nfedraws;
            init_w  = ones(Ndraws,1)/Ndraws; % initial weight vector, weights are uniform, so all 1/Ndraws

            dates         = CGMmcmc.dates;
            Zdata         = CGMmcmc.Zdata;
            Znanny        = CGMmcmc.Znanny; % used only for inclusion in result file
            Ylabel        = CGMmcmc.Ylabel; % used only for inclusion in result file
            Yfuture       = CGMmcmc.Yfuture;


            T          = CGMmcmc.T;
            Nhorizons  = CGMmcmc.Nhorizons;
            Nbar       = CGMmcmc.Nbar;
            Nquantiles = CGMmcmc.Nquantiles;
            if ~isequal(quantileP, CGMmcmc.quantileP)
                error('mismatch w/assumed values of quantileP')
            end


            % find subsample where ET is applied
            firstT = max(find(~isnan(SPFprob_act.probs(:,1)),1,'first'),CGMmcmc.Tstart);
            lastT  = min(find(~isnan(SPFprob_act.probs(:,1)),1,'last'), CGMmcmc.T);
            numT   = lastT - firstT + 1;
            first_last_Date = dates([firstT;lastT]);

            first_last_DateLabel = datestr(first_last_Date, 'yyyyqq');



            %% clear CGMmcmc
            clear CGMmcmc

            %% load ET results
            ET_results = matfile(fullfile(ET_results_folder,sprintf('kensington-ET%s-Ndraws%d-%s%s.mat', ETlabel, MCMCdraws, datalabel, modellabel)));
            ET_w_all   = permute(ET_results.ET_w_all, [2 1]);
            ET_w_all   = ET_w_all(:,firstT:lastT);

            % check if ET subsample is the same
            if first_last_Date ~= ET_results.first_last_Date
                error('ET dates not matching, check if matching CGM and ET file loaded!')
            end

            %% collect RELESS
            allweights = sort(ET_w_all);
            reless     = 1 ./ sum(allweights.^2, 1) ./ Ndraws;
            ess        = 1 ./ sum(allweights.^2, 1);
            RELESS(firstT:lastT, ndxET, m, thisD) = reless;
            ESS(firstT:lastT, ndxET, m, thisD)    = ess;
            small_weight_threshold = (1/Ndraws)*1e-03;

            %% plot individual figure
            if doPlotSingle
                ETdates     = dates(firstT:lastT);
                thisfig     = figure;
                plot(ETdates, reless)
                ylim([0 1])
                xtickdates(ETdates)
                title(sprintf('%s\n%s %s', datalabel, ETlabel, modellabel))
                wrapthisfigure(thisfig, sprintf('RELESS-%s-%s-%s', datalabel, modellabel, ETlabel), wrap);
            end

            %% calculate omega_x as in Cogley, Morozov and Sargent (2005, JEDC, p.1917)

            twist_x_vec = [10,10*Ndraws/5000]';
            % CMS had 5000 draws and calculated relative weight of 10 largest, we have Ndraws draws, so it might make sense to scale up

            for twist_x_idx = 1 : length(twist_x_vec)

                twist_x = twist_x_vec(twist_x_idx);
                omega_x(twist_x_idx,ndxET,firstT:lastT,m,thisD) = (Ndraws/twist_x) * sum((allweights(end-twist_x+1:end,:) * Ndraws).^2) ./ sum((allweights * Ndraws).^2);

            end

            %% calculate KLIC, issue warnings for zero and large weights
            log_w_unif = -log(Ndraws);
        
            for t_idx = 1 : lastT - firstT + 1
                w_tmp = allweights(:,t_idx);
                log_w_diff_tmp = (log(w_tmp)-log_w_unif);
                KLIC_tmp = 0;
                zero_weight_counter     = 0;
                large_weight_counter    = 0;
                large_weight_sum        = 0;
                small_weight_counter    = 0;
                for ii = 1 : Ndraws
                    if w_tmp(ii,1) == 0 % while weights cannot be zero theoretically, numerically it might happen, then we disregard that draw and weight and skip to the next one
                        zero_weight_counter = zero_weight_counter + 1;
                        continue
                    else
                        KLIC_tmp = KLIC_tmp + w_tmp(ii,1)*log_w_diff_tmp(ii,1);
                    end
                    
                    if w_tmp(ii,1) > large_weight_threshold
                        large_weight_counter = large_weight_counter + 1;
                        large_weight_sum     = large_weight_sum + w_tmp(ii,1);
                    end

                    if w_tmp(ii,1) < small_weight_threshold
                        small_weight_counter = small_weight_counter + 1;
                    end


                end

                if zero_weight_counter > 0
                    warning('%d zero weights %s %s %s %s ',zero_weight_counter, datestr(dates(firstT+t_idx-1),'yyyyQQ'), datalabel, modellabel, ETlabel);
                end

                zero_weight_all(firstT+t_idx-1, ndxET, m, thisD) = zero_weight_counter;

                if large_weight_counter > 0
                    warning('%d large weights with sum %1.2f %s %s %s %s ',large_weight_counter, large_weight_sum, datestr(dates(firstT+t_idx-1),'yyyyQQ'), datalabel, modellabel, ETlabel);
                end

                large_weight_counter_all(firstT+t_idx-1, ndxET, m, thisD) = large_weight_counter;
                large_weight_sum_all(firstT+t_idx-1, ndxET, m, thisD) = large_weight_sum;

                small_weight_counter_all(firstT+t_idx-1, ndxET, m, thisD) = small_weight_counter;

                KLIC(firstT+t_idx-1, ndxET, m, thisD) = KLIC_tmp;
    
            end



        end
    end
end

max_large_weight_sum_all = max(large_weight_sum_all,[],'all');


%% plot comparison across ET modes
ETlabel = {'binsOnly', 'binsAndMeans'};
for thisD = 1 : length(DATALABELS)
    datalabel = DATALABELS{thisD};
    for m = 1 : length(models)
        
        modellabel  = models{m};

        theseRELESS = RELESS(:,:,m,thisD);
        theseESS    = ESS(:,:,m,thisD);
        firstT      = find(any(~isnan(theseRELESS), 2), 1, 'first');
        ETdates     = dates(firstT:end);

        % relative ESS
        thisfig     = figure;
        hold on
        plot(ETdates, theseRELESS(firstT:end,1), '-', 'linewidth', 2);
        plot(ETdates, theseRELESS(firstT:end,2), '-.', 'linewidth', 3);
        legend(ETlabel)
        ylim([0 1])
        xtickdates(ETdates)
        title(sprintf('%s \n(%s)', datalabel, modellabel))
        wrapthisfigure(thisfig, sprintf('RELESS-%s-%s', datalabel, modellabel), wrap);
        
        % ESS
        thisfig     = figure;
        hold on
        plot(ETdates, theseESS(firstT:end,1), '-', 'linewidth', 2);
        plot(ETdates, theseESS(firstT:end,2), '-.', 'linewidth', 3);
        %yline(ESS_threshold,':k');
        small_ESS_1 = theseESS(firstT:end,1)<ESS_threshold;
        small_ESS_2 = theseESS(firstT:end,2)<ESS_threshold;
        if any(small_ESS_1 ~= 0)
            xline(ETdates(small_ESS_1), '-', 'linewidth', 2);
        end
        if any(small_ESS_2 ~= 0)
            xline(ETdates(small_ESS_2), '-.', 'linewidth', 3);
        end
        legend(ETlabel)
        xtickdates(ETdates)
        title(sprintf('%s \n(%s)', datalabel, modellabel))
        wrapthisfigure(thisfig, sprintf('ESS-%s-%s', datalabel, modellabel), wrap);

        % omega
        theseomega_x = omega_x(:,:,:,m,thisD);

        for twist_x_idx = 1 : length(twist_x_vec)

            thisfig     = figure;
            hold on
            plot(ETdates, squeeze(theseomega_x(twist_x_idx,1,firstT:end)), '-', 'linewidth', 2);
            plot(ETdates, squeeze(theseomega_x(twist_x_idx,2,firstT:end)), ':', 'linewidth', 2);
            %legend({strcat(ETlabel{1},' \omega_{',num2str(twist_x_vec(twist_x_idx)),'}'),...
            %        strcat(ETlabel{2},' \omega_{',num2str(twist_x_vec(twist_x_idx)),'}')})
            legend(ETlabel)
            xtickdates(ETdates)
            title(sprintf('%s \n(%s) \\omega_{%d}', datalabel, modellabel,twist_x_vec(twist_x_idx)))
            wrapthisfigure(thisfig, sprintf('Omega-%d-%s-%s', twist_x_vec(twist_x_idx), datalabel, modellabel), wrap);

        end


        % large weights
        thisfig     = figure;
        sgtitle(sprintf('%s \n(%s)', datalabel, modellabel))
        subplot(1,2,1)
        hold on
        plot(ETdates, large_weight_counter_all(firstT:end,1,m,thisD), '-s', 'linewidth', 1);
        plot(ETdates, large_weight_counter_all(firstT:end,2,m,thisD), '-.s', 'linewidth', 1.5);
        xtickdates(ETdates)
        title('# of large weights');
        hold off
        subplot(1,2,2)
        hold on
        plot(ETdates, large_weight_sum_all(firstT:end,1,m,thisD), '-d', 'linewidth', 1);
        plot(ETdates, large_weight_sum_all(firstT:end,2,m,thisD), '-.d', 'linewidth', 1.5);
        title('sum of large weights');
        hold off
        legend(ETlabel)
        ylim([0 0.2])
        xtickdates(ETdates)
        wrapthisfigure(thisfig, sprintf('largeweight-%s-%s', datalabel, modellabel), wrap);

        % zero weights
        thisfig     = figure;
        hold on
        title(sprintf('%s \n(%s)', datalabel, modellabel))
        plot(ETdates, zero_weight_all(firstT:end,1,m,thisD), '-s', 'linewidth', 1);
        plot(ETdates, zero_weight_all(firstT:end,2,m,thisD), '-.s', 'linewidth', 1.5);
        xtickdates(ETdates)
        legend(ETlabel)
        wrapthisfigure(thisfig, sprintf('zeroweight-%s-%s', datalabel, modellabel), wrap);

        % small weights
        thisfig     = figure;
        hold on
        title(sprintf('%s \n(%s)', datalabel, modellabel))
        plot(ETdates, small_weight_counter_all(firstT:end,1,m,thisD), '-s', 'linewidth', 1);
        plot(ETdates, small_weight_counter_all(firstT:end,2,m,thisD), '-.s', 'linewidth', 1.5);
        xtickdates(ETdates)
        legend(ETlabel)
        wrapthisfigure(thisfig, sprintf('smallweight-%s-%s', datalabel, modellabel), wrap);


    end
end

%% finish
dockAllFigures
finishwrap
finishscript