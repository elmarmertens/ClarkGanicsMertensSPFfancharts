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
models            = {'STATEtrendgapSV', 'STATEconst' }; % ,  'STATEscaleSV'}; 'STATEtrendgapSVnoise2'
DATALABELS        = {'RGDP', 'PGDP', 'UNRATE'};
ET_types = {'binsOnly', 'binsAndMeans',...
            'NMeansOnly', 'NMeansAndVars',...   % normal means only, normal means and variances
            'GBMeansOnly', 'GBMeansAndVars',... % generalized beta means only, generalized beta means and variances
            'GBMeansVarsAndSkew',...            % generalized beta means, variances and skewness
            'GBVarsSkewAndSPFAnnPoints',...     % generalized beta variance and skewness while keeping the SPF annual point forecast as means             
            'SVBins',...                        % tilting CONST to SV probabilities assigned to SPF bins
            'binsOnlyNextYearOnly',...          % same as binsOnly, but only uses next-year SPF bins
            'GBMeansVarsAndSkewNextYearOnly'};  % same as GBMeansVarsAndSkew, but only uses next-year fitted GB moments

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

%% start loop over models
for m = 1 : length(models)

    modellabel = models{m};

    %% start loop over variables
    for thisD = 1  : length(DATALABELS)

        datalabel = DATALABELS{thisD};
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

        % only evaluate ET starting with 1992Q1 / 2009Q2 for UNRATE
  
        if strcmpi(datalabel,'UNRATE') ~= 1
            firstT = find(dates == datenum(1992,01,01),1,'first');
        else
            firstT = find(dates == datenum(2009,04,01),1,'first');
        end

        %firstT = max(find(~isnan(SPFprob_act.probs(:,1)),1,'first'),CGMmcmc.Tstart);
        lastT  = min(find(~isnan(SPFprob_act.probs(:,1)),1,'last'), CGMmcmc.T);
        numT   = lastT - firstT + 1;
        first_last_Date = dates([firstT;lastT]);

        first_last_DateLabel = datestr(first_last_Date, 'yyyyqq');

        % todo: maybe store firstT lastT in ET_results
        %         if ~isequal(firstT, ET_results.firstT)
        %             error('firstT mismatch')
        %         end
        %         if ~isequal(lastT, ET_results.lastT)
        %             error('lastT mismatch')
        %         end

        % load YdensityDraws
        YdensityDraws = permute(CGMmcmc.YdensityDraws, [2 3 1]);
        YdensityDraws = YdensityDraws(:,:,firstT:lastT);
        % load YBARdensityDraws
        YBARdensityDraws = permute(CGMmcmc.YBARdensityDraws, [2 3 1]);
        YBARdensityDraws = YBARdensityDraws(:,:,firstT:lastT);
        % load YCURRENTdensityDraws
        YCURRENTdensityDraws = permute(CGMmcmc.YCURRENTdensityDraws, [2 3 1]);
        YCURRENTdensityDraws = YCURRENTdensityDraws(:,:,firstT:lastT);

        %% see if YAVGdraws exist
        if any(any(cellfun(@(x) strcmp(x, 'YAVGdensityDraws'), who(CGMmcmc))))
            YAVGdensityDraws = permute(CGMmcmc.YAVGdensityDraws, [2 3 1]);
            YAVGdensityDraws = YAVGdensityDraws(:,:,firstT:lastT);
        else
            YdensityDraws    = permute(YdensityDraws, [2 3 1]);
            YAVGdensityDraws = NaN(size(YdensityDraws));
            for h = Kavg : size(YAVGdensityDraws, 3)
                YAVGdensityDraws(:,:,h) = sum(YdensityDraws(:,:,h+(-(Kavg-1):0)), 3) / Kavg;
            end
            YAVGdensityDraws = permute(YAVGdensityDraws, [3 1 2]);
            YdensityDraws    = permute(YdensityDraws, [3 1 2]);
        end
        if any(any(cellfun(@(x) strcmp(x, 'YAVGfuture'), who(CGMmcmc))))
            YAVGfuture       = CGMmcmc.YAVGfuture;
        else
            YAVGfuture = NaN(size(Yfuture));
            for h = Kavg : size(YAVGfuture, 2)
                YAVGfuture(:,h) = sum(Yfuture(:,h+(-(Kavg-1):0)), 2) / Kavg;
            end
        end

        %% clear CGMmcmc
        clear CGMmcmc

        %% loop over ET modes
        if strcmpi(datalabel,'RGDP')
            ET_types_idx_vec = 1:length(ET_types);
        else
            ET_types_idx_vec = 1:length(ET_types)-2;
        end

        for thisET = ET_types_idx_vec

            ETlabel = ET_types{thisET};

            if strcmpi(ETlabel,'SVBins') && strcmpi(modellabel,'STATEtrendgapSV')
                fprintf('skipping %s for %s and %s \n',datalabel,modellabel,ETlabel);
                continue
            end

            diaryfilename=fullfile(localstore, sprintf('kensington-EVAL-ET%s-Ndraws%d-%s%s.log', ETlabel, MCMCdraws, datalabel, modellabel));
            if exist(diaryfilename, 'file')
                delete(diaryfilename);
            end
            diary(diaryfilename);
            
            fprintf('now doing %s for %s and %s \n',datalabel,modellabel,ETlabel);

            %% load ET results
            ET_results = matfile(fullfile(ET_results_folder,sprintf('kensington-ET%s-Ndraws%d-%s%s.mat', ETlabel, MCMCdraws, datalabel, modellabel)));
            ET_w_all   = permute(ET_results.ET_w_all, [2 1]);
            ET_w_all   = ET_w_all(:,firstT:lastT);

            % check if ET subsample is the same
            if first_last_Date ~= ET_results.first_last_Date 
                error('ET dates not matching, check if matching CGM and ET file loaded!')
            end


            %% allocate memory

            % in keeping with the CMM template, allocate memory for T jumpoffs,
            % even though we have only T-Tstart+1 jumpoffs

            fcstYhat              = NaN(T,Nhorizons); % predictive mean
            fcstYcrps             = NaN(T,Nhorizons);
            fcstYquantiles        = NaN(T,Nhorizons,Nquantiles);
            fcstYpercentiles      = NaN(T,Nhorizons,Npercentiles);
            fcstYvol              = NaN(T,Nhorizons);
            fcstYlogscore         = NaN(T,Nhorizons);

            fcstYAVGhat              = NaN(T,Nhorizons); % predictive mean
            fcstYAVGcrps             = NaN(T,Nhorizons);
            fcstYAVGquantiles        = NaN(T,Nhorizons,Nquantiles);
            fcstYAVGpercentiles      = NaN(T,Nhorizons,Npercentiles);
            fcstYAVGvol              = NaN(T,Nhorizons);
            fcstYAVGlogscore         = NaN(T,Nhorizons);

            fcstYBARhat           = NaN(T,Nbar);
            fcstYBARvol           = NaN(T,Nbar);
            fcstYBARquantiles     = NaN(T,Nbar,Nquantiles);
            fcstYBARpercentiles   = NaN(T,Nbar,Npercentiles);

            fcstYCURRENThat          = NaN(T,1);
            fcstYCURRENTvol          = NaN(T,1);
            fcstYCURRENTquantiles    = NaN(T,1,Nquantiles);
            fcstYCURRENTpercentiles  = NaN(T,1,Npercentiles);


            %% compute RECESS for RGDP
            if strcmpi(datalabel, 'RGDP')
                theseDraws  = log(1 + YdensityDraws(1:5,:,:) ./ 100) .* 100;
                theseDraws  = permute(theseDraws, [3 1 2]);
                fcstYrecess = NaN(T, 5);
                fcstYrecess(firstT:lastT,:) = sum((theseDraws < 0) .* permute(ET_w_all, [2 3 1]), 3) .* 100;
                clear theseDraws
            end

            %% start loop over jump-offs
            for thisT=firstT:lastT

                thisToffset = thisT - firstT + 1;

                thisDate      = dates(thisT);
                thisDateLabel = datestr(thisDate, 'yyyyqq');
                fprintf('********************************** \n');
                fprintf('now doing %s (%s) in %s for %s \n',datalabel,ETlabel,thisDateLabel,modellabel);

                weights = ET_w_all(:,thisToffset);

                theseYdensityDraws        = transpose(YdensityDraws(:,:,thisToffset));
                theseYAVGdensityDraws     = transpose(YAVGdensityDraws(:,:,thisToffset));
                theseYBARdensityDraws     = transpose(YBARdensityDraws(:,:,thisToffset));
                theseYCURRENTdensityDraws = transpose(YCURRENTdensityDraws(:,:,thisToffset));

                % predictive densities
                thisYhat                  = sum(theseYdensityDraws .* weights, 1);
                fcstYhat(thisT,:)         = thisYhat;
                fcstYvol(thisT,:)         = sqrt(sum((theseYdensityDraws - thisYhat).^2 .* weights, 1));

                thisYhat                  = sum(theseYAVGdensityDraws .* weights, 1);
                fcstYAVGhat(thisT,:)      = thisYhat;
                fcstYAVGvol(thisT,:)      = sqrt(sum((theseYAVGdensityDraws - thisYhat).^2 .* weights, 1));

                thisYhat                  = sum(theseYBARdensityDraws .* weights, 1);
                fcstYBARhat(thisT,:)      = thisYhat;
                fcstYBARvol(thisT,:)      = sqrt(sum((theseYBARdensityDraws - thisYhat).^2 .* weights, 1));

                thisYhat                  = sum(theseYCURRENTdensityDraws .* weights, 1);
                fcstYCURRENThat(thisT,:)  = thisYhat;
                fcstYCURRENTvol(thisT,:)  = sqrt(sum((theseYCURRENTdensityDraws - thisYhat).^2 .* weights, 1));


                %% sort draws for quantiles and CRPS
                [sortedY, ndx]        = sort(theseYdensityDraws, 1);
                sortedYweights        = weights(ndx);

                [sortedYAVG, ndx]     = sort(theseYAVGdensityDraws, 1);
                sortedYAVGweights     = weights(ndx);

                [sortedYBAR, ndx]     = sort(theseYBARdensityDraws, 1);
                sortedYBARweights     = weights(ndx);
                
                [sortedYCURRENT, ndx] = sort(theseYCURRENTdensityDraws, 1);
                sortedYCURRENTweights = weights(ndx);

                %% quantiles 
                fcstYquantiles(thisT,:,:) = transpose(wprctile1(sortedY, quantileP, sortedYweights, true));
                fcstYpercentiles(thisT,:,:) = transpose(wprctile1(sortedY, percentiles, sortedYweights, true));

                fcstYAVGquantiles(thisT,:,:)   = transpose(wprctile1(sortedYAVG, quantileP, sortedYAVGweights, true));
                fcstYAVGpercentiles(thisT,:,:) = transpose(wprctile1(sortedYAVG, percentiles, sortedYAVGweights, true));

                fcstYBARquantiles(thisT,:,:)   = transpose(wprctile1(sortedYBAR, quantileP, sortedYBARweights));
                fcstYBARpercentiles(thisT,:,:) = transpose(wprctile1(sortedYBAR, percentiles, sortedYBARweights));

                if all(~isnan(theseYCURRENTdensityDraws)) % NaN can happen for certain NIPA jump off due to NaN in past data
                    fcstYCURRENTquantiles(thisT,:,:)   = transpose(wprctile1(sortedYCURRENT, quantileP, sortedYCURRENTweights));
                    fcstYCURRENTpercentiles(thisT,:,:) = transpose(wprctile1(sortedYCURRENT, percentiles, sortedYCURRENTweights));
                end
                %% CRPS and logscore

                for hh = 1 : Nhorizons
                    fcstYcrps(thisT,hh)      = crpsDrawsWeighted(Yfuture(thisT,hh), sortedY(:,hh), sortedYweights(:,hh), true);
                    fcstYlogscore(thisT,hh) = log(ksdensity(theseYdensityDraws(:,hh),Yfuture(thisT,hh),'Weights',weights));
                end
                for hh = 1 : Nhorizons
                    fcstYAVGcrps(thisT,hh)      = crpsDrawsWeighted(YAVGfuture(thisT,hh), sortedYAVG(:,hh), sortedYAVGweights(:,hh), true);
                    if ~all(isnan(theseYAVGdensityDraws(:,hh)))
                        fcstYAVGlogscore(thisT,hh) = log(ksdensity(theseYAVGdensityDraws(:,hh), YAVGfuture(thisT,hh),'Weights',weights));
                    end
                end


                fprintf('DONE WITH %s (%s) in %s for %s \n',datalabel,ETlabel,thisDateLabel,modellabel);
                fprintf('********************************** \n');

            end % thisT

            %% collect median, uncertainty etc
            fcstYmedian      = fcstYquantiles(:,:,ndxQ2);
            fcstYmederror    = fcstYmedian - Yfuture;
            fcstYhaterror    = fcstYhat    - Yfuture;

            fcstYAVGmedian      = fcstYAVGquantiles(:,:,ndxQ2);
            fcstYAVGmederror    = fcstYAVGmedian - YAVGfuture;
            fcstYAVGhaterror    = fcstYAVGhat    - YAVGfuture;

            % TODO: construct YBARfuture
            fcstYBARmedian   = fcstYBARquantiles(:,:,ndxQ2);
            %             fcstYBARmederror    = fcstYBARmedian - YBARfuture;
            %             fcstYBARhaterror    = fcstYBARhat    - YBARfuture;
            fcstYCURRENTmedian   = fcstYCURRENTquantiles(:,:,ndxQ2);
            %             fcstYCURRENTmederror    = fcstYCURRENTmedian - YCURRENTfuture;
            %             fcstYCURRENThaterror    = fcstYCURRENThat    - YCURRENTfuture;

            fcstYuncertainty      = range(fcstYquantiles(:,:,ndx68), 3);
            fcstYAVGuncertainty   = range(fcstYAVGquantiles(:,:,ndx68), 3);
            YBARuncertainty       = range(fcstYBARquantiles(:,:,ndx68), 3);
            YCURRENTuncertainty   = range(fcstYCURRENTquantiles(:,:,ndx68), 3);

            % Bowley Skew
            fcstYskew        = (fcstYquantiles(:,:,ndxQ3) - 2 .* fcstYquantiles(:,:,ndxQ2) + fcstYquantiles(:,:,ndxQ1)) ...
                ./ (fcstYquantiles(:,:,ndxQ3) - fcstYquantiles(:,:,ndxQ1));
            fcstYAVGskew        = (fcstYAVGquantiles(:,:,ndxQ3) - 2 .* fcstYAVGquantiles(:,:,ndxQ2) + fcstYAVGquantiles(:,:,ndxQ1)) ...
                ./ (fcstYAVGquantiles(:,:,ndxQ3) - fcstYAVGquantiles(:,:,ndxQ1));
            YBARskew = (fcstYBARquantiles(:,:,ndxQ3) - 2 .* fcstYBARquantiles(:,:,ndxQ2) + fcstYBARquantiles(:,:,ndxQ1)) ...
                ./ (fcstYBARquantiles(:,:,ndxQ3) - fcstYBARquantiles(:,:,ndxQ1));
            YCURRENTskew = (fcstYCURRENTquantiles(:,:,ndxQ3) - 2 .* fcstYCURRENTquantiles(:,:,ndxQ2) + fcstYCURRENTquantiles(:,:,ndxQ1)) ...
                ./ (fcstYCURRENTquantiles(:,:,ndxQ3) - fcstYCURRENTquantiles(:,:,ndxQ1));


            %% here we save the fc results after ET
            save(fullfile(localstore, sprintf('kensington-EVAL-ET%s-Ndraws%d-%s%s.mat', ETlabel, MCMCdraws, datalabel, modellabel)), ...
                '-regexp','fcst\w*', 'YAVGfuture', 'firstT','lastT','dates','Zdata','Yfuture','Nhorizons', ...
                'quantileP', 'percentiles', ...
                'Znanny', 'Ylabel', ...
                '-v7.3');
            fprintf('DONE %s for %s and %s \n',datalabel,modellabel,ETlabel);

            diary off;
        end
    end
end
