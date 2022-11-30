%% plotting cumulative histograms and fitted normal and generalized beta CDFs, along with moments

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
addpath ../kensingtonToolbox

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*ASGLU>
%#ok<*NASGU>
%#ok<*UNRCH>
%#ok<*SAGROW>

%% clear workspace
clear variables

close all
fclose all;
clc


%% parameters
DATALABELS = {'RGDP', 'PGDP', 'UNRATE'};
Nvariables = length(DATALABELS);
fontsize   = 18;
doBinEdgeTicks = true;

datadir     = fullfile('..', 'kensingtonDataMatfiles');
resultdir   = localresultsMCMC;
modellabel  = 'STATEtrendgapSV';
Ndraws      = 3000;

%% line colors
SPFcolor   = [0 0.4470 0.7410]; % dark blue
Ncolor    = Colors4Plots(7); % dark red
Gbetacolor = Colors4Plots(3); % dark yellow
Modelcolor = 'k';
RealizedColor = Colors4Plots(5);

%% latexwrapper
doPlot = false;
wrap   = [];
if doPlot
    titlename = sprintf('figuresFitCDFSPFplots');
else
    titlename = sprintf('figuresFitNOCDFSPFplots');
end
initwrap

%% load fitted values
FitCDFSPF_all = load(fullfile(datadir,'kensingtonFitDistSPF.mat'));



%% loop over variables
for d = 2 %: length(DATALABELS)

thislabel = DATALABELS{d};

thisDATA = FitCDFSPF_all.(thislabel);

SPFdata        = matfile(fullfile(datadir,sprintf('kensington%sdata',thislabel)));
YBARfuture     = SPFdata.YBARfuture;
YCURRENTfuture = SPFdata.YCURRENTfuture;
thisRealiz     = [YCURRENTfuture,YBARfuture];
Nbar           = SPFdata.Nbar;
SPFdates       = SPFdata.dates;
Zdata          = SPFdata.Zdata;
Ybarforecast   = SPFdata.YactualForecasts(:,end-Nbar+1:end); % forecasts for calendar years starting with next year


% check if dates are the same in fitting results and data file
if ~isequal(SPFdates,thisDATA.dates)
    warning('dates do NOT match between fitting results and master data!');
end

matfilename = sprintf('slimCGMmcmc-%s%s-Ndraws%d',thislabel, modellabel, Ndraws);
CGMmat = matfile(fullfile(resultdir, matfilename));

% check if dates are the same in fitting results and MCMC results file
if ~isequal(SPFdates,CGMmat.dates)
    warning('dates do NOT match between fitting results and MCMC results!');
end

% only do plotting starting with 1992Q1 / 2009Q2 for UNRATE
if strcmp(thislabel,'UNRATE') ~= 1
    thisfirstT = find(thisDATA.dates == datenum(1992,01,01),1,'first');
else
    thisfirstT = find(thisDATA.dates == datenum(2009,04,01),1,'first');
end
thislastT =  find(~isnan(thisDATA.bins(:,1)),1,'last');

legtext_resnorm = {'Current year','Next year','2 years ahead','3 years ahead'};

marker_vec = {'none','s','d','o'};
linestylevec = {'-','--',':','-.'};

this_betameanstdskew = NaN(thislastT,4,3); % time, hz, [mean, variance, 3rd moment skewness]
this_normalmeanstd   = NaN(thislastT,4,2); % time, hz, [mean, variance]

betaSPFdrps = NaN(thislastT,1+Nbar);
normSPFdrps = NaN(thislastT,1+Nbar);

Ycurrentforecast = NaN(size(Ybarforecast,1),1);
thisForecast   = [Ycurrentforecast, Ybarforecast];
thisVol        = [CGMmat.fcstYCURRENTvol, CGMmat.fcstYBARvol];
thisSkew       = [CGMmat.fcstYCURRENTskew, CGMmat.fcstYBARskew];

SPFHISTdata        = matfile(fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS',thislabel)));


%% loop over jump offs

    for thisT = thisfirstT : thislastT

        Nhz             = thisDATA.bins(thisT,1);
        allSPFbinEdges  = thisDATA.bins(thisT,2:find(~isnan(thisDATA.bins(thisT,1:end)),1,'last'));
        SPFbinEdges     = allSPFbinEdges(1:length(allSPFbinEdges)/Nhz);

        SPFcdf           = reshape(thisDATA.cdf(thisT,1:Nhz*(length(SPFbinEdges)+1)),[],Nhz) * 100; % DON'T drop current year

        for hh = 1 : Nhz


            thisBetaparams      = vec(thisDATA.params_B(thisT,hh,:));
            thisNormalparams    = vec(thisDATA.params_N(thisT,hh,:));

            betamean_tmp        = thisBetaparams(1) + (thisBetaparams(2) - thisBetaparams(1)) * (thisBetaparams(3) / (thisBetaparams(3) + thisBetaparams(4)));
            betastd_tmp         = sqrt(thisBetaparams(3) * thisBetaparams(4) * (thisBetaparams(2) - thisBetaparams(1))^2 / ((thisBetaparams(3) + thisBetaparams(4))^2 * (thisBetaparams(3) + thisBetaparams(4) + 1)));
            
            betaskew_tmp        = 2 * (thisBetaparams(4) - thisBetaparams(3)) * sqrt(thisBetaparams(3) + thisBetaparams(4) + 1) / ((thisBetaparams(3) + thisBetaparams(4) + 2) * sqrt(thisBetaparams(3) * thisBetaparams(4)));

            this_betameanstdskew(thisT,hh,1:3) = [betamean_tmp;betastd_tmp;betaskew_tmp];
            this_normalmeanstd(thisT,hh,1:2)   = [thisNormalparams(1);thisNormalparams(2)];

            % calculate CDF, implied by fitted distribution
            betacdf_tmp     = gen_beta_cdf(thisBetaparams,SPFbinEdges')';
            normcdf_tmp     = normcdf(SPFbinEdges',thisNormalparams(1),thisNormalparams(2))';
                if ~isnan(thisRealiz(thisT,hh))
                    thisYbarhit             = double(thisRealiz(thisT,hh) < SPFbinEdges);
                    betaSPFdrps(thisT,hh) = sum((betacdf_tmp - thisYbarhit).^2);
                    normSPFdrps(thisT,hh) = sum((normcdf_tmp - thisYbarhit).^2);
                end

            % SPF + N + GBeta in one
        if doPlot
            thisfig = figure;
            set(gca, 'FontSize', fontsize)

            hold on
            
            hSPF  = plot(SPFbinEdges, SPFcdf(1:end-1,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            hN    = plot(SPFbinEdges(1):0.01:SPFbinEdges(end), 100*normcdf((SPFbinEdges(1):0.01:SPFbinEdges(end))',thisDATA.params_N(thisT,hh,1),thisDATA.params_N(thisT,hh,2)), '-', 'color', Ncolor, 'linewidth', 1);
            hGbeta= plot(SPFbinEdges(1):0.01:SPFbinEdges(end), 100*gen_beta_cdf(thisBetaparams,(SPFbinEdges(1):0.01:SPFbinEdges(end))'), '-', 'color', Gbetacolor, 'linewidth', 1);
            
            yticks(0 : 25 : 100);
            ylim([0 100]);
            xlim([SPFbinEdges(1) SPFbinEdges(end)]);
            if doBinEdgeTicks
                xticks(SPFbinEdges);
            end
            grid on
            orient landscape
            wrapthisfigure(thisfig, sprintf('SPFFitcdfplot-%s-%d-%s', thislabel, hh, datestr(thisDATA.dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
            hl = legend([hSPF; hN; hGbeta], 'SPF' ,'Normal', 'Gen. beta', 'location', 'northwest');
            wrapthisfigure(thisfig, sprintf('SPFFitcdfplot-%s-%d-%s-WITHLEGEND', thislabel, hh, datestr(thisDATA.dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
            delete(hl)
            % add realized value
            hData = plot([thisRealiz(thisT,hh) thisRealiz(thisT,hh)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
            wrapthisfigure(thisfig, sprintf('SPFFitcdfplot-%s-%d-%s-WITHDATA', thislabel, hh, datestr(thisDATA.dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
            legend([hSPF; hN; hGbeta; hData], 'SPF' ,'Normal', 'Gen. beta', 'Realization', 'location', 'northwest')
            wrapthisfigure(thisfig, sprintf('SPFFitcdfplot-%s-%d-%s-WITHDATALEGEND', thislabel, hh, datestr(thisDATA.dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);

            hold off
            close(thisfig)
        end
        end % end hh (horizons)
    end % end jump offs
%%
    max_resnorm = max(vec([thisDATA.resnorm_N(thisfirstT:thislastT,2:end),thisDATA.resnorm_B(thisfirstT:thislastT,2:end)]));

    if doPlot
    thisfig = figure;
    hold on
    for hz = 1:4
        plot(thisDATA.dates,thisDATA.resnorm_N(:,hz),'Marker',marker_vec{hz});
    end

    ylim([0 max_resnorm]);

    %xline(thisDATA.dates(thisDATA.dates==datenum('1992Q1','yyyyQQ')),'r','LineWidth',1.5);
    
    xlim([thisDATA.dates(find(~isnan(thisDATA.resnorm_N(:,2)),1,"first")) thisDATA.dates(thislastT)]);
    datetick('x','yyyy','keeplimits');
  
    legend(legtext_resnorm{1:Nbar+1},'Box','on','Location','northwest');
    wrapthisfigure(thisfig, sprintf('SPFFitcdfresnormN-%s',thislabel), wrap, [], [], [], [], true);
    close(thisfig)

    thisfig = figure;
    hold on
    for hz = 1:Nbar+1
        plot(thisDATA.dates,thisDATA.resnorm_B(:,hz),'Marker',marker_vec{hz});
    end

    ylim([0 max_resnorm]);

    %xline(thisDATA.dates(thisDATA.dates==datenum('1992Q1','yyyyQQ')),'r','LineWidth',1.5);
    
    xlim([thisDATA.dates(find(~isnan(thisDATA.resnorm_B(:,2)),1,"first")) thisDATA.dates(thislastT)]);
    datetick('x','yyyy','keeplimits');
    legend(legtext_resnorm{1:Nbar+1},'Box','on','Location','northwest');
    wrapthisfigure(thisfig, sprintf('SPFFitcdfresnormGenbeta-%s',thislabel), wrap, [], [], [], [], true);
    close(thisfig)

    end
    
    for hh = 1:Nbar+1

        % plot means from N, GB and SPF point forecasts, then add realization
        thisfig = figure;
        hold on

        if hh > 1 % annual forecasts still pending for current year
            plot(thisDATA.dates,thisForecast(:,hh),'Color',SPFcolor,'LineWidth',1.5,'LineStyle',linestylevec{1});
        end
        plot(thisDATA.dates,this_normalmeanstd(:,hh,1),'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
        plot(thisDATA.dates,this_betameanstdskew(:,hh,1),'Color',Gbetacolor,'LineWidth',1.5,'LineStyle',linestylevec{3});
        
        xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        if hh > 1
            legend('SPF','Normal','Gen. beta','Box','on','Location','best');
        else
            legend('Normal','Gen. beta','Box','on','Location','best');
        end
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitmeanplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        
        % add realized value
        hData = plot(thisDATA.dates,thisRealiz(:,hh), '-.', 'color', RealizedColor, 'linewidth', 1.5);
        if hh > 1
            legend('SPF','Normal','Gen. beta','Realization','Box','on','Location','best');
        else
            legend('Normal','Gen. beta','Realization','Box','on','Location','best');
        end
        wrapthisfigure(thisfig, sprintf('SPFFitmeanplot-WITHDATALEGEND-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        close(thisfig)

        % plot volatilities from N, GB and SV model
        thisfig = figure;
        hold on

        plot(thisDATA.dates,thisVol(:,hh),'Color',SPFcolor,'LineWidth',1.5,'LineStyle',linestylevec{1});
        plot(thisDATA.dates,this_normalmeanstd(:,hh,2),'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
        plot(thisDATA.dates,this_betameanstdskew(:,hh,2),'Color',Gbetacolor,'LineWidth',1.5,'LineStyle',linestylevec{3});
        
        xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        legend('SV model','Normal','Gen. beta','Box','on','Location','best');
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitvolplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        close(thisfig)

        % plot standardized 3rd moment skewness from GB and SV model
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        
        %plot(thisDATA.dates,thisSkew(:,hh),'Color',SPFcolor,'LineWidth',1.5,'LineStyle',linestylevec{1});
        plot(thisDATA.dates,this_betameanstdskew(:,hh,3),'Color',Colors4Plots(8),'LineWidth',5,'LineStyle',linestylevec{1});
        yline(0,'Color',SPFcolor,'LineWidth',1,'LineStyle',linestylevec{3});

        xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        %ylim([-1 1]);
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitskewplot-%s-%d',thislabel,hh), wrap, [], [], [], [], false);
        close(thisfig)

        

        % plot DRPS, SPF only for h>1 for now
        
        thisfig = figure;
        hold on

        startDRPS_tmp = find(~isnan(normSPFdrps(:,hh)),1,'first');
        endDRPS_tmp   = find(~isnan(normSPFdrps(:,hh)),1,'last');
        numpadend     = thislastT - endDRPS_tmp;
        numpadstart   = startDRPS_tmp-1;
        recmean_fact  = (1 : endDRPS_tmp-startDRPS_tmp+1)';


        if hh > 1
            %lastNotNaN   = find(~isnan(SPFHISTdata.fcstSPFdrps(:,hh-1)),1,'last');
            %firstNotNaN  = find(~isnan(SPFHISTdata.fcstSPFdrps(:,hh-1)),1,'first');
            %numpadend    = thislastT-lastNotNaN;
            %recmean_fact = (1 : (lastNotNaN-thisfirstT+1))';
            recmeanSPF = [NaN(startDRPS_tmp-1,1);cumsum(SPFHISTdata.fcstSPFdrps(startDRPS_tmp:endDRPS_tmp,hh-1))./recmean_fact;NaN(numpadend,1)];
            plot(thisDATA.dates,recmeanSPF,'Color',SPFcolor,'LineWidth',1.5,'LineStyle',linestylevec{1});
        end
        %lastNotNaN   = find(~isnan(normSPFdrps(:,hh)),1,'last');
        %numpadend    = thislastT-lastNotNaN;
        %recmean_fact = [NaN(thisfirstT-1,1);(1 : (lastNotNaN-thisfirstT+1))'];
        recmeanNorm  = [NaN(startDRPS_tmp-1,1);cumsum(normSPFdrps(startDRPS_tmp:endDRPS_tmp,hh))./recmean_fact;NaN(numpadend,1)];
        recmeanGBeta = [NaN(startDRPS_tmp-1,1);cumsum(betaSPFdrps(startDRPS_tmp:endDRPS_tmp,hh))./recmean_fact;NaN(numpadend,1)];
        plot(thisDATA.dates,recmeanNorm,'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
        plot(thisDATA.dates,recmeanGBeta,'Color',Gbetacolor,'LineWidth',1.5,'LineStyle',linestylevec{3});
        
        xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        if hh > 1
            legend('SPF','Normal','Gen. beta','Box','on','Location','best');
        else
             legend('Normal','Gen. beta','Box','on','Location','best');
        end
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitDRPSplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        close(thisfig)

        % calculate relative DRPS, both differences and ratios

         thisfig = figure;
         hold on
         plot(thisDATA.dates,recmeanNorm - recmeanGBeta,'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
         xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        
        legend('Normal - Gen. beta','Box','on','Location','best');
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitDRPSreldiffNvsGBplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        close(thisfig)

         thisfig = figure;
         hold on
         plot(thisDATA.dates,recmeanNorm./recmeanGBeta,'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
         xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        
        legend('Normal/Gen. beta','Box','on','Location','best');
        datetick('x','yyyy','keeplimits');
        wrapthisfigure(thisfig, sprintf('SPFFitDRPSrelratioNvsGBplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
        close(thisfig)

        if hh > 1

            thisfig = figure;
            hold on
                   
          
            plot(thisDATA.dates,recmeanNorm-recmeanSPF,'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
            plot(thisDATA.dates,recmeanGBeta-recmeanSPF,'Color',Gbetacolor,'LineWidth',1.5,'LineStyle',linestylevec{3});
        
            xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        
            legend('Normal','Gen. beta','Box','on','Location','best');
            datetick('x','yyyy','keeplimits');
            wrapthisfigure(thisfig, sprintf('SPFFitDRPSreldiffvsSPFplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
            close(thisfig)

            thisfig = figure;
            hold on

            plot(thisDATA.dates,recmeanNorm./recmeanSPF,'Color',Ncolor,'LineWidth',1.5,'LineStyle',linestylevec{2});
            plot(thisDATA.dates,recmeanGBeta./recmeanSPF,'Color',Gbetacolor,'LineWidth',1.5,'LineStyle',linestylevec{3});
        
            xlim([thisDATA.dates(thisfirstT) thisDATA.dates(thislastT)]);
        
            legend('Normal','Gen. beta','Box','on','Location','best');
            datetick('x','yyyy','keeplimits');
            wrapthisfigure(thisfig, sprintf('SPFFitDRPSrelratiovsSPFplot-%s-%d',thislabel,hh), wrap, [], [], [], [], true);
            close(thisfig)


        end

    end % end horizons
end % end DATALABELS

%% finish / clean up
finishwrap
finishscript
dockAllFigures


