%% compute and plot bias in SPF-consistent expectations of VAReta

% todo:
% - load predictive means directly (once stored with MCMC)
% - in MCMC: compute predictive means analytically (instead of MC)
% - compute bias directly in MCMC to generate also uncertainty around bias

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*ASGLU>
%#ok<*NASGU>
%#ok<*UNRCH> 

%% clear workspace
clear variables

close all
fclose all;
clc


%% some parameters


modeltype   = 'VARprior10STATEtrendgapSV';
modelpretty = 'non-MDS';
% modeltype   = 'VARprior5STATEtrendgapSV';
% modelpretty = 'non-MDS (tighter prior SV)';
% modeltype   = 'VARprior5STATEconst';
% modelpretty = 'non-MDS (tighter prior; CONST)';
% Ndraws      = 3e3;
% modeltype   = 'VARprior50STATEtrendgapSV';
% modelpretty = 'non-MDS (loose prior; SV)';
Ndraws      = 3e3;

resultdir = localresultsMCMC;
resultdir = '~/jam/lager/KENSINGTON/kensingtonresults3';

doMedian = false; % does not change results much

fontsize = 18;

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};

%% initwrap
wrap = [];
titlename = sprintf('SPFbias-%s', modeltype);
if doMedian
    titlename = strcat(titlename, '-Median');
end
initwrap

%% loop over datalabels

for d = 1 : length(DATALABELS)

    close all
    datalabel   = DATALABELS{d};
    modellabel  = strcat(datalabel, modeltype);

    %% prepare latexwrapper
    %     wrap = [];
    %     titlename = sprintf('SPFbias-%s-%s', datalabel, modeltype);
    %     initwrap

    fprintf('Processing %s ... \n', datalabel)

    %% load data
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    if doMedian
        Yhat  = mat.fcstYmedian;
    else
        Yhat  = mat.fcstYhatRB; 
    end
    Ydraws        = mat.Ydraws;

    Yfuture       = mat.Yfuture;
    Nz            = mat.Nz;
    Zdata         = mat.Zdata;
    Znanny        = mat.Znanny;

    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    Tstart        = mat.Tstart;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    horizonLabels = arrayfun(@(x) sprintf('h=%d', x), 0:Nhorizons-1, 'uniformoutput', false);

    T2019 = find(dates == datenum(2019,10,1));

    %% compute bias
    if doMedian
        spfYhat = median(Ydraws(2:end,:,:),3)';
    else
        spfYhat = mean(Ydraws(2:end,:,:),3)';
    end
    spfBias = Yhat - spfYhat;

    %% plot G-eigenvalues
    %     Gpriorlambda = mat.Gpriorlambda;
    %     Glambda      = mat.Glambda;
    %
    %     thisfig = figure;
    %     hold on
    %     set(gca, 'FontSize', fontsize)
    %     plot(dates, median(Gpriorlambda,1), 'k--', 'linewidth', 2)
    %     plot(dates, prctile(Gpriorlambda, [25 75], 1), 'k--')
    %
    %     plot(dates, median(Glambda,1), 'b-', 'linewidth', 2)
    %     plot(dates, prctile(Glambda, [25 75], 1), 'b-')
    %
    %
    %     ylim([min([ylim, 0]) max(ylim)])
    %     xtickdates(dates([Tstart T]))
    %     wrapthisfigure(thisfig, sprintf('Geigenvalues-%s', modellabel), wrap)

    for doPreCOVID = [true false]
        if doPreCOVID
            covidLabel = '-preCOVID';
            Tend       = T2019;
        else
            covidLabel = '';
            Tend       = T;
        end

        if doMedian
            covidLabel = strcat(covidLabel, '-Median');
        end

        %% plot bias; all in one, to collect YLIM
        %         thisfig = figure;
        %         h = plot(dates(Tstart:end) , spfBias(Tstart:end, :), 'linewidth', 2);
        %         legend(h, horizonLabels, 'Location', 'best')
        %         xtickdates(dates(Tstart:Tend))
        %         YLIM = ylim;
        %         close(thisfig)
        %         %     wrapthisfigure(thisfig, sprintf('SPFbiasAll-%s-%s-WITHLEGEND', ...
        %         %         datalabel, modeltype), wrap, [], [], [], [], true);
        %
        %         %% plot bias for groups of horizons
        %         hgroups = {[0 3 7 11 15], 0:3, 4:7, 8:11, 12:15}; % hardcoded for SPF
        %         for hh = 1 : length(hgroups)
        %             hndx = hgroups{hh} + 1;
        %             hndx = hndx(hndx <= Nhorizons);
        %             if ~isempty(hndx) % to catch cases where we have less than 16 Horizons
        %                 thisfig = figure;
        %                 hold on
        %                 h = plot(dates(Tstart:Tend) , spfBias(Tstart:Tend, hndx), 'linewidth', 2);
        %                 ylim(YLIM);
        %                 box off
        %                 set(gca, 'FontSize', fontsize)
        %                 nbershades(dates(Tstart:Tend))
        %                 xticks(datenum(1970:5:2020,1,1))
        %                 xtickdates(dates(Tstart:Tend), 'keepticks')
        %                 legend(h, horizonLabels(hndx), 'Location', 'best')
        %                 wrapthisfigure(thisfig, sprintf('SPFbias-group%d-%s-%s%s-WITHLEGEND', hh-1, ...
        %                     datalabel, modeltype, covidLabel), wrap, [], [], [], [], true);
        %                 % xticks(datenum(1970:10:2020,1,1))
        %                 %             xtickdates(dates(Tstart:T2019), 'keepticks')
        %                 %             wrapthisfigure(thisfig, sprintf('SPFbias-group%d-%s-%s-preCOVID-WITHLEGEND', hh-1, ...
        %                 %                 datalabel, modeltype), wrap, [], [], [], [], false);
        %             end
        %         end

        %% plot selected
        horizons      = 0:Nhorizons-1;
        HorizonLabels = arrayfun(@(x) sprintf('h=%d', x), horizons, 'uniformoutput', false);

        ndxHnear = [1 4];
        ndxH     = unique([8 : 4 : Nhorizons, Nhorizons]); % note: unique also does an automatic sort

        % full sample
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        box off
        hanniNear = plot(dates(Tstart:Tend), spfBias(Tstart:Tend,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(Tstart:Tend), spfBias(Tstart:Tend,ndxH), 'linewidth', 2);
        nbershades(dates(Tstart:Tend))
        xticks(datenum(1970:5:2020,1,1))
        xtickdates(dates(Tstart:Tend), 'keepticks')
        wrapthisfigure(thisfig, sprintf('SPFbias-%s%s', modellabel, covidLabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('SPFbias-%s%s-WITHLEGEND', modellabel, covidLabel), wrap)
    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
