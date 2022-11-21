%% plot skew estimates from model densities

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>

%% clear workspace
clear variables

close all
fclose all;
clc

%% some parameters

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};
MODELTYPES = {'trendgapSV', 'const'};

% doETAVAR = false;

doShowCURRENT = true;

Ndraws    = 3e3;

resultdir = localresultsMCMC;

fontsize  = 18;

%% prepare latexwrapper
wrap = [];
% initwrap

if isempty(wrap) && ~isdesktop
    initwrap
end

%% loop over variables
for m = 1 % : length(MODELTYPES)
    modeltype = MODELTYPES{m};
    close all
    for d = 1 : length(DATALABELS)
        %% prepare things


        datalabel = DATALABELS{d};


        % construct a model-label indicating dataset and important parameters, to
        % append to picture names
        modellabel = datalabel;
        %#ok<*UNRCH>

        %         if doETAVAR
        %             modellabel = strcat(modellabel, 'VARSTATE', modeltype);
        %         else
        modellabel = strcat(modellabel, 'STATE', modeltype);
        %         end

        %% load data
        fprintf('Processing %s ... \n', modellabel)

        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
        mat = matfile(fullfile(resultdir, matfilename));


        %     YdensityDraws = mat.YdensityDraws;
        Ydraws        = mat.Ydraws;
        Yfuture       = mat.Yfuture;
        Nz            = mat.Nz;
        Zdata         = mat.Zdata;
        Znanny        = mat.Znanny;
        Nbar          = mat.Nbar;

        Yskew              = mat.fcstYskew;
        YBARskew           = mat.fcstYBARskew;

        if doShowCURRENT
            Nbar              = Nbar + 1;

            YCURRENTskew       = mat.fcstYCURRENTskew;
            YBARskew           = cat(2, YCURRENTskew, YBARskew);
        end

        Nhorizons     = mat.Nhorizons;
        dates         = mat.dates;
        T             = mat.T;
        datesQ        = mat.datesQ;
        doNIPA        = mat.doNIPA;

        Tstart        = mat.Tstart;
        horizons      = 0:Nhorizons-1;


        %% Y: plot selected 2D
        HorizonLabels = arrayfun(@(x) sprintf('h=%d', x), horizons, 'uniformoutput', false);

        ndxHnear = [1 2 5];
        ndxH     = unique([9 : 4 : Nhorizons, Nhorizons]); % note: unique also does an automatic sort

        % full sample
        datendx = Tstart : T;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanniNear = plot(dates(datendx), Yskew(datendx,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(datendx), Yskew(datendx,ndxH), 'linewidth', 2);
        xtickdates(dates(datendx))
        ylim([-1 1])
        title(datalabel)
        wrapthisfigure(thisfig, sprintf('skew-%s', modellabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('skew-%s-WITHLEGEND', modellabel), wrap)

        % pre covid
        T2019 = find(dates == datenum(2019,10,1));
        datendx = Tstart : T2019;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanniNear = plot(dates(datendx), Yskew(datendx,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(datendx), Yskew(datendx,ndxH), 'linewidth', 2);
        xtickdates(dates(datendx))
        title(datalabel)
        ylim([-1 1])
        wrapthisfigure(thisfig, sprintf('skew-%s-preCovid', modellabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('skew-%s-WITHLEGEND-preCovid', modellabel), wrap)

        % since covid
        datendx = T2019 - 7 : T;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanniNear = plot(dates(datendx), Yskew(datendx,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(datendx), Yskew(datendx,ndxH), 'linewidth', 2);
        xtickdates(dates(datendx))
        title(datalabel)
        ylim([-1 1])
        wrapthisfigure(thisfig, sprintf('skew-%s-sinceCovid', modellabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('skew-%s-WITHLEGEND-sinceCovid', modellabel), wrap)

        %% YBAR
        HorizonLabels = arrayfun(@(x) sprintf('y=%d', x), 0:Nbar-1, 'uniformoutput', false);

        % full sample
        datendx = Tstart : T;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanni = plot(dates(datendx), YBARskew(datendx,:), 'linewidth', 2);
        xtickdates(dates(datendx))
        title(datalabel)
        ylim([-1 1])
        wrapthisfigure(thisfig, sprintf('skewBAR-%s', modellabel), wrap, [], [], [], [], true);
        legend(hanni, HorizonLabels, 'Location','best')
        wrapthisfigure(thisfig, sprintf('skewYBAR-%s-WITHLEGEND', modellabel), wrap)

        % pre covid
        T2019 = find(dates == datenum(2019,10,1));
        datendx = Tstart : T2019;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanni = plot(dates(datendx), YBARskew(datendx,:), 'linewidth', 2);
        xtickdates(dates(datendx))
        title(datalabel)
        wrapthisfigure(thisfig, sprintf('skewBAR-%s-preCovid', modellabel), wrap, [], [], [], [], true);
        legend(hanni, HorizonLabels, 'Location','best')
        wrapthisfigure(thisfig, sprintf('skewYBAR-%s-WITHLEGEND-preCovid', modellabel), wrap)

        % since covid
        datendx = T2019 - 7 : T;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanni = plot(dates(datendx), YBARskew(datendx,:), 'linewidth', 2);
        xtickdates(dates(datendx))
        title(datalabel)
        ylim([-1 1])
        wrapthisfigure(thisfig, sprintf('skewBAR-%s-sinceCovid', modellabel), wrap, [], [], [], [], true);
        legend(hanni, HorizonLabels, 'Location','best')
        wrapthisfigure(thisfig, sprintf('skewYBAR-%s-WITHLEGEND-sinceCovid', modellabel), wrap)

    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
