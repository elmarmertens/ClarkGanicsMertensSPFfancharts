%% plots term structure of uncertainty

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*UNRCH> 

%% clear workspace
clear variables

close all
fclose all;
clc

%% some parameters

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};
MODELTYPES = {'trendgapSV', 'trendgapSVnoise2', 'const'};
MODELPRETTY = {'SV', 'SV (w/outliers)', 'CONST'};

Ndraws    = 3e3;

resultdir = localresultsMCMC;

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx90     = [2 8];
fontsize  = 18;

doTitle = false;
doYBAR  = false;

%% prepare latexwrapper
wrap = [];
initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% loop over variables
for m = 1 % [1 3] % : 2 % length(MODELTYPES)

    close all

    modeltype   = MODELTYPES{m};
    modelpretty = MODELPRETTY{m};

    for d =  1 : length(DATALABELS)
        %% prepare things


        datalabel = DATALABELS{d};


        modellabel = strcat(datalabel, 'STATE', modeltype);

        %% load data
        fprintf('Processing %s ... \n', modellabel)
        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
        mat = matfile(fullfile(resultdir, matfilename));


        Yquantiles    = mat.fcstYquantiles;
        Nz            = mat.Nz;

        YBARquantiles     = mat.fcstYBARquantiles;
        YBARskew          = mat.fcstYBARskew;
        Nbar              = mat.Nbar;

        YCURRENTquantiles = mat.fcstYCURRENTquantiles;
        YBARquantiles     = cat(2, YCURRENTquantiles, YBARquantiles);
        YCURRENTskew      = mat.fcstYCURRENTskew;
        YBARskew          = cat(2, YCURRENTskew, YBARskew);
        Nbar              = Nbar + 1;

        Ylabel            = mat.Ylabel;
        YBARlabel         = cat(1, 'year=0', mat.YBARlabel);


        Nhorizons     = mat.Nhorizons;
        dates         = mat.dates;
        T             = mat.T;
        datesQ        = mat.datesQ;
        doNIPA        = mat.doNIPA;

        Tstart        = mat.Tstart;
        horizons      = 0:Nhorizons-1;


        dateTicks     = datenum(1985:5:2022,1,1);

        %% compute term structures of uncertainty

        Yuncertainty    = range(Yquantiles(:,:,ndx68), 3);
        YBARuncertainty = range(YBARquantiles(:,:,ndx68), 3);

        %% plot selected
        HorizonLabels = arrayfun(@(x) sprintf('h=%d', x), horizons, 'uniformoutput', false);

        ndxHnear = [1 4];
        ndxH     = unique([8 : 4 : Nhorizons, Nhorizons]); % note: unique also does an automatic sort

        % full sample
        datendx = Tstart : T;
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hanniNear = plot(dates(datendx), Yuncertainty(datendx,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(datendx), Yuncertainty(datendx,ndxH), 'linewidth', 2);
        xticks(dateTicks)
        xtickdates(dates(datendx), 'keepticks')
        if doTitle
            title(datalabel)
        end
        wrapthisfigure(thisfig, sprintf('Yuncertainty-%s', modellabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('Yuncertainty-%s-WITHLEGEND', modellabel), wrap)


        %% YBAR
        if doYBAR

            HorizonLabels = arrayfun(@(x) sprintf('y=%d', x), 0:Nbar-1, 'uniformoutput', false); 

            % full sample
            datendx = Tstart : T;
            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)
            hanni = plot(dates(datendx), YBARuncertainty(datendx,:), 'linewidth', 2);
            xtickdates(dates(datendx))
            if doTitle
                title(datalabel)
            end
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s', modellabel), wrap, [], [], [], [], true);
            legend(hanni, HorizonLabels, 'Location','best')
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s-WITHLEGEND', modellabel), wrap)

            % pre covid
            T2019 = find(dates == datenum(2019,10,1));
            datendx = Tstart : T2019;
            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)
            hanni = plot(dates(datendx), YBARuncertainty(datendx,:), 'linewidth', 2);
            xtickdates(dates(datendx))
            if doTitle
                title(datalabel)
            end
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s-preCovid', modellabel), wrap, [], [], [], [], true);
            legend(hanni, HorizonLabels, 'Location','best')
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s-WITHLEGEND-preCovid', modellabel), wrap)

            % since covid
            datendx = T2019 - 7 : T;
            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)
            hanni = plot(dates(datendx), YBARuncertainty(datendx,:), 'linewidth', 2);
            xtickdates(dates(datendx))
            if doTitle
                title(datalabel)
            end
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s-sinceCovid', modellabel), wrap, [], [], [], [], true);
            legend(hanni, HorizonLabels, 'Location','best')
            wrapthisfigure(thisfig, sprintf('YBARuncertainty-%s-WITHLEGEND-sinceCovid', modellabel), wrap)
        end

    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
