%% SEP-style annual fan charts from SV model

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

DATALABELS = {'RGDP', 'UNRATE', 'CPI'};


doDateAxis = true;
doTitle    = false;

modeltype = 'trendgapSV';
modelpretty = 'SV';

Ndraws    = 3e3;

resultdir   = localresultsMCMC;

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
fontsize  = 18;


colorMCMC = Colors4Plots(1);
colorData = Colors4Plots(7);

%% loop over datalabels
for d = 1 : 3

    % tndx = 158 : 4 : 185;  % 2011 - 2014
    % tndx = 202 : 4 : 214;  % 2019 - 2022
    tndx = 158 : 4 : 214;    % 2011 - 2022

    tndx = 158 : 215;        % 2011 - end of sample (2022Q2)

    datalabel = DATALABELS{d};

    modellabel = strcat(datalabel, 'STATE', modeltype);

    %% load data

    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    Tstart        = mat.Tstart;

    %% pick outputs
    switch datalabel
        case {'RGDP', 'PGDGP', 'CPI'}
            Yhat          = mat.fcstYAVGhat;
            Yquantiles    = mat.fcstYAVGquantiles;
            Yfuture       = mat.YAVGfuture;
        case {'UNRATE'}
            Yhat          = mat.fcstYhat;
            Yquantiles    = mat.fcstYquantiles;
            Yfuture       = mat.Yfuture;
        otherwise
            error('datalabel <<%s>> not yet implemented')

    end


    %% pad extra dates at end
    xdates = genrQdates(year(dates(end)), year(dates(end))+4, 1);
    dates  = union(dates, xdates);

    %% prepare latexwrapper
    wrap = [];
    titlename = strcat('SEPfancharts-', datalabel, '-', modeltype);
    if isempty(wrap) && ~isdesktop
        initwrap
    end

    fprintf('Processing %s ... \n', modellabel)

    %% plot predictive density of Y

    if isempty(tndx)
        if ismac && isdesktop
            tndx = 166; % T - 12;
        else
            tndx = Tstart : T;
        end
    end

    for thisT = tndx

        showQuarters    = (3 : 4 : Nhorizons) - (datesQ(thisT) - 1);
        ndxShowQuarters = showQuarters + 1;

        if doDateAxis
            showQuarters    = dates(thisT+showQuarters);
        end

        mid     = Yhat(thisT,ndxShowQuarters)';
        tails68 = squeeze(Yquantiles(thisT,ndxShowQuarters,ndx68));

        yfuture = Yfuture(thisT,ndxShowQuarters)';

        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density
        hMCMC = plot(showQuarters,mid, '-', 'color', colorMCMC, 'LineWidth', 4);
        plot(showQuarters,tails68, '-', 'color', colorMCMC, 'LineWidth', 2)

        xticks(showQuarters)
        xlim(showQuarters([1 end]))
        if doDateAxis
            datetick('x', 'yyyyQQ', 'keeplimits', 'keepticks')
        end

        ylimits = ylim;
        switch upper(datalabel)
            case 'UNRATE'
                ylimits = [0 12];
            case 'RGDP'
                ylimits = [-2 8];
        end
        ylim(ylimits)

        grid on

        if doTitle
            title(sprintf('%s per %s', datalabel, datestr(dates(thisT), 'yyyyQQ')))
        end

        wrapthisfigure(thisfig, sprintf('Ysepfanchart-%s-%s', datestr(dates(thisT), 'yyyyqq'), modellabel), wrap, [], [], [], [], true);

        % add realized data
        if ~all(isnan(Yfuture(thisT,:)))
            ylim('auto')
            hFuture = plot(showQuarters, yfuture, 's:', 'color', colorData, 'linewidth', 3);
            ylimits = ylim;
            ylimits = ylim;
            switch upper(datalabel)
                case 'UNRATE'
                    ylimits = [0 12];
                case 'RGDP'
                    ylimits = [-2 8];
            end
            ylim(ylimits)
            wrapthisfigure(thisfig, sprintf('Ysepfanchart-%s-%s-WITHDATA', datestr(dates(thisT), 'yyyyqq'), modellabel), wrap, [], [], [], [], true);
            hl = legend([hMCMC hFuture], modelpretty, 'Realized values', 'location', 'best');
            wrapthisfigure(thisfig, sprintf('Ysepfanchart-%s-%s-WITHLEGENDDATA', datestr(dates(thisT), 'yyyyqq'), modellabel), wrap, [], [], [], [], false);
        end

    end

    %% finish datalabel loop
    finishwrap


end

%% finish / clean up
finishscript
dockAllFigures
