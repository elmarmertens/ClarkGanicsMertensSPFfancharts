%% STATE-SCALE-SV model
% uses MDS assumption

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
doTitle    = true;

doCharts   = false;

modeltype = 'trendgapSV';
modelpretty = 'SV';
% modeltype = 'const';
% % modeltype = 'trendgapCONST';
% modelpretty = 'CONST';
Ndraws    = 3e3;

resultdir   = localresultsMCMC;
resultdir = '~/jam/lager/KENSINGTON/kensingtonSPF2022Q2';

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
fontsize  = 16;


colorMCMC = Colors4Plots(1);
colorData = Colors4Plots(7);

%% import SEP table
SEPdata = importdata('../kensingtonDataUSSPF/SEPTable2Data.xlsx');

%% prepare latexwrapper
wrap = [];
titlename = strcat('SEPfanchartsUncertainty-', modeltype);
initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% loop over datalabels
for d = 1 : 3

    tndx = 158 : 4 : 210;  % 2008 - 2014
    %     tndx = 74 : 4 : 210;   % 1990 - 2014

    tndx = 158 : 215;        % 2011 - end of sample (2022Q2)

    datalabel = DATALABELS{d};

    modellabel = strcat(datalabel, 'STATE', modeltype);

    %% read SEP uncertainty data
    switch datalabel
        case 'RGDP'
            sheetname = 'RealGDP';
        case 'UNRATE'
            sheetname = 'UnemploymentRate';
        case 'CPI'
            sheetname = 'TotalConsumerPrices';
            warning('Matching CPI uncertainty with PCE from SEP')
        otherwise
            error('datalabel <<%s>> not yet implemented', datalabel)
    end

    if ispc
        sepdates = datenum(SEPdata.textdata.(sheetname)(2:end,1));
        sepRMSE  = SEPdata.data.(sheetname);
    else
        sepdates = x2mdate(SEPdata.data.(sheetname)(:,1));
        sepRMSE  = SEPdata.data.(sheetname)(:,2:end);
    end

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

    fprintf('Processing %s ... \n', modellabel)

    %% collect uncertainty data
    maxShowHorizons = Nhorizons / 4; % 4 or 3
    uncertaintyMCMC = NaN(length(tndx), maxShowHorizons);
    for tt = 1 : length(tndx)
        thisT = tndx(tt);
        showQuarters    = (3 : 4 : Nhorizons) - (datesQ(thisT) - 1);
        ndxShowQuarters = showQuarters + 1;
        uncertaintyMCMC(tt,:) =  range(Yquantiles(thisT,ndxShowQuarters,ndx68),3);
    end

    %% create labels
    % legendLabelsXXX    = {'this year', 'next year', 'two years ahead', 'three years ahead'};
    legendLabelsXXX    = {'this year', 'next year', 'two years', 'three years'};
    legendLabelsMCMC   = cellfun(@(x) sprintf('%s (%s)', x, modelpretty), legendLabelsXXX, 'UniformOutput', false);


    % match SEP dates
    sepndx          = ismember(sepdates, dates(tndx));
    uncertaintySEP  = 2.* sepRMSE(sepndx, :);
    legendLabelsSEP = cellfun(@(x) sprintf('%s (SEP)', x), legendLabelsXXX, 'UniformOutput', false);
    sepdates        = sepdates(sepndx);

    % drop three-year ahead
    %     legendLabelsMCMC = legendLabelsMCMC(1:3);
    %     legendLabelsSEP  = legendLabelsSEP(1:3);
    %     uncertaintyMCMC  = uncertaintyMCMC(:,1:3);
    %     uncertaintySEP   = uncertaintySEP(:,1:3);

    forecastOriginDates = dates(tndx);

    %% plot width of uncertainty bands -- SEP and MODELS
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hMCMC  = plot(forecastOriginDates,uncertaintyMCMC, '-', 'LineWidth', 2);
    % patch this year
    nanNdx = isnan(uncertaintyMCMC(:,1));
    plot(forecastOriginDates(~nanNdx),uncertaintyMCMC(~nanNdx), '-', 'LineWidth', 2, 'color', Colors4Plots(1));
    sepNaN = all(isnan(uncertaintySEP), 2);
    if doCharts
        sepLine = 'd';
    else
        sepLine = ':';
    end
    hSEP   = plot(sepdates(~sepNaN),uncertaintySEP(~sepNaN,1:maxShowHorizons), sepLine, 'LineWidth', 4);

    for nn = 1 : maxShowHorizons
        hMCMC(nn).Color = Colors4Plots(nn);
        hSEP(nn).Color  = Colors4Plots(nn);
    end

    switch datalabel
        case 'RGDP'
            ylim([0 12])
        case 'UNRATE'
            ylim([0 8])
        case 'CPI'
            ylim([0 7])
    end

    xticks(forecastOriginDates(1:4:end))
    xlim(forecastOriginDates([1 end]))
    if doDateAxis
        datetick('x', 'yyyy', 'keeplimits', 'keepticks')
    end

    grid on

    if doTitle
        title(sprintf('%s', datalabel))
    end

    legendLabels    = cat(1, legendLabelsMCMC, legendLabelsSEP);
    hanni           = transpose([hMCMC hSEP]);
    hl = legend(gca, hanni(:), legendLabels(:), 'location', 'southoutside', 'fontsize', 14, 'box', 'off');
    hl.NumColumns = maxShowHorizons;

    orient landscape
    wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-%s-WITHLEGEND', modellabel), wrap, [], [], [], [], false);
    delete(hl)
    wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-%s', modellabel), wrap, [], [], [], [], true);
    delete(hSEP)
    wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-%s-NOSEP', modellabel), wrap, [], [], [], [], true);
    hl = legend(gca, hMCMC, legendLabelsMCMC, 'location', 'southoutside', 'fontsize', 14, 'box', 'off');
    hl.NumColumns = maxShowHorizons;
    wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-%s-NOSEP-WITHLEGEND', modellabel), wrap, [], [], [], [], false);


    %% plot width of uncertainty bands -- only SEP
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    sepNaN = all(isnan(uncertaintySEP), 2);
    if doCharts
        sepLine = 'd';
    else
        sepLine = ':';
    end
    hSEP   = plot(sepdates(~sepNaN),uncertaintySEP(~sepNaN,:), sepLine, 'LineWidth', 4);

    for nn = 1 : maxShowHorizons
        hSEP(nn).Color  = Colors4Plots(nn);
    end

    switch datalabel
        case 'RGDP'
            ylim([0 12])
        case 'UNRATE'
            ylim([0 8])
        case 'CPI'
            ylim([0 7])
    end

    xticks(forecastOriginDates(1:4:end))
    xlim(forecastOriginDates([1 end]))
    if doDateAxis
        datetick('x', 'yyyy', 'keeplimits', 'keepticks')
    end

    grid on

    if doTitle
        title(sprintf('%s', datalabel))
    end

    legendLabels    = cat(1, legendLabelsSEP);
    hanni           = transpose(hSEP);
    hl = legend(gca, hanni(:), legendLabels(:), 'location', 'southoutside', 'fontsize', 14, 'box', 'off');
    hl.NumColumns = maxShowHorizons;

    orient landscape
    wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfanchartsSEPONLY-%s-WITHLEGEND', modellabel), wrap, [], [], [], [], false);


    %% Loop over horizons
    for hh = 1 : maxShowHorizons
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        nanNdx = isnan(uncertaintyMCMC(:,hh));
        hMCMC = plot(forecastOriginDates(~nanNdx),uncertaintyMCMC(~nanNdx,hh), '-d', 'LineWidth', 2);
        sepNaN = all(isnan(uncertaintySEP(:,hh)), 2);
        if doCharts
            sepLine = 'd';
        else
            sepLine = ':';
        end
        hSEP   = plot(sepdates(~sepNaN),uncertaintySEP(~sepNaN,hh), sepLine, 'LineWidth', 4);

        %     switch datalabel
        %         case 'RGDP'
        %             ylim([0 12])
        %         case 'UNRATE'
        %             ylim([0 8])
        %         case 'CPI'
        %             ylim([0 7])
        %     end
        ylim([0 max(ylim)])

        xticks(forecastOriginDates(1:4:end))
        xlim(forecastOriginDates([1 end]))
        if doDateAxis
            datetick('x', 'yyyy', 'keeplimits', 'keepticks')
        end

        grid on

        if doTitle
            title(sprintf('%s, h=%d', datalabel, hh-1))
        end

        legendLabels    = cat(1, legendLabelsMCMC(hh), legendLabelsSEP(hh));
        hanni           = transpose([hMCMC hSEP]);
        hl = legend(gca, hanni(:), legendLabels(:), 'location', 'southoutside', 'fontsize', 14, 'box', 'off');
        hl.NumColumns = 2;

        orient landscape
        wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-hh%d-%s-WITHLEGEND', hh, modellabel), wrap, [], [], [], [], false);
    end

end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
