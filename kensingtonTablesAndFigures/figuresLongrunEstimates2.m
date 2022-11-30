%% shifting endpoints of term structures of expectations, SV and CONST models

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

fontsize = 18;
doAVG5 = false; % CPI only

%% some parameters

modeltype   = 'STATEtrendgapSV';
modelpretty = 'SV';
Ndraws      = 3e3;

modeltype2   = 'STATEconst';
modelpretty2 = 'CONST';
Ndraws2      = 3e3;

resultdir = localresultsMCMC;

fontsize  = 24;
doFinal   = false;

DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'TBILL', 'CPI'};

%% set up wrapper
titlename = sprintf('longrunestimates2-%s-vs-%s', modeltype, modeltype2);
wrap = [];
initwrap

%% loop over datalabels

for d = 1 : length(DATALABELS)

    %     close all
    datalabel   = DATALABELS{d};
    modellabel  = strcat(datalabel, modeltype);
    modellabel2 = strcat(datalabel, modeltype2);


    fprintf('Processing %s ... \n', datalabel)

    %% load data
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel2, Ndraws2);
    mat2 = matfile(fullfile(resultdir, matfilename));



    Ydraws        = mat.Ydraws;
    YFINALdraws   = mat.YFINALdraws;
    yrealized     = mat.Yfuture(:,1);

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

    % model2
    Ydraws2        = mat2.Ydraws;
    YFINALdraws2   = mat2.YFINALdraws;
    if ~isequaln(yrealized, mat2.Yfuture(:,1))
        error('yrealized mismatch')
    end
    if ~isequaln(Zdata, mat2.Zdata)
        warning('Zdata mismatch')
    end

    %% import RGDP10 / CPI10 etc. (if appropriate)
    switch upper(datalabel)
        case 'RGDP'
            spfdata10 = importdata('../kensingtonDataUSSPF/Mean_RGDP10_Level.xlsx');
            spfdates  = datenum(spfdata10.data(:,1), (spfdata10.data(:,2)-1)*3+1, 1);
            AVG10    = log(1 + spfdata10.data(:,3)/100) * 100;
            clear spfdata10
            AVG5  = [];
        case 'TBILL'
            spfdata10 = importdata('../kensingtonDataUSSPF/Mean_BILL10_Level.xlsx');
            spfdates  = datenum(spfdata10.data(:,1), (spfdata10.data(:,2)-1)*3+1, 1);
            AVG10     = spfdata10.data(:,3);
            clear spfdata10
            AVG5  = [];
        case 'CPI'
            spfdata10 = importdata('../kensingtonDataUSSPF/Mean_CPI10_Level.xlsx');
            spfdates  = datenum(spfdata10.data(:,1), (spfdata10.data(:,2)-1)*3+1, 1);
            AVG10     = spfdata10.data(:,3);
            clear spfdata10
            if doAVG5
                spfdata5  = importdata('../kensingtonDataUSSPF/Mean_CPI5YR_Level.xlsx');
                spfdates5 = datenum(spfdata5.data(:,1), (spfdata5.data(:,2)-1)*3+1, 1);
                AVG5      = spfdata5.data(:,3);
                clear spfdata5
                checkdiff(spfdates5, spfdates);
                clear spfdates5
            else
                AVG5  = [];
            end
        otherwise
            AVG10 = [];
            AVG5  = [];
    end


    %% compare qrt estimates
    yrealizedMA4 = sumK(yrealized,4,1) / 4;

    ystardraws = squeeze(Ydraws(end,:,:));
    mid        = median(ystardraws,2);
    tails      = prctile(ystardraws, normcdf([-1 1]) * 100, 2);

    ystardraws = squeeze(Ydraws2(end,:,:));
    mid2       = median(ystardraws,2);
    tails2     = prctile(ystardraws, normcdf([-1 1]) * 100, 2);

    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    h    = plot(dates, mid, 'k-', 'linewidth', 2);
    %     plot(dates, tails, 'k', 'linewidth', 1);

    h2    = plot(dates, mid2, 'r-.', 'linewidth', 2);
    %     plot(dates, tails2, 'r--', 'linewidth', 1);

    hrealized = plot(dates, yrealizedMA4, 'b-', 'linewidth', 1);
    if ~isempty(AVG10)
        h10 = plot(spfdates(~isnan(AVG10)), AVG10(~isnan(AVG10)), '-', 'color', Colors4Plots(5), 'linewidth', 2);
        %         if iscompact(AVG10)
        %             h10 = plot(spfdates, AVG10, '--', 'color', Colors4Plots(5), 'linewidth', 2);
        %         else
        %             h10 = plot(spfdates(~isnan(AVG10)), AVG10(~isnan(AVG10)), 'd', 'color', Colors4Plots(5), 'linewidth', 2);
        %         end
    end
    if ~isempty(AVG5)
        if iscompact(AVG5)
            h5 = plot(spfdates, AVG5, ':', 'color', Colors4Plots(3), 'linewidth', 2);
        else
            h5 = plot(spfdates(~isnan(AVG5)), AVG5(~isnan(AVG5)), 'd', 'color', Colors4Plots(3), 'linewidth', 2);
        end
    end
    xticks(datenum(1970:10:2020,1,1))
    xtickdates(dates, 'keepticks')
    wrapthisfigure(thisfig, sprintf('ystarQRT-%s-%s-vs-%s', ...
        datalabel, modeltype, modeltype2), wrap, [], [], [], [], true);
    if ~isempty(AVG5)
        legend([h h2 hrealized h10 h5], sprintf('y* %s', modelpretty), sprintf('y* %s', modelpretty2), ...
            'y (MA4)', 'AVG10', 'AVG05', 'Location', 'best')
    elseif ~isempty(AVG10)
        legend([h h2 hrealized h10], sprintf('y* %s', modelpretty), sprintf('y* %s', modelpretty2), ...
            'y (MA4)', 'AVG10', 'Location', 'best')
    else
        legend([h h2 hrealized], sprintf('y* %s', modelpretty), sprintf('y* %s', modelpretty2), ...
            'y (MA4)', 'Location', 'best')
    end
    wrapthisfigure(thisfig, sprintf('ystarQRT-%s-%s-vs-%s-WITHLEGEND', ...
        datalabel, modeltype, modeltype2), wrap, [], [], [], [], false);
    %     delete(hrealized)
    %     wrapthisfigure(thisfig, sprintf('ystarQRT-%s-%s-vs-%s-NODATA-WITHLEGEND', ...
    %         datalabel, modeltype, modeltype2), wrap, [], [], [], [], true);

    %% compare final estimates
    if doFinal
        yrealizedMA4 = sumK(yrealized,4,1) / 4;

        ystardraws = squeeze(YFINALdraws(end,:,:));
        mid        = median(ystardraws,2);
        tails      = prctile(ystardraws, normcdf([-1 1]) * 100, 2);

        ystardraws = squeeze(YFINALdraws2(end,:,:));
        mid2       = median(ystardraws,2);
        tails2     = prctile(ystardraws, normcdf([-1 1]) * 100, 2);

        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        h    = plot(dates, mid, 'k-', 'linewidth', 2);
        %     plot(dates, tails, 'k', 'linewidth', 1);

        h2    = plot(dates, mid2, 'r-.', 'linewidth', 2);
        %     plot(dates, tails2, 'r--', 'linewidth', 1);

        hrealized = plot(dates, yrealizedMA4, 'b-', 'linewidth', 1);

        if ~isempty(AVG10)

            h10 = plot(spfdates, AVG10, 'd', 'color', Colors4Plots(5), 'linewidth', 1);
            legend([h h2 hrealized h10], sprintf('y* %s', modelpretty), sprintf('y* %s', modelpretty2), ...
                'y (MA4)', 'RGDP10', 'Location', 'best')
        else
            legend([h h2 hrealized], sprintf('y* %s', modelpretty), sprintf('y* %s', modelpretty2), ...
                'y (MA4)', 'Location', 'best')
        end
        xticks(datenum(1970:10:2020,1,1))
        xtickdates(dates, 'keepticks')
        wrapthisfigure(thisfig, sprintf('ystarfinal-%s-%s-vs-%s-WITHLEGEND', ...
            datalabel, modeltype, modeltype2), wrap, [], [], [], [], false);
        delete(hrealized)
        wrapthisfigure(thisfig, sprintf('ystarfinal-%s-%s-vs-%s-NODATA-WITHLEGEND', ...
            datalabel, modeltype, modeltype2), wrap, [], [], [], [], false);
    end

end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
