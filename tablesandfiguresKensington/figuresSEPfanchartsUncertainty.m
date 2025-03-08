%% Plot evoluation of Uncertainty from CGM model and compare against SEP
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
%#ok<*DATIC>
%#ok<*DATNM>

%% clear workspace
clear variables

close all
fclose all;
clc

%% some parameters

DATALABELS = {'RGDP', 'UNRATE', 'CPI'};


doDateAxis = true;
doTitle    = false;

modeltype   = 'trendHScycleSVt2blockNoiseHS-y1q4-NgapBOP-samStart1968Q4';
Ndraws      = 3e3;

resultdir = '../mcmcKensington/foo/';

%% process some parameters

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx50     = 5;
fontsize  = 22;

%% loop over MDStypes

for MDStype = {'MDS'} % {'MDS', 'VAR0'}

    switch MDStype{:}
        case 'MDS'
            MDSpretty = 'MDS';
        case 'VAR0'
            MDSpretty = 'VAR';
        otherwise
            error('MDStype <<%s>> not yet implemented', MDStype{:})
    end

    modelpretty = strcat('CGM-', MDSpretty);

    colorMCMC = colors4plots(1);
    colorData = colors4plots(7);

    %% import SEP tables
    SEPrmseData          = importdata('../rawdataKensington/SEPTable2Data.xlsx');

    SEPabserrorsTable    = readtable('../rawdataKensington/SEPforecasterrors_absvalues.xlsx','Sheet','abs errors');
    SEPabserrorsDates    = datenum(SEPabserrorsTable.Var1, 'yyyy:QQ');
    SEPabserrorsVarnames = SEPabserrorsTable.Properties.VariableNames;

    SEPerrorsTable    = readtable('../rawdataKensington/SEPforecasterrors_absvalues.xlsx','Sheet','raw errors');
    SEPerrorsDates    = datenum(SEPabserrorsTable.Var1, 'yyyy:QQ');
    SEPerrorsVarnames = SEPabserrorsTable.Properties.VariableNames;

    if ~isequal(SEPabserrorsDates, SEPerrorsDates)
        error('dates in abs errors and raw errors do not match')
    end
    if ~isequal(SEPabserrorsVarnames, SEPerrorsVarnames)
        error('varnames in abs errors and raw errors do not match')
    end

    %% prepare latexwrapper
    wrap = [];
    titlename = strcat('SEPfanchartsUncertainty-', modelpretty, modeltype);
    initwrap
    if isempty(wrap) && ~isdesktop
        initwrap
    end

    %% loop over datalabels
    for d = 1 : 3

        tndx = 158 : 221;        % 2011 - end of sample

        datalabel = DATALABELS{d};

        modellabel = strcat(datalabel, '-', MDStype{:}, modeltype);

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
            SEPbandsDates = datenum(SEPrmseData.textdata.(sheetname)(2:end,1));
            SEPrmse       = SEPrmseData.data.(sheetname);
        else
            SEPbandsDates = x2mdate(SEPrmseData.data.(sheetname)(:,1)); %#ok<XMDATE>
            SEPrmse       = SEPrmseData.data.(sheetname)(:,2:end);
        end

        %% read SEP absolute errors data
        switch datalabel
            case 'RGDP'
                varNdx = contains(SEPabserrorsVarnames, 'GDP');
            case 'UNRATE'
                varNdx = contains(SEPabserrorsVarnames, 'UNRATE');
            case 'CPI'
                varNdx = contains(SEPabserrorsVarnames, 'PCEINFL');
            otherwise
                error('datalabel <<%s>> not recognized', datalabel)
        end
        SEPabserrors = SEPabserrorsTable{:,varNdx};
        SEPerrors    = SEPerrorsTable{:,varNdx};

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
                Yerrors       = mat.fcstYAVGhaterror;
            case {'UNRATE'}
                Yhat          = mat.fcstYhat;
                Yquantiles    = mat.fcstYquantiles;
                Yfuture       = mat.Yfuture;
                Yerrors       = mat.fcstYhaterror;
            otherwise
                error('datalabel <<%s>> not yet implemented')
        end

        % set maxNhorizons
        maxAnnualHorizons = 4;

        %% pad extra dates at end
        xdates = genrQdates(year(dates(end)), year(dates(end))+4, 1);
        dates  = union(dates, xdates);

        fprintf('Processing %s ... \n', modellabel)

        %% collect uncertainty data
        uncertaintyMCMC = NaN(length(tndx), maxAnnualHorizons);
        errorbandMCMC   = NaN(length(tndx), maxAnnualHorizons, 2);
        errorsMCMC      = NaN(length(tndx), maxAnnualHorizons);
        for tt = 1 : length(tndx)
            thisT           = tndx(tt);
            showQuarters    = (3 : 4 : Nhorizons) - (datesQ(thisT) - 1);
            ndxShowQuarters = showQuarters + 1;
            thisErrorBand         = Yquantiles(thisT,ndxShowQuarters,ndx68);
            uncertaintyMCMC(tt,:) = range(thisErrorBand, 3);
            thisErrorBand         = thisErrorBand - Yquantiles(thisT,ndxShowQuarters,ndx50);
            errorbandMCMC(tt,:,:) =  thisErrorBand;
            % checkdiff(uncertaintyMCMC(tt,:), range(thisErrorBand, 3));

            errorsMCMC(tt,:) = - Yerrors(thisT,ndxShowQuarters); % negative, since MCMC stored Yhat minus Yfuture
                        
        end

        %% create labels
        legendLabelsXXX    = {'this year', 'next year', 'two years', 'three years'};
        legendLabelsMCMC   = cellfun(@(x) sprintf('%s (%s)', x, modelpretty), legendLabelsXXX, 'UniformOutput', false);


        % match SEP dates
        sepndx          = ismember(SEPbandsDates, dates(tndx));
        uncertaintySEP  = 2.* SEPrmse(sepndx, :);
        legendLabelsSEP = cellfun(@(x) sprintf('%s (SEP)', x), legendLabelsXXX, 'UniformOutput', false);
        SEPbandsDates   = SEPbandsDates(sepndx);

        forecastOriginDates = dates(tndx);



        %% SEP vs CGM -- separate plots for each horizon
        for yy = 1 : maxAnnualHorizons
            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)
            nanNdx  = isnan(uncertaintyMCMC(:,yy));
            hMCMC   = plot(forecastOriginDates(~nanNdx),uncertaintyMCMC(~nanNdx,yy), '-', 'LineWidth', 4);
            % sepNaN  = all(isnan(uncertaintySEP(:,yy)), 2);
            sepLine = ':';
            hSEP    = plot(SEPbandsDates,uncertaintySEP(:,yy), sepLine, 'LineWidth', 4);

            ylim([0 max(ylim)])

            nbershades(forecastOriginDates)
            xticks(forecastOriginDates(1:12:end))
            h = gca;
            h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
            h.XAxis.MinorTickValues = forecastOriginDates(1:4:end);
            xlim(forecastOriginDates([1 end]))
            if doDateAxis
                datetick('x', 'yyyy', 'keeplimits', 'keepticks')
            end

            if doTitle
                title(sprintf('%s, y=%d', datalabel, yy-1))
            end

            % store without legend
            wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-hh%d-%s', yy, modellabel), wrap, [], [], [], [], true);

            % add legend
            legendLabels    = cat(1, legendLabelsMCMC(yy), legendLabelsSEP(yy));
            hanni           = transpose([hMCMC hSEP]);
            hl = legend(gca, hanni(:), 'CGM', 'SEP', 'location', 'best', 'fontsize', fontsize, 'box', 'on');
            wrapthisfigure(thisfig, sprintf('UNCERTAINTYsepfancharts-hh%d-%s-WITHLEGEND', yy, modellabel), wrap, [], [], [], [], false);
        end

        %% SEP vs CGM -- rmse-bands with errors
        rmsebandSEP  = .5 .* uncertaintySEP;
        rmsebandMCMC = .5 .* uncertaintyMCMC;

        for yy = 1 : maxAnnualHorizons
            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)
            nanNdx  = isnan(rmsebandMCMC(:,yy));
            hMCMC   = plot(forecastOriginDates(~nanNdx),rmsebandMCMC(~nanNdx,yy), '-', 'LineWidth', 4);
            % sepNaN  = all(isnan(rmsebandSEP(:,yy)), 2);
            sepLine = ':';
            hSEP    = plot(SEPbandsDates,rmsebandSEP(:,yy), sepLine, 'LineWidth', 4);

            hSEPerror = plot(SEPabserrorsDates, SEPabserrors(:,yy), 'd', 'color', 'black', 'linewidth', 4);

            ylim([0 max(ylim)])

            nbershades(forecastOriginDates)
            xticks(forecastOriginDates(1:12:end))
            h = gca;
            h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
            h.XAxis.MinorTickValues = forecastOriginDates(1:4:end);
            xlim(forecastOriginDates([1 end]))
            if doDateAxis
                datetick('x', 'yyyy', 'keeplimits', 'keepticks')
            end

            if doTitle
                title(sprintf('%s, y=%d', datalabel, yy-1))
            end

            % store without legend
            wrapthisfigure(thisfig, sprintf('RMSEABSERRsepfancharts-hh%d-%s', yy, modellabel), wrap, [], [], [], [], true);

            % add legend
            legendLabels    = cat(1, legendLabelsMCMC(yy), legendLabelsSEP(yy));
            hanni           = transpose([hMCMC hSEP hSEPerror]);
            hl = legend(gca, hanni(:), 'CGM', 'SEP', 'abs. SEP error', 'location', 'best', 'fontsize', fontsize, 'box', 'on');
            wrapthisfigure(thisfig, sprintf('RMSEABSERRsepfancharts-hh%d-%s-WITHLEGEND', yy, modellabel), wrap, [], [], [], [], false);
        end

        %% SEP vs CGM -- erorrs with bands
        errorbandSEP  = cat(3, -.5 .* uncertaintySEP, .5 .* uncertaintySEP);

        for yy = 1 : maxAnnualHorizons
            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)
            nanNdx  = any(isnan(errorbandMCMC(:,yy,:)),3); 
            hMCMC   = plot(forecastOriginDates(~nanNdx),squeeze(errorbandMCMC(~nanNdx,yy,:)), '-', 'color', colors4plots('blue'), 'LineWidth', 4);
            sepLine = ':';
            hSEP    = plot(SEPbandsDates,squeeze(errorbandSEP(:,yy,:)), sepLine, 'color', colors4plots('orange'), 'LineWidth', 4);

            hSEPerror = plot(SEPabserrorsDates, SEPerrors(:,yy), 'd', 'color', 'black', 'linewidth', 4);

            hSPFerror = plot(forecastOriginDates, errorsMCMC(:,yy), 'o', 'color', colors4plots('green'), 'linewidth', 2);

            nbershades(forecastOriginDates)
            xticks(forecastOriginDates(1:12:end))
            h = gca;
            h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
            h.XAxis.MinorTickValues = forecastOriginDates(1:4:end);
            xlim(forecastOriginDates([1 end]))
            if doDateAxis
                datetick('x', 'yyyy', 'keeplimits', 'keepticks')
            end

            if doTitle
                title(sprintf('%s, y=%d', datalabel, yy-1))
            end

            % store without legend
            wrapthisfigure(thisfig, sprintf('ERRORBANDsepfancharts-hh%d-%s', yy, modellabel), wrap, [], [], [], [], true);

            % add legend
            legendLabels    = cat(1, legendLabelsMCMC(yy), legendLabelsSEP(yy));
            hanni           = transpose([hMCMC(1) hSEP(1) hSEPerror hSPFerror]);
            hl = legend(gca, hanni(:), 'CGM', 'SEP', 'SEP error', 'SPF error', 'location', 'best', 'NumColumns', 2, 'fontsize', fontsize, 'box', 'on');
            wrapthisfigure(thisfig, sprintf('ERRORBANDsepfancharts-hh%d-%s-WITHLEGEND', yy, modellabel), wrap, [], [], [], [], false);
        end

    end

    finishwrap

end % MDS type
%% finish / clean up
finishscript
dockAllFigures
