%% SPF moments, bounds on mean


%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/

%#ok<*UNRCH>

%% clear workspace
clear variables

close all
fclose all;
clc


%% parameters
datadir     = fullfile('..', 'kensingtonDataMatfiles');
DATALABELS  = {'RGDP','PGDP','UNRATE'};

yearlabels  = {'current year','next year','2 years ahead','3 years ahead'};

fontsize    = 18;

doStart1992 = true;
doTitle     = true;
%% setup latexwrapper

titlename = 'SPFmoments';

wrap = [];
initwrap;

quantilesP = [normcdf(-1) .25 .50 .75 normcdf(1)];
Nquantiles = length(quantilesP);
ndx68      = [1 5];
ndx15      = 1;
ndx85      = 5;
ndx25      = 2;
ndx50      = 3;
ndx75      = 4;

%% import mean-of-growth-rates data from Todd -- RGDP only, and only for hz=2 and hz=3
mg            = importdata(fullfile('..', 'kensingtonDataUSSPF', 'meangrowthSPF.txt'));
mg.dates      = datenum(mg.textdata(2:end,1), 'yyyy:QQ');
mg.meangrowth = mg.data(:,[1 3]);
mg.growthmean = mg.data(:,[2 4]);

doShowMG = false;
doShowOutlier = false;

%% loop over variables
for d = 1 % : 2 % length(DATALABELS)

    datalabel              = DATALABELS{d};
    matPoint               = matfile(fullfile(datadir,sprintf('kensington%sdata',datalabel)));
    matProb                = matfile(fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS',datalabel)));

    SPFactualPointforecast = matPoint.YactualForecasts;
    SPFhistograms          = matProb.SPFhistograms;

    dates     = matProb.dates;
    if dates ~= matPoint.dates
        error('dates mismatch')
    end
    T         = size(dates,1);
    Nbar      = matPoint.Nbar;

    %% data check vs mg file
    if strcmpi(datalabel, 'RGDP')
        mg.ndx   = ismember(dates, mg.dates);
        delta = SPFactualPointforecast(mg.ndx,5+(2:3)) - mg.growthmean;
        checkdiff(delta);
    end
    %% hard code LB/UB for support
    % RGDP: [-5; 8]
    % PGDP: [-1; 10]
    % UNRATE: [2; 16]

    switch datalabel
        case 'RGDP'
            binSupportLB = -5;
            binSupportUB = 8;
        case 'PGDP'
            binSupportLB = -1;
            binSupportUB = 10;
        case 'UNRATE'
            binSupportLB = 2;
            binSupportUB = 16;
    end

    %% pick start of sample
    if doStart1992
        ndxStart = find(dates == datenum(1992,1,1));
    else
        ndxStart = 1;
    end

    %% compute lower/upper mean
    [meanMid, meanLB, meanUB ] = deal(NaN(T,Nbar+1));

    for t = 1 : T

        if ~isempty(SPFhistograms(t).probs)

            binEdgesLeft  = cat(1, binSupportLB, SPFhistograms(t).binEdges(1:end-1));

            binEdgesRight = cat(1, SPFhistograms(t).binEdges(1:end-1), binSupportUB);

            theseNhz              = SPFhistograms(t).Nhz;
            meanMid(t,1:theseNhz) = sum(SPFhistograms(t).binCenter .* SPFhistograms(t).probs);
            meanLB(t,1:theseNhz)  = sum(binEdgesLeft .* SPFhistograms(t).probs);
            meanUB(t,1:theseNhz)  = sum(binEdgesRight .* SPFhistograms(t).probs);

        end
    end


    %% plot lower/upper mean
    for hz = 1 + (1 : Nbar) % skip current year
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        hMid = plot(dates, meanMid(:,hz), '--', 'Color',  Colors4Plots(7), 'LineWidth', 2);
        plot(dates, meanLB(:,hz), '-', 'Color',  Colors4Plots(7), 'LineWidth', 2);
        plot(dates, meanUB(:,hz), '-', 'Color',  Colors4Plots(7), 'LineWidth', 2);

        SPFpoint = SPFactualPointforecast(:,5+hz-1);
        hPoint   = plot(dates, SPFpoint, ':', 'Color',  Colors4Plots(6), 'LineWidth', 3);

        if doShowOutlier
            ndx = (SPFpoint < meanLB(:,hz)) | (SPFpoint > meanUB(:,hz));
            plot(dates(ndx), SPFpoint(ndx), 's', 'Color',  Colors4Plots(1), 'LineWidth', 3, 'MarkerSize', 10);
        end

        hl = legend([hMid,hPoint],{'Mid-range of SPF bins', 'SPF point forecast'}, 'location', 'best');

        if hz > 2
            xlim([datenum(2009,1,1) dates(end)]);
            xticks(datenum(1990:2:2022,1,1))
        else
            xlim([dates(ndxStart) dates(end)]);
            xticks(datenum(1995:5:2020,1,1))
        end

        if doShowMG && strcmpi(datalabel, 'RGDP') && hz > 2
            ndxUB = meanUB(mg.ndx,hz) < mg.meangrowth(:,hz-2);

            hMG = plot(mg.dates, mg.meangrowth(:,hz-2), '-.', 'Color',  Colors4Plots(8), 'LineWidth', 2);
            plot(mg.dates(ndxUB), mg.meangrowth(ndxUB,hz-2), 'd', 'Color',  Colors4Plots(8), 'LineWidth', 3, 'MarkerSize', 12);
            hl = legend([hMid,hPoint, hMG],{'Mid-range of SPF bins', 'SPF point forecast', 'SPF growth of means'}, 'location', 'best');
        end

        if doTitle
            title(yearlabels{hz})
        end
        datetick('x','YYYY','keeplimits', 'keepticks');
        wrapthisfigure(thisfig, sprintf('%s-mean-WITHPOINT-hz%d-WITHLEGEND', datalabel, hz-1), wrap);
        delete(hl)
        wrapthisfigure(thisfig, sprintf('%s-mean-WITHPOINT-hz%d', datalabel, hz-1), wrap);
        title(yearlabels{hz})
        wrapthisfigure(thisfig, sprintf('%s-mean-WITHPOINT-hz%d-WITHTITLE', datalabel, hz-1), wrap);
        hl = legend([hMid,hPoint],{'Mid-range of SPF bins', 'SPF point forecast'}, 'location', 'best');
        wrapthisfigure(thisfig, sprintf('%s-mean-WITHPOINT-hz%d-WITHTITLELEGEND', datalabel, hz-1), wrap);

        ndx = SPFpoint < meanLB(:,hz);
        ndx(dates < datenum(1992,1,1)) = false;
        if any(ndx)
            hrulefill
            fprintf('%s (h=%d): %d SPF point forecasts below bin range:\n', ...
                datalabel, hz-1, sum(ndx));
            delta = abs(SPFpoint - meanLB(:,hz));
            thesedelta = delta(ndx);
            thesedates = dates(ndx);
            fprintf('mean deviation: %f \n', mean(thesedelta));
            for nn = 1 : length(thesedates)
                fprintf('LB violation in %s (by %6.4f) \n', datestr(thesedates(nn), 'yyyyQQ'), thesedelta(nn));
            end
        end

        ndx = SPFpoint > meanUB(:,hz);
        ndx(dates < datenum(1992,1,1)) = false;
        if any(ndx)
            hrulefill
            fprintf('%s (h=%d): %d SPF point forecasts above bin range:\n', ...
                datalabel, hz-1, sum(ndx));
            delta = abs(SPFpoint - meanUB(:,hz));
            thesedelta = delta(ndx);
            thesedates = dates(ndx);
            fprintf('mean deviation: %f \n', mean(thesedelta));
            for nn = 1 : length(thesedates)
                fprintf('UB violation in %s (by %6.4f) \n', datestr(thesedates(nn), 'yyyyQQ'), thesedelta(nn));
            end
        end
    end

end

%% finish
dockAllFigures
finishwrap



