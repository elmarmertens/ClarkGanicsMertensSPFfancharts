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

doStart1981 = false;


%% loop over variables
for d = 1 : length(DATALABELS)



    %% setup latexwrapper

    titlename = sprintf('SPFcdfRanges-%s', DATALABELS{d});

    wrap = [];
    initwrap;

    %% load data
    matPoint               = matfile(fullfile(datadir,sprintf('kensington%sdata',DATALABELS{d})));
    matProb                = matfile(fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS',DATALABELS{d})));

    SPFactualPointforecast = matPoint.YactualForecasts;
    SPFhistograms          = matProb.SPFhistograms;

    dates     = matProb.dates;
    if dates ~= matPoint.dates
        error('dates mismatch')
    end
    T         = size(dates,1);
    Nbar      = matPoint.Nbar;

    %% pick start of sample
    if doStart1981
        ndxStart = find(dates == datenum(1981,7,1));
    else
        ndxStart = find(dates == datenum(1992,1,1));
    end
    ndx2009Q2 = find(dates == datenum(2009,4,1));


    %% plot CDF
    hz = 2; % next year only
    for t = ndxStart : T
        if ~isempty(SPFhistograms(t).Nbins)

            %% plot
            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)

            hMid    = plot(SPFhistograms(t).binEdges, SPFhistograms(t).cdf(:,hz), 'kd','LineWidth', 2);
            hUpper  = plot(SPFhistograms(t).upperEdges, SPFhistograms(t).upperCDF(:,hz), 'b--','LineWidth', 2);
            hLower  = plot(SPFhistograms(t).lowerEdges, SPFhistograms(t).lowerCDF(:,hz), 'r--','LineWidth', 2);

            title(sprintf('y=%d', hz-1))
            ylim([0 1])
            yticks(0 : .25 : 1)
            theseTicks = round(SPFhistograms(t).binEdges, 1);
            xticks(theseTicks)
            xlim(theseTicks([1 end]))
            grid on
            title(sprintf('\\bf %s per %s (next year)',  DATALABELS{d}, datestr(dates(t), 'yyyyQQ')), 'FontSize', fontsize)
            orient landscape
            wrapthisfigure(thisfig, sprintf('%s-%s-cdf-nextyear', DATALABELS{d}, datestr(dates(t), 'yyyyQQ')), wrap);
            close(thisfig)
        end
    end

    for t = ndx2009Q2 : T
        if SPFhistograms(t).Nhz > 2
            thisfig = figure;
            for hz = 3 : SPFhistograms(t).Nhz % skip current year
                subplot(1,SPFhistograms(t).Nhz-2,hz-2)
                hold on
                set(gca, 'FontSize', fontsize)

                hMid    = plot(SPFhistograms(t).binEdges, SPFhistograms(t).cdf(:,hz), 'kd','LineWidth', 2);
                hUpper  = plot(SPFhistograms(t).upperEdges, SPFhistograms(t).upperCDF(:,hz), 'b--','LineWidth', 2);
                hLower  = plot(SPFhistograms(t).lowerEdges, SPFhistograms(t).lowerCDF(:,hz), 'r--','LineWidth', 2);

                title(sprintf('y=%d', hz-1))
                ylim([0 1])
                yticks(0 : .25 : 1)
                theseTicks = round(SPFhistograms(t).binEdges, 1);
                xticks(theseTicks)
                xlim(theseTicks([1 end]))
                grid on
            end
            sgtitle(sprintf('\\bf %s per %s',  DATALABELS{d}, datestr(dates(t), 'yyyyQQ')), 'FontSize', fontsize)
            orient landscape
            wrapthisfigure(thisfig, sprintf('%s-%s-cdf', DATALABELS{d}, datestr(dates(t), 'yyyyQQ')), wrap);
            close(thisfig)
        end
    end

    finishwrap

end % d

%% finish
dockAllFigures



