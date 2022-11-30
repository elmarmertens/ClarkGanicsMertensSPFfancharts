%% compute and plot bias in SPF-consistent expectations of VAReta

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

Ndraws      = 3e3;

resultdir = localresultsMCMC;

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
