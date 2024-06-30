%% compares terms structure of expectations from two models

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
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc

%% some parameters



modeltype   = 'trendHScycleSVt2blockNoiseHS-y1q4-NgapBOP-samStart1968Q4';
Ndraws      = 3e3;

resultdir = '~/jam/lager/kensingtonStateSpaceDraws/';


DATALABELS = {'RGDP', 'UNRATE', 'CPI', 'PGDP'};

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx90     = [2 8];
fontsize  = 18;

samStart  = datenum(1990,1,1);

forecastOrigins = [datenum(2019,10,1); datenum(2024,1,1)];

%% loop over MDStype

for MDSTYPES = {'MDS'} % {'MDS', 'VAR0'}

    MDStype = MDSTYPES{1};

    %% prepare latexwrapper
    wrap = [];
    titlename = sprintf('termstructures-%s-%s', MDStype, modeltype);
    initwrap
    if isempty(wrap) && ~isdesktop
        initwrap
    end

    %% loop over datalabels
    for d = 1 : length(DATALABELS)

        datalabel   = DATALABELS{d};


        modellabel  = strcat(datalabel, '-', MDStype, modeltype);



        fprintf('Processing %s ... \n', datalabel)

        %% load data

        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
        mat = matfile(fullfile(resultdir, matfilename));


        Yhat           = mat.fcstYhat;
        Yquantiles     = mat.fcstYquantiles;

        Yfuture       = mat.Yfuture;
        Nz            = mat.Nz;
        Zdata         = mat.Zdata;
        Znanny        = mat.Znanny;

        Nhorizons     = mat.Nhorizons;
        dates         = mat.dates;
        datesQ        = mat.datesQ;
        T             = mat.T;
        doNIPA        = mat.doNIPA;
        Tstart        = mat.Tstart;

        Tstart = max(Tstart, find(dates == samStart));

        % extend dates by Nhorizons
        [lastY, lastM, lastD] = datevec(dates(end));
        datesX = [dates; datenum(lastY, lastM + 3 + 3 *(0:Nhorizons-1)', lastD)];

        ndxFixedHorizon     = 1 + 1 + (0:4);


        %% collect statistics
        Ytail68     = Yquantiles(:,:,ndx68);
        Ytail90     = Yquantiles(:,:,ndx90);


        %% loop over forecastOrigins
        Tset = transpose(find(ismember(dates, forecastOrigins)));
        colorSPF1 = colors4plots('black');
        colorSPF2 = colors4plots('blue');
        
        colorCGM = colors4plots('orange');

        for thisT = Tset

            thisDate = datestr(dates(thisT), 'yyyyqq');

            mid     = transpose(Yhat(thisT,:));
            tail68  = squeeze(Ytail68(thisT,:,:));
            tail90  = squeeze(Ytail90(thisT,:,:));

            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)
            % plot density
            hDense = plot(0:Nhorizons-1, mid, '-', 'color', colorCGM, 'LineWidth', 3);
            plot(0:Nhorizons-1,tail68, '-.', 'color', colorCGM, 'LineWidth', 2)
            plot(0:Nhorizons-1,tail90, '--', 'color', colorCGM, 'LineWidth', 1)
            % add SPF
            if doNIPA
                MAweights       = -2 : 4;
            else
                MAweights       = 1 : 4;
            end
            ndxB            = 2 + 4 + 1;
            horizonsB       = (4 - datesQ(thisT)) + MAweights;
            ndxC            = 2 + 4 + 2;
            horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
            ndxD            = 2 + 4 + 3;
            horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
            hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'd', 'color', colorSPF1, 'linewidth', 3, 'MarkerSize', 15);
            if ~Znanny(thisT,ndxB)
                hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), 's', 'color', colorSPF2, 'linewidth', 3, 'MarkerSize', 15);
                plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':', 'color', colorSPF2, 'linewidth', 4);
            end
            if ndxC <= Nz && ~Znanny(thisT,ndxC)
                hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), 's', 'color', colorSPF2, 'linewidth', 3, 'MarkerSize', 15);
                plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':', 'color', colorSPF2, 'linewidth', 4);
            end
            if ndxD <= Nz && ~Znanny(thisT,ndxD)
                hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), 's', 'color', colorSPF2, 'linewidth', 3, 'MarkerSize', 15);
                plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':', 'color', colorSPF2, 'linewidth', 4);
            end

            ndxCalendarYears = (4 - datesQ(thisT)) + 1:4:12;
            xline(ndxCalendarYears, ':')
            % axis settings
            xmax = Nhorizons - 1;
            xticks(0 : 2 : xmax)
            xlim([0 xmax])
            wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s', modellabel, thisDate), wrap, [], [], [], [], true)
            legend([hDense hSPF hSPFb(1)], 'Predictive Density', 'SPF fixed-horizon data', 'SPF fixed-event data', 'location', 'best');
            wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGEND', modellabel, thisDate), wrap, [], [], [], [], true)
            title(sprintf('%s: %s per %s', MDStype, datalabel, thisDate))
            wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGENDTITLE', modellabel, thisDate), wrap)


        end % thisT

    end % datalabels


    dockAllFigures
    finishwrap

end % MDStype

%% finish / clean up
finishscript
dockAllFigures
