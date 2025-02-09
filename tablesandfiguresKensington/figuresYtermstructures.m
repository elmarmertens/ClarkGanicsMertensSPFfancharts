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
fontsize  = 18;

samStart  = datenum(1990,1,1);

%% loop over MDStype

for MDSTYPES = {'MDS', 'VAR0'}

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


        Ydraws         = mat.Ydraws;
        Yquantiles     = mat.fcstYquantiles;

        Yfuture       = mat.Yfuture;
        Nz            = mat.Nz;
        Zdata         = mat.Zdata;
        Znanny        = mat.Znanny;

        Nhorizons     = mat.Nhorizons;
        dates         = mat.dates;
        T             = mat.T;
        doNIPA        = mat.doNIPA;
        Tstart        = mat.Tstart;

        Tstart = max(Tstart, find(dates == samStart));

        % extend dates by Nhorizons
        [lastY, lastM, lastD] = datevec(dates(end));
        datesX = [dates; datenum(lastY, lastM + 3 + 3 *(0:Nhorizons-1)', lastD)];

        %% patch in LR SPF
        switch upper(datalabel)
            case 'RGDP'
                spfdata10  = importdata('../rawdataKensington/mean_RGDP10_level.xlsx');
                SPFdates   = datenum(spfdata10.data(:,1), (spfdata10.data(:,2)-1)*3+1, 1);
                SPFlongrun = log(1 + spfdata10.data(:,3)/100) * 100;
                clear spfdata10
            case 'CPI'
                spfdata10  = importdata('../rawdataKensington/mean_CPI10_level.xlsx');
                SPFdates   = datenum(spfdata10.data(:,1), (spfdata10.data(:,2)-1)*3+1, 1);
                SPFlongrun = spfdata10.data(:,3);
                clear spfdata10
            otherwise
                SPFlongrun = [];
                SPFdates   = [];
        end

        %% collect statistics
        SPFhat      = transpose(mean(Ydraws, 3)); % T x H
        SPFhatBand  = permute(prctile(Ydraws, normcdf11, 3), [2 3 1]); % T x Bands x H
        band68width = Yquantiles(:,:,ndx68(2)) - Yquantiles(:,:,ndx68(1));

        % pad SPFhat with trend (if needed)
        if size(SPFhat, 2) < Nhorizons
            SPFhat = [SPFhat, repmat(SPFhat(:,end), 1, Nhorizons - size(SPFhat, 2))]; %#ok<AGROW>
        end

        %% plot term structure of expectations separately by quarter INCL LR SPF
        if isempty(SPFlongrun)
            % T1985 = find(dates == datenum(1985,1,1));
            % T2020 = find(dates == datenum(2020,1,1));



            % plot endpoint of TS vs LR SPF
            tickDates = dates(Tstart : 20 : T);
            thisfig = figure;
            hold on
            hCGM   = plot(dates(Tstart:T), SPFhat(Tstart:T,end), 'linewidth', 2, 'color', colors4plots('darkblue'));
            hCGM68 = plot(dates(Tstart:T), SPFhatBand(Tstart:T,:,end), 'LineWidth', 1, 'color', colors4plots('darkblue'));
            % if iscompact(SPFlongrun(Tstart:T))
            %     linestyle = '-.';
            % else
            %     linestyle = 'd';
            % end
            % hSPF = plot(dates(Tstart:T), SPFlongrun(Tstart:T), linestyle, 'linewidth', 2, 'color', colors4plots('darkred'));
            plotOrigin
            % xticks(tickDates)
            % xtickdates(dates([Tstart T]), 'keepticks')
            nbershades(dates([Tstart T]), [], [], tickDates)

            set(gca, 'fontsize', fontsize)
            wrapthisfigure(thisfig, sprintf('termstructureLR-%s-%s-EndPoints', MDStype, datalabel), wrap);
            legend([hCGM, hCGM68(1)], 'Endpoints of fitted term structures', '68% band for estimated endpoints', 'location', 'southwest')
            wrapthisfigure(thisfig, sprintf('termstructureLR-%s-%s-EndPoints-WITHLEGEND', MDStype, datalabel), wrap);

        else % add SPFLR
            % T1985 = find(dates == datenum(1985,1,1));
            % T2020 = find(dates == datenum(2020,1,1));



            % plot endpoint of TS vs LR SPF
            tickDates = dates(Tstart : 20 : T);
            thisfig = figure;
            hold on
            hCGM   = plot(dates(Tstart:T), SPFhat(Tstart:T,end), 'linewidth', 2, 'color', colors4plots('darkblue'));
            hCGM68 = plot(dates(Tstart:T), SPFhatBand(Tstart:T,:,end), 'LineWidth', 1, 'color', colors4plots('darkblue'));
            if iscompact(SPFlongrun(Tstart:T))
                linestyle = '-.';
            else
                linestyle = 'd';
            end
            hSPF = plot(dates(Tstart:T), SPFlongrun(Tstart:T), linestyle, 'linewidth', 2, 'color', colors4plots('darkred'));
            plotOrigin
            % xticks(tickDates)
            % xtickdates(dates([Tstart T]), 'keepticks')
            nbershades(dates([Tstart T]), [], [], tickDates)

            set(gca, 'fontsize', fontsize)
            wrapthisfigure(thisfig, sprintf('termstructureLR-%s-%s-EndPoints', MDStype, datalabel), wrap);
            legend([hCGM, hCGM68(1), hSPF], 'Endpoints of fitted term structures', '68% band for estimated endpoints', 'SPF long-run forecasts', 'location', 'southwest')
            wrapthisfigure(thisfig, sprintf('termstructureLR-%s-%s-EndPoints-WITHLEGEND', MDStype, datalabel), wrap);


        end % SPF

        %% plot term structure of expectations separately by quarter (without LR SPF)
        % T1985 = find(dates == datenum(1985,1,1));
        T2020 = find(dates == datenum(2020,1,1));
        for Tstop = [find(dates == datenum(2019,10,1)) T]
            tickDates = dates(Tstart : 20 : max(T2020, Tstop));
            for offset = 0 : 3
                thisfig = figure;
                hold on
                ii = 0;
                for tt = Tstart + offset : 4 :  Tstop
                    ii = ii + 1;
                    thesedates = datesX(tt + (0 : Nhorizons - 1));
                    hanni(ii) = plot(thesedates, SPFhat(tt,:), 'linewidth', 2);
                end
                hRealized = plot(dates(Tstart+offset:Tstop), Yfuture(Tstart+offset:Tstop,1), 'k-.', 'LineWidth', 1);
                plotOrigin
                % xticks(tickDates)
                % xtickdates(dates([Tstart+offset max(T2020, Tstop)]), 'keepticks')
                nbershades(dates([Tstart max(T2020, Tstop)]), [], [], tickDates)

                set(gca, 'fontsize', fontsize)
                wrapthisfigure(thisfig, sprintf('termstructure-%s-%s-Q%d-Until%s', MDStype, datalabel, quarter(dates(Tstart+offset)), datestr(dates(Tstop), 'yyyyQQ')), wrap);
                legend([hRealized, hanni(1)], 'Realized', 'Fitted term structures', 'location', 'best')
                wrapthisfigure(thisfig, sprintf('termstructure-%s-%s-Q%d-Until%s-WITHLEGEND', MDStype, datalabel, quarter(dates(Tstart+offset)), datestr(dates(Tstop), 'yyyyQQ')), wrap);
            end
        end


        %% plot selected lines for uncertainty bands


        Tstop = T;
        Horizons = [0, 1, 4 : 4 : Nhorizons-1];
        hLabel = arrayfun(@(x) sprintf('h=%d', x), Horizons, 'uniformoutput', false);
        thisfig = figure;
        hanni = plot(dates(Tstart:Tstop), band68width(Tstart:Tstop,Horizons+1), 'linewidth', 2);
        LineStyles = {'-', '-', '--', ':', '-.','-'};
        Colors     = {'blue', 'orange', 'yellow', 'purple', 'green', 'lightblue'};
        MarkerStyles = {'d', 'none', 'none', 'none', 'none', 's'};
        for i = 1:length(hanni)
            hanni(i).LineStyle = LineStyles{i};
            hanni(i).Color = colors4plots(Colors{i});
            if ~strcmpi(MarkerStyles{i}, 'none')
            hanni(i).Marker = MarkerStyles{i};
            hanni(i).MarkerSize = 8;
            hanni(i).LineWidth  = 1;
            end
        end
        % xticks(tickDates)
        % xtickdates(dates([Tstart Tstop]), 'keepticks')
        nbershades(dates([Tstart Tstop]), [], [], tickDates)

        legend(hanni, hLabel, 'Location','best')
        set(gca, 'fontsize', fontsize)
        set(gca, 'box', 'off')
        wrapthisfigure(thisfig, sprintf('uncertainty-%s-%s', MDStype, datalabel), wrap);

    end % datalabels


    dockAllFigures
    finishwrap

end % MDStype

%% finish / clean up
finishscript
dockAllFigures
