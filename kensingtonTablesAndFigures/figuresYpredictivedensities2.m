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

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};

modeltype   = 'STATEtrendgapSV';
modelpretty = 'SV';
Ndraws      = 3e3;

modeltype2   = 'STATEtrendgaplongSV';
modelpretty2 = 'SV-AVG10';
Ndraws2      = 3e3;

resultdir = localresultsMCMC;

doTight = true; % tight axis
doTitle = false;
doPanel2 = false;


quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx90     = [2 8];
fontsize  = 18;




for d = 2 % 2 4] % 1 : length(DATALABELS)

    %% prepare things

    tndx = [];

    tndx = 156 : 215;
    %          tndx = 208 % 206 : 214; % since 2020Q1

    close all
    datalabel = DATALABELS{d};

    modellabel  = strcat(datalabel, modeltype);
    modellabel2 = strcat(datalabel, modeltype2);

    %% load data

    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel2, Ndraws2);
    mat2 = matfile(fullfile(resultdir, matfilename));


    %     YdensityDraws = mat.YdensityDraws;
    Yhat          = mat.fcstYhatRB;
    Yquantiles    = mat.fcstYquantiles;
    Ydraws        = mat.Ydraws;

    % YdensityDraws2 = mat2.YdensityDraws;
    Yhat2         = mat2.fcstYhatRB;
    Yquantiles2   = mat2.fcstYquantiles;
    Ydraws2       = mat2.Ydraws;
    Nhorizons2    = mat2.Nhorizons;

    Yfuture       = mat.Yfuture;
    Nz            = mat.Nz;
    Zdata         = mat.Zdata;
    Znanny        = mat.Znanny;

    avg10         = mat2.Zdata(:,end);
    
    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    Tstart        = max(mat.Tstart, mat2.Tstart);

    % checks
    if ~isequal(dates, mat2.dates)
        error('mismatch between input files for model1 and model2')
    end
    %     if ~isequaln(Zdata, mat2.Zdata)
    %         error('mismatch between input files for model1 and model2')
    %     end
    %     if ~isequaln(Yfuture, mat2.Yfuture)
    %         error('mismatch between input files for model1 and model2')
    %     end

    %% prepare latexwrapper
    wrap = [];
    titlename = sprintf('fancharts-%s-vs-%s', modellabel, modeltype2);
    initwrap

    fprintf('Processing %s ... \n', modellabel)

    %% plot predictive density of Y


    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end

    ndxFixedHorizon     = 1 + 1 + (0:4);


    if isempty(tndx)
        if ismac && isdesktop
            tndx = 166; % T - 12;
        else
            tndx = Tstart : T;
        end
    end

    for thisT = tndx


        ndxB            = 2 + 4 + 1;
        horizonsB       = (4 - datesQ(thisT)) + MAweights;
        ndxC            = 2 + 4 + 2;
        horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
        ndxD            = 2 + 4 + 3;
        horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;


        mid     = Yhat(thisT,:);
        tails   = squeeze(Yquantiles(thisT,:,union(ndx68,ndx90)));
        tails68 = squeeze(Yquantiles(thisT,:,ndx68));
        tails90 = squeeze(Yquantiles(thisT,:,ndx90));

        mid2     = Yhat2(thisT,:);
        tails268 = squeeze(Yquantiles2(thisT,:,ndx68));
        tails290 = squeeze(Yquantiles2(thisT,:,ndx90));

        %% TWO PANEL PLOT, two models side-by-side
        if doPanel2
            thisfig = figure;

            subplot(1,2,1)
            lhs = gca;
            hold on
            set(gca, 'fontsize', fontsize)
            % plot uncertainty around mid points
            [~,hY] = plotCI([], prctile(squeeze(Ydraws(2:end,thisT,:)), normcdf([-1 1]) * 100, 2), ...
                0:Nhorizons-1, [], 'b-.', 'LineWidth', 2);
            % plot density
            hsim = plot(0:Nhorizons-1,mid, 'r-', 'LineWidth', 4);
            plot(0:Nhorizons-1,tails68, 'r-.', 'LineWidth', 2)
            plot(0:Nhorizons-1,tails90, 'r:', 'LineWidth', 2)

            % plot SPF obs as available
            % thisData = [Zdata(thisT, 1) Yfuture(thisT,1)];
            thisData = Zdata(thisT, 1:2); % connect with nowcast
            hData = plot(-1:0, thisData , 'd-', 'color', [0 .8 0], 'linewidth', 2);
            hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'kd', 'linewidth', 2);
            if ~Znanny(thisT,ndxB)
                hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), 'k:o', 'linewidth', 2);
            end
            if ndxC <= Nz && ~Znanny(thisT,ndxC)
                hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), 'k:o', 'linewidth', 2);
            end
            if ndxD <= Nz && ~Znanny(thisT,ndxD)
                hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), 'k:o', 'linewidth', 2);
            end
            xticks([-1 0 : 2 : Nhorizons-1])
            xlim([-1 Nhorizons-1])


            ylimits = ylim;
            switch upper(datalabel)
                case {'CPI'}
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    elseif all(ylimits < 0)
                        ylimits(2) = 0;
                    end
                    if ~doTight && dates(thisT) >= datenum(2020,1,1)
                        ylimits = [-4 10];
                    end
                case 'UNRATE'
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    elseif all(ylimits < 0)
                        ylimits(2) = 0;
                    end
                case 'RGDP'
                    if ~doTight && dates(thisT) > datenum(2020,1,1)
                        ylimits = [-40 40];
                        %                 else
                        %                     ylimits = [-10 10];
                    end
            end
            ylim(ylimits)
            lhslim = ylim;

            subplot(1,2,2)
            rhs = gca;
            hold on
            set(gca, 'fontsize', fontsize)
            % note: using the mean (to match observed)
            [~,hY] = plotCI([], prctile(squeeze(Ydraws2(2:end,thisT,:)), normcdf([-1 1]) * 100, 2), ...
                0:Nhorizons2-1, [], 'b-.', 'LineWidth', 2);
            % plot density
            hsim = plot(0:Nhorizons2-1,mid2, 'r-', 'LineWidth', 4);
            plot(0:Nhorizons2-1,tails268, 'r-.', 'LineWidth', 2)
            plot(0:Nhorizons2-1,tails290, 'r:', 'LineWidth', 2)

            % plot SPF obs as available
            % thisData = [Zdata(thisT, 1) Yfuture(thisT,1)];
            thisData = Zdata(thisT, 1:2); % connect with nowcast
            hData = plot(-1:0, thisData , 'd-', 'color', [0 .8 0], 'linewidth', 2);
            hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'kd', 'linewidth', 2);
            if ~Znanny(thisT,ndxB)
                hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), 'k:o', 'linewidth', 2);
            end
            if ndxC <= Nz && ~Znanny(thisT,ndxC)
                hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), 'k:o', 'linewidth', 2);
            end
            if ndxD <= Nz && ~Znanny(thisT,ndxD)
                hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), 'k:o', 'linewidth', 2);
            end
            xticks([-1 0 : 2 : Nhorizons-1])
            xlim([-1 Nhorizons-1])


            ylimits = ylim;
            switch upper(datalabel)
                case {'CPI'}
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    elseif all(ylimits < 0)
                        ylimits(2) = 0;
                    end
                    if ~doTight && dates(thisT) >= datenum(2020,1,1)
                        ylimits = [-4 10];
                    end
                case 'UNRATE'
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    elseif all(ylimits < 0)
                        ylimits(2) = 0;
                    end
                case 'RGDP'
                    if ~doTight && dates(thisT) > datenum(2020,1,1)
                        ylimits = [-40 40];
                        %                 else
                        %                     ylimits = [-10 10];
                    end
            end
            ylim(ylimits)
            rhslim = ylim;

            alllims = [min(lhslim(1), rhslim(1)), max(lhslim(2), rhslim(2))];

            ylim(lhs, alllims);
            ylim(rhs, alllims);

            title(lhs, modelpretty)
            title(rhs, modelpretty2)
            sgtitle(sprintf('\\bf %s', datestr(dates(thisT), 'yyyyqq')))

            orient landscape
            wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-%s-vs-%s-2panel', datestr(dates(thisT), 'yyyyqq'), ...
                modellabel, modeltype2), wrap, [], [], [], [], false);
        end

        %% SINGLE PANEL PLOT (omitting mid point distribution)
        SPFcolor  =  Colors4Plots(3);

        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density

        hsim = plotCI(mid, tails, 0:Nhorizons-1, [], 'k-', 'LineWidth', 4);


        hsim2 = plot(0:Nhorizons2-1,mid2, '-.', 'color', Colors4Plots(1), 'LineWidth', 4);
        plot(0:Nhorizons2-1,tails268, '-.', 'color', Colors4Plots(1), 'LineWidth', 2)
        plot(0:Nhorizons2-1,tails290, ':', 'color', Colors4Plots(1), 'LineWidth', 2)

        % plot SPF obs as available
        % thisData = [Zdata(thisT, 1) Yfuture(thisT,1)];
        %         thisData = Zdata(thisT, 1:2); % connect with nowcast
        %         hData = plot(-1:0, thisData , 'd-', 'color', lightblue, 'linewidth', 2);
        hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'd', 'color', SPFcolor, 'linewidth', 2);
        if ~Znanny(thisT,ndxB)
            hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', SPFcolor, 'linewidth', 2);
        end
        if ndxC <= Nz && ~Znanny(thisT,ndxC)
            hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', SPFcolor,  'linewidth', 2);
        end
        if ndxD <= Nz && ~Znanny(thisT,ndxD)
            hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', SPFcolor, 'linewidth', 2);
        end
        avg10last = avg10(thisT);
        if isnan(avg10last)
            thatT = thisT - datesQ(thisT) + 1;
            avg10last = avg10(thatT);
        end

        hAVG10  = plot([0 39], [avg10last avg10last], '-', 'color', SPFcolor, 'linewidth', 2);



        xticks([0 3 : 4 : Nhorizons2-1])
        xlim([0 Nhorizons2-1])


        ylimits = ylim;
        switch upper(datalabel)
            case {'CPI'}
                if all(ylimits > 0)
                    ylimits(1) = 0;
                elseif all(ylimits < 0)
                    ylimits(2) = 0;
                end
                if ~doTight && dates(thisT) >= datenum(2020,1,1)
                    ylimits = [-4 10];
                end
            case 'UNRATE'
                if all(ylimits > 0)
                    ylimits(1) = 0;
                elseif all(ylimits < 0)
                    ylimits(2) = 0;
                end
            case 'RGDP'
                if ~doTight && dates(thisT) > datenum(2020,1,1)
                    ylimits = [-40 40];
                    %                 else
                    %                     ylimits = [-10 10];
                end
        end
        if doTitle
            title(sprintf('\\bf %s', datestr(dates(thisT), 'yyyyqq')))
        end

        orient landscape
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-%s-vs-%s', datestr(dates(thisT), 'yyyyqq'), ...
            modellabel, modeltype2), wrap, [], [], [], [], true);
        hl = legend([hsim hsim2 hSPF hAVG10], modelpretty, modelpretty2, 'SPF', 'SPF-AVG10', 'location', 'best');
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-%s-vs-%s-WITHLEGEND', datestr(dates(thisT), 'yyyyqq'), ...
            modellabel, modeltype2), wrap, [], [], [], [], false);


        %% clean up
        if ispc
            close(thisfig);
        end

    end

    %% finish datalabel loop
    finishwrap


end

%% finish / clean up
finishscript
dockAllFigures
