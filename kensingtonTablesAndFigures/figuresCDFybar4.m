%% collect and tabulate eval stats

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
%#ok<*SAGROW>

%% clear workspace
clear variables

close all
fclose all;
clc


%% parameters

DATALABELS = {'RGDP', 'PGDP', 'UNRATE'};
fontsize   = 18;

Nvariables = length(DATALABELS);

quantileP   = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
percentiles = 1 : 99;
ndxMedian   = percentiles == 50;

datadir     = fullfile('..', 'kensingtonDataMatfiles');

T = 215;
Npercentiles = 99;

fontsize = 18;

doBinEdgeTicks = true;

%% define set of models to compare
Ndraws = 3e3;

m = 0;

m = m + 1;
models(m).type    = 'STATEtrendgapSV';
models(m).ETlabel = 'none';
models(m).pretty  = 'SV';

m = m + 1;
models(m).type    = 'STATEtrendgapSV';
models(m).ETlabel = 'binsOnly';
models(m).pretty  = 'ET (SV)';

m = m + 1;
models(m).type    = 'STATEconst';
models(m).pretty  = 'CONST';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'none';

m = m + 1;
models(m).type    = 'STATEconst';
models(m).pretty  = 'ET (CONST)';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsOnly';

%% set modeldir
resultdir = localresultsMCMC;


%% latexwrapper
wrap   = [];
titlename = sprintf('figuresCDFplots-YBAR4');

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% line colors
SPFcolor   = [0 0.4470 0.7410]; % dark blue

SVcolor    = Colors4Plots(7); % dark red
CONSTcolor = Colors4Plots(3); % dark yellow
RealizedColor = Colors4Plots(5);

%% collect stats for every variable

for d = 1 % : length(DATALABELS)

    %% prepare things


    datalabel = DATALABELS{d};
    bindata   = matfile(fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS',DATALABELS{d})));
    % patch in YBARfuture (if not yet included in mcmc output)
    SPFdata    = matfile(fullfile(datadir,sprintf('kensington%sdata',DATALABELS{d})));
    YBARfuture = SPFdata.YBARfuture;

    ybarQuants = NaN(Npercentiles,T,length(models));
    hh = 1;
    %% load matfiles for each model
    mat = cell(length(models), 1);
    for mm = 1 : length(models)

        modellabel{mm} = models(mm).type;
        ETlabel        = models(mm).ETlabel;

        thismodel = strcat(datalabel, models(mm).type);

        switch ETlabel
            case 'none'
                matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', thismodel, Ndraws);
            case {'binsOnly', 'binsAndMeans'}
                matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel, Ndraws, thismodel);
                modellabel{mm} = strcat(modellabel{mm}, '-ET');
            otherwise
                error('ETlabel <<%s>> not recognized', ETlabel0)
        end
        mat{mm} = matfile(fullfile(resultdir, matfilename));

        thisquant          = permute(mat{mm}.fcstYBARpercentiles, [3 2 1]);
        ybarQuants(:,:,mm) = thisquant(:,hh,:);

        if mm == 1
            dates = mat{mm}.dates;
        else
            if dates ~= mat{mm}.dates
                error('date mismatch')
            end
        end
    end


    %% loop over jump offs
    SPFhistograms = bindata.SPFhistograms;

    for thisT = 156 % [156 178] % 156 % : T

        close all

        SPFbinEdges      = SPFhistograms(thisT).binEdges;
        SPFcdf           = SPFhistograms(thisT).cdf(:,2:end) * 100; % drop current year

        upperSPFbinEdges = SPFhistograms(thisT).upperEdges;
        upperSPFcdf      = SPFhistograms(thisT).upperCDF(:,2:end) * 100; % drop current year
        lowerSPFbinEdges = SPFhistograms(thisT).lowerEdges;
        lowerSPFcdf      = SPFhistograms(thisT).lowerCDF(:,2:end) * 100; % drop current year

        Nhz              = SPFhistograms(thisT).Nhz - 1;

        xtickEdges       = round(SPFbinEdges);

        %% all four in one
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = 1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-', 'color', SVcolor, 'linewidth', 3);
        mm = 2;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', SVcolor,  'linewidth', 3);
        mm = 3;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-',  'color', CONSTcolor, 'linewidth', 3);
        mm = 4;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', CONSTcolor, 'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR4-%s-%s-h%d', datalabel, datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);
        hl = legend([hanni; hSPF], models(1).pretty, models(2).pretty, models(3).pretty, models(4).pretty, ...
            'SPF', 'location', 'northwest');
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR4-%s-%s-h%d-WITHLEGEND', datalabel, datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], false);
        delete(hl)
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        XLIM4 = xlim;
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR4-%s-%s-h%d-WITHDATA', datalabel, datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);
        legend([hanni; hSPF; hData], models(1).pretty, models(2).pretty, models(3).pretty, models(4).pretty, ...
            'SPF', 'Realization', 'location', 'northwest')
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR4-%s-%s-h%d-WITHDATALEGEND', datalabel, datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% SPF only
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-SPF-%s', datalabel, ...
            datestr(dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-SPF-%s-h%d-WITHDATA', datalabel, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);
        legend([hSPF; hData], 'SPF', 'Realization', 'location', 'northwest')
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-SPF-%s-h%d-WITHDATALEGEND', datalabel, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% SV vs CONST
        m1 = 1;
        m2 = 3;
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = m1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-', 'color', SVcolor, 'linewidth', 3);
        %         mm = 2;
        %         hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', SVcolor,  'linewidth', 3);
        mm = m2;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-',  'color', CONSTcolor, 'linewidth', 3);
        %         mm = 4;
        %         hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', CONSTcolor, 'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        legend([hanni([m1 m2]); hSPF], models(m1).pretty, models(m2).pretty, 'SPF', 'location', 'northwest')
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        legend([hanni([m1 m2]); hSPF; hData], models(m1).pretty, models(m2).pretty, 'SPF', 'Realization', 'location', 'northwest')
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s-h%d-WITHDATA', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% SV 
        m1 = 1;
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = m1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-', 'color', SVcolor, 'linewidth', 3);
        %         mm = 2;
        %         hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', SVcolor,  'linewidth', 3);
        %         mm = m2;
        %         hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-',  'color', CONSTcolor, 'linewidth', 3);
        %         mm = 4;
        %         hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', CONSTcolor, 'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        plothorzline(50, [], '--', 'color', SPFcolor, 'linewidth', 2)
        yticks(0 : 25 : 100)
        ylim([0 100])
        legend([hanni(m1); hSPF], models(m1).pretty,'SPF', 'location', 'northwest')
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-%s', datalabel, ...
            modellabel{m1}, ...
            datestr(dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
        % add realized value
        %         hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        %         legend([hanni(m1); hSPF; hData], models(m1).pretty, 'SPF', 'Realization', 'location', 'northwest')
        %         wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-%s-h%d-WITHDATA', datalabel, ...
        %             modellabel{m1}, ...
        %             datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% SV vs SV-ET
        m1 = 1;
        m2 = 2;
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = m1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-', 'color', SVcolor, 'linewidth', 3);
        mm = m2;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', SVcolor,  'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        legend([hanni([m1 m2]); hSPF], models(m1).pretty, models(m2).pretty, 'SPF', 'location', 'northwest')
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        legend([hanni([m1 m2]); hSPF; hData], models(m1).pretty, models(m2).pretty, 'SPF', 'Realization', 'location', 'northwest')
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s-h%d-WITHDATA', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% CONST vs CONST-ET
        m1 = 3;
        m2 = 4;
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = m1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, '-', 'color', CONSTcolor, 'linewidth', 3);
        mm = m2;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', CONSTcolor,  'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        legend([hanni([m1 m2]); hSPF], models(m1).pretty, models(m2).pretty, 'SPF', 'location', 'northwest')
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ')), wrap, [], [], [], [], true);
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        legend([hanni([m1 m2]); hSPF; hData], models(m1).pretty, models(m2).pretty, 'SPF', 'Realization', 'location', 'northwest')
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s-h%d-WITHDATA', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

        %% SV-ET vs CONST-ET
        m1 = 2;
        m2 = 4;
        hanni = NaN(length(models),1);
        thisfig = figure;
        set(gca, 'FontSize', fontsize)
        hold on
        mm = m1;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':', 'color', SVcolor, 'linewidth', 3);
        mm = m2;
        hanni(mm) = plot(ybarQuants(:,thisT,mm), percentiles, ':',  'color', CONSTcolor, 'linewidth', 3);
        if hh <= Nhz
            hSPF = plot(SPFbinEdges, SPFcdf(:,hh), 'd', 'color', SPFcolor, 'linewidth', 2);
            plot(upperSPFbinEdges, upperSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
            plot(lowerSPFbinEdges, lowerSPFcdf(:,hh), '-', 'color', SPFcolor, 'linewidth', 1)
        end
        yticks(0 : 25 : 100)
        ylim([0 100])
        legend([hanni([m1 m2]); hSPF], models(m1).pretty, models(m2).pretty, 'SPF', 'location', 'northwest')
        if doBinEdgeTicks
            xticks(xtickEdges);
        end
        grid on
        orient landscape
        xlim(XLIM4)
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s-h%d', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], false);
        % add realized value
        hData = plot([YBARfuture(thisT,1) YBARfuture(thisT,1)], [0 100], '-', 'color', RealizedColor, 'linewidth', 2);
        legend([hanni([m1 m2]); hSPF; hData], models(m1).pretty, models(m2).pretty, 'SPF', 'Realization', 'location', 'northwest')
        wrapthisfigure(thisfig, sprintf('cdfplot-YBAR-%s-%s-vs-%s-%s-h%d-WITHDATA', datalabel, ...
            modellabel{m1}, modellabel{m2}, ...
            datestr(dates(thisT), 'yyyyQQ'), hh), wrap, [], [], [], [], true);

    end % thisT
end % DATALABEL


%% finish / clean up
finishwrap
finishscript
dockAllFigures




