%% plot width of SPF bins

clear; close all; clc;

%#ok<*UNRCH>

%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/

%% parameters
DATALABELS = {'RGDP','PGDP','UNRATE'};

wrap = [];
initwrap;

for doShowPre1992 = [true false]


    %% bin endpoints
    % PGDP
    PRPGDP1_bins=[-3,-2.05:1:9.95]; % 1968Q4-1973Q1
    PRPGDP2_bins=[-1,-0.05:1:11.95]; % 1973Q2-1974Q3
    PRPGDP3_bins=[3,3.95:1:15.95]; % 1974Q4-1981Q2
    PRPGDP4_bins=[4,5.95:2:11.95]; % 1981Q3-1985Q1  next year got included
    PRPGDP5_bins=[2,3.95:2:9.95]; % 1985Q2-1991Q4
    PRPGDP6_bins=[0,0.95:1:7.95]; % 1992Q1-2013Q4
    PRPGDP7_bins=[0,0.45:0.5:3.95]; % 2014Q1-present

    bins.PGDP = cat(1, {PRPGDP1_bins}, {PRPGDP2_bins}, {PRPGDP3_bins}, {PRPGDP4_bins}, ...
        {PRPGDP5_bins}, {PRPGDP6_bins}, {PRPGDP7_bins});
    binDates.PGDP = datetime({'1968Q4', '1973Q2', '1974Q4', '1981Q3', ...
        '1985Q2', '1992Q1', '2014Q1'}, 'InputFormat', 'yyyyQQQ');


    % RGDP
    PRGDP1_bins=[-3,-2.05:1:9.95]; % 1968Q4-1973Q1
    PRGDP2_bins=[-1,-0.05:1:11.95]; % 1973Q2-1974Q3
    PRGDP3_bins=[3,3.95:1:15.95]; % 1974Q4-1981Q2
    PRGDP4_bins=[-2,-0.05:2:5.95]; % 1981Q3-1991Q4  next year got included
    PRGDP5_bins=[-2,-1.05:1:5.95]; % 1992Q1-2009Q1
    PRGDP6_bins=[-3,-2.05:1:5.95]; % 2009Q2-2020Q1  years 3 and 4 got included
    PRGDP7_bins=[-12,-6.05,-3.05,-0.5,1.45,2.45,3.95,6.95,9.95,15.95]; % 2020Q2-present

    bins.RGDP = cat(1, {PRGDP1_bins}, {PRGDP2_bins}, {PRGDP3_bins}, {PRGDP4_bins}, ...
        {PRGDP5_bins}, {PRGDP6_bins}, {PRGDP7_bins});
    binDates.RGDP = datetime({'1968Q4', '1973Q2', '1974Q4', '1981Q3', ...
        '1992Q1', '2009Q2', '2020Q1'}, 'InputFormat', 'yyyyQQQ');


    % UNEMP
    PRUNEMP1_bins=[6,6.95,7.45,7.95,8.45,8.95,9.45,9.95,10.95]; % 2009Q2-2013Q4     next, years 3 and 4 included
    PRUNEMP2_bins=[4,4.95,5.45,5.95,6.45,6.95,7.45,7.95,8.95]; % 2014Q1-2020Q1
    PRUNEMP3_bins=[3,3.95,4.95,5.95,6.95,7.95,9.95,11.95,14.95]; % 2020Q2-present

    bins.UNRATE     = cat(1, {PRUNEMP1_bins}, {PRUNEMP2_bins}, {PRUNEMP3_bins});
    binDates.UNRATE = datetime({'2009Q2', '2014Q1', '2020Q1'}, 'InputFormat', 'yyyyQQQ');


    %% prune bin dates prior to 1981Q3 (which pertain only to the current year)

    plotLabel = '';

    if doShowPre1992
        plotLabel       = strcat(plotLabel, '-pre1992');

        % copy to separate array (to control plot colors)
        binDatesPre1992 = binDates;
        binsPre1992     = bins;

        for d = 1 : length(DATALABELS)
            this = DATALABELS{d};
            ndx92 = binDatesPre1992.(this) < datetime('1992Q1', 'InputFormat', 'yyyyQQQ');
            ndx81 = binDatesPre1992.(this) >= datetime('1981Q3', 'InputFormat', 'yyyyQQQ');
            ndx = ndx81 & ndx92;
            binDatesPre1992.(this) = binDatesPre1992.(this)(ndx);
            binsPre1992.(this)     = binsPre1992.(this)(ndx);

        end
    end

    % whack out the pre 1992 bins
    for d = 1 : length(DATALABELS)
        this = DATALABELS{d};
        ndx  = binDates.(this) < datetime('1992Q1', 'InputFormat', 'yyyyQQQ');

        binDates.(this) = binDates.(this)(~ndx);
        bins.(this)     = bins.(this)(~ndx);

    end

    %% plot
    for d = 1 % : length(DATALABELS)
        this = DATALABELS{d};
        Ndates = length(binDates.(this));

        matpoints = matfile(fullfile('..', 'kensingtonDataMatfiles', ...
            sprintf('kensington%sdata', this)));
        ndx            = [1 6 matpoints.Nforecasts];
        ndx            = unique(ndx); % Nforecasts can be equal to 6
        Nforecasts     = length(ndx);
        Yforecasts     = matpoints.YactualForecasts;
        Yforecasts     = Yforecasts(:,ndx);
        YforecastLabel = matpoints.YforecastLabel;
        YforecastLabel = YforecastLabel(ndx,1);

        Ydates         = datetime(matpoints.dates, 'convertfrom', 'datenum');
        YforecastLabel = cellfun(@(x) sprintf('Point forecast %s', x), ...
            YforecastLabel, 'Uniformoutput', false);


        %% PLOT WITH BINS
        thisfig = figure;
        set(gca, 'fontsize', 18)
        hold on

        % plot point forecasts
        greyscales = [.5 .25 0];
        theseLines = {'-', '-.', '--'};
        hmid = NaN(Nforecasts, 1);
        for nn = 1 % : Nforecasts
            hmid(nn) = plot(Ydates, Yforecasts(:,nn), theseLines{nn}, 'color', greyscales(nn) * [1 1 1], 'linewidth', 2);
        end

        % plot diamonds for when new bin widths were introduced
        for nn = 1 : Ndates
            Nbins = length(bins.(this){nn});
            hbin = plot(repmat(binDates.(this)(nn), Nbins, 1), bins.(this){nn}, 'd', 'linewidth', 2, 'color', Colors4Plots(nn));
        end

        % y axis limits
        YLIM = ylim;
        if YLIM(1) > 0
            YLIM(1) = 0;
        end
        if isodd(YLIM(2))
            YLIM(2) = YLIM(2) + 1;
        else
            YLIM(2) = YLIM(2) + 2;
        end
        ylim(YLIM)

        % x axis limits
        XLIM = xlim;
        xlim([dateshift(XLIM(1), 'start', 'year'), datetime('2022-12-1')])
        XLIM = xlim;

        % vertical lines for years when bin edges introduced
        for nn = 1 : Ndates
            that = plot(repmat(binDates.(this)(nn), 2, 1), YLIM, '--', 'color', Colors4Plots(nn), 'linewidth', 1);
            if nn == 1
                hbinyear = that;
            end
        end

        % vertical lines for years when bin edges introduced
        for nn = 1 : Ndates
            that = plot(repmat(binDates.(this)(nn), 2, 1), YLIM, '--', 'color', Colors4Plots(nn), 'linewidth', 1);
            if nn == 1
                hbinyear = that;
            end
        end
        % plot horizontal lines for bin edges
        for nn = 1 : Ndates - 1
            that = plot(binDates.(this)([nn nn+1]), [bins.(this){nn}; bins.(this){nn}], ':', 'color', Colors4Plots(nn), 'linewidth', 2);
            if nn == 1
                hbin = that;
            end
        end
        nn = Ndates;
        plot([binDates.(this)(nn) XLIM(2)], [bins.(this){nn}; bins.(this){nn}], ':', 'color', Colors4Plots(nn), 'linewidth', 2);

        if doShowPre1992 && ~isempty(binDatesPre1992.(this))
            for nn = 1 : length(binDatesPre1992.(this))
                % plot diamonds for when new bin widths were introduced
                thiscolor = Colors4Plots(5); % Colors4Plots(Ndates + nn);
                Nbins = length(binsPre1992.(this){nn});
                plot(repmat(binDatesPre1992.(this)(nn), Nbins, 1), binsPre1992.(this){nn}, 'd', 'linewidth', 2, 'color', thiscolor);
                % vertical lines for years when bin edges introduced
                plot(repmat(binDatesPre1992.(this)(nn), 2, 1), YLIM, '--', 'color', thiscolor, 'linewidth', 1);
            end

            % plot horizontal lines for bin edges
            for nn = 1 : length(binDatesPre1992.(this)) - 1
                thiscolor = Colors4Plots(5); Colors4Plots(Ndates + nn);
                plot(binDatesPre1992.(this)([nn nn+1]), [binsPre1992.(this){nn}; binsPre1992.(this){nn}], ':', 'color', thiscolor, 'linewidth', 2);
            end
            nn = length(binDatesPre1992.(this));
            plot([binDatesPre1992.(this)(nn) binDates.(this)(1)], [binsPre1992.(this){nn}; binsPre1992.(this){nn}], ':', 'color', thiscolor, 'linewidth', 2);
            % shade pre 1992 bins
            %             shadedates = [binDatesPre1992.(this)(1) binDates.(this)(1)];
            %             YLIM = ylim;
            %             shadecolor = .7 * [1 1 1];
            %             hanni1 = bar(shadedates, [min(YLIM) min(YLIM)], 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);
            %             hanni2 = bar(shadedates, [max(YLIM) max(YLIM)], 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);
            %             uistack([hanni1, hanni2], 'bottom')
            %             set(gca, 'layer', 'top')
        end

        plothorzline(0, [], 'k-')
        %         if d == 1
        %             legendLocation = 'southwest';
        %         else
        %             legendLocation = 'best';
        %         end
        wrapthisfigure(thisfig,sprintf('binwidths-%s%s', this, plotLabel), wrap, [], [], [], [], true)
        %         hl = legend(hmid, YforecastLabel, 'Location', legendLocation);
        %         wrapthisfigure(thisfig,sprintf('binwidths-%s%s-WITHLEGEND', this, plotLabel), wrap, [], [], [], [], true)
        title(this)
        %         wrapthisfigure(thisfig,sprintf('binwidths-%s%s-WITHTITLELEGEND', this, plotLabel), wrap)
        %         delete(hl)
        wrapthisfigure(thisfig,sprintf('binwidths-%s%s-WITHTITLE', this, plotLabel), wrap, [], [], [], [], true)

        %% REDO PLOT WITHOUT BINS
        thisfig = figure;
        set(gca, 'fontsize', 18)
        hold on

        % plot point forecasts
        greyscales = [.5 .25 0];
        theseLines = {'-', '-.', '--'};
        hmid = NaN(Nforecasts, 1);
        nn = 1;
        hmid(nn) = plot(Ydates, Yforecasts(:,nn), theseLines{nn}, 'color', greyscales(nn) * [1 1 1], 'linewidth', 2);
        for nn = 2 : Nforecasts
            hmid(nn) = plot(Ydates, Yforecasts(:,nn), theseLines{nn}, 'color', Colors4Plots(nn), 'linewidth', 2);
        end

        % plot vertical lines for first available forecast
        for nn = 2 : Nforecasts
            thisNdx = find(~isnan(Yforecasts(:,nn)), 1, 'first');
            plot(Ydates([thisNdx thisNdx]), YLIM, theseLines{nn}, 'color', Colors4Plots(nn), 'linewidth', 2);
            fprintf('%s\n', datestr(Ydates(thisNdx)))
        end

        if d == 1 
            YforecastLabel = {'Quarterly forecasts since 1968Q4', 'Next-year since 1981Q3', 'Two- and three-years since 2009Q2'};
        end

        % y axis limits sam as with bins
        ylim(YLIM)

        % x axis limits
        XLIM = xlim;
        xlim([dateshift(XLIM(1), 'start', 'year'), datetime('2022-12-1')])
        XLIM = xlim;

        plothorzline(0, [], 'k-')
        if d == 1
            legendLocation = 'southwest';
        else
            legendLocation = 'best';
        end
        wrapthisfigure(thisfig,sprintf('pointforecasts-%s%s', this, plotLabel), wrap, [], [], [], [], true)
        hl = legend(hmid, YforecastLabel, 'Location', legendLocation);
        wrapthisfigure(thisfig,sprintf('pointforecasts-%s%s-WITHLEGEND', this, plotLabel), wrap, [], [], [], [], true)
        title(this)
        wrapthisfigure(thisfig,sprintf('pointforecasts-%s%s-WITHTITLELEGEND', this, plotLabel), wrap)
        delete(hl)
        wrapthisfigure(thisfig,sprintf('pointforecasts-%s%s-WITHTITLE', this, plotLabel), wrap, [], [], [], [], true)
    end

end

%% finish
dockAllFigures
finishwrap
finishscript