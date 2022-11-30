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

datadir = fullfile('..', 'kensingtonDataMatfiles');

% set modeldir
resultdir = localresultsMCMC;

%% define set of models to compare against SPF

m = 1;
models(m).type   = 'STATEtrendgapSV';
models(m).pretty = 'SV';
models(m).Ndraws = 3e3;

m = m + 1;
models(m).type   = 'STATEconst';
models(m).pretty = 'CONST';
models(m).Ndraws = 3e3;


%% define set of eval windows
s = 0;

% s = s + 1;
% samples(s).start  = [];
% samples(s).stop   = [];
% samples(s).label  = 'fullsample';
% samples(s).pretty = 'full sample';

s = s + 1;
samples(s).start  = datenum(1992,1,1);
samples(s).stop   = [];
samples(s).label  = 'since1992';
samples(s).pretty = 'since 1992';

s = s + 1;
samples(s).start  = datenum(1992,1,1);
samples(s).stop   = datenum(2016,12,1);
samples(s).label  = 'preCOVIDsince1992';
samples(s).pretty = 'pre COVID (1992-2016)';

s = s + 1;
samples(s).start  = datenum(2009,1,1);
samples(s).stop   = [];
samples(s).label  = 'since2009';
samples(s).pretty = 'since 2009';

% s = s + 1;
% samples(s).start  = [];
% samples(s).stop   = datenum(2016,12,1);
% samples(s).label  = 'preCOVID';
% samples(s).pretty = 'pre COVID (2016)';

s = s + 1;
samples(s).start  = datenum(2009,1,1);
samples(s).stop   = datenum(2016,12,1);
samples(s).label  = 'preCOVIDsince2009';
samples(s).pretty = 'pre COVID (2009-2016)';



%% latexwrapper
wrap   = [];

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% allocate memory
maxNbar = 3;
T       = 215;

[ DRPS0, DRPS1, dmTstat, firstObs, lastObs ] = deal(NaN(maxNbar, length(samples), length(DATALABELS), length(models)));
[ LOSS0, LOSS1, AVGLOSS0, AVGLOSS1 ] = deal(NaN(T, maxNbar, length(samples), length(DATALABELS), length(models)));


%% loop over models

for mm = 1 : length(models)

    modeltype   = models(mm).type;
    modelpretty = models(mm).pretty;
    Ndraws      = models(mm).Ndraws;

    %% collect stats for every variable

    for dd = 1 : length(DATALABELS)

        %% prepare things


        datalabel = DATALABELS{dd};
        fprintf('Processing %s ... \n', datalabel)


        %% load SPF data
        PROBdata = matfile(fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS',DATALABELS{dd})));


        %% load MCMC data
        modellabel0  = strcat(datalabel, modeltype);

        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel0, Ndraws);
        mat0        = matfile(fullfile(resultdir, matfilename));
        dates       = mat0.dates;
        theseNbar   = mat0.Nbar;
        Tstart      = mat0.Tstart;
        T           = mat0.T;

        for ss = 1 : length(samples)

            if isempty(samples(ss).start)
                samStart = Tstart;
            else
                samStart = find(dates >= samples(ss).start, 1, 'first');
            end

            if isempty(samples(ss).stop)
                samStop = T;
            else
                samStop = find(dates <= samples(ss).stop, 1, 'last');
            end

            sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));

            comparisonNote = sprintf('Evaluation window from %s through %s (and as far as realized values are available).', ...
                datestr(dates(samStart), 'yyyyqq'), datestr(dates(samStop), 'yyyyqq'));


            loss0 = PROBdata.fcstSPFdrps;
            loss1 = mat0.fcstYBARdrps;

            % establish common sample (needed for reporting loss levels; dmtest would catch it)
            nanny        = (isnan(loss0) | isnan(loss1)) | ~sammy; % note: automatic array expansion w.r.t. sammy
            loss0(nanny) = NaN;
            loss1(nanny) = NaN;

            for hh = 1 : theseNbar
                if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                    [~, dmTstat(hh,ss,dd,mm)] = dmtest(loss0(:,hh), loss1(:,hh), hh * 4);
                else
                    dmTstat(hh,ss,dd,mm) = 0;
                end
            end

            DRPS0(1:theseNbar, ss, dd, mm) = mean(loss0,1, 'omitnan');
            DRPS1(1:theseNbar, ss, dd, mm) = mean(loss1,1, 'omitnan');

            for hh = 1 : theseNbar
                thesedates                = dates(~nanny(:,hh));
                firstObs(hh, ss, dd, mm)  = min(thesedates);
                lastObs(hh, ss, dd, mm)   = max(thesedates);
            end

            % collect losses (for comparison plots)
            LOSS0(:, 1:theseNbar, ss, dd, mm) = loss0;
            LOSS1(:, 1:theseNbar, ss, dd, mm) = loss1;

            for hh = 1 : theseNbar
                nanny           = isnan(loss0(samStart:samStop, hh));
                Tndx            = cumsum(~nanny);
                avgloss0        = cumsum(loss0(samStart:samStop, hh), 'omitnan') ./ Tndx;
                avgloss0(nanny) = NaN;

                if ~isequal(nanny, isnan(loss1(samStart:samStop, hh)))
                    error('sample mismatch') % the DRPS calculations above should have enforced a common sample
                end
                avgloss1        = cumsum(loss1(samStart:samStop, hh), 'omitnan') ./ Tndx;
                avgloss1(nanny) = NaN;
                AVGLOSS0(samStart:samStop, hh, ss, dd, mm) = avgloss0;
                AVGLOSS1(samStart:samStop, hh, ss, dd, mm) = avgloss1;
            end
        end % sample
    end % datalabel
end % models

relDRPS = DRPS1 ./ DRPS0;

% obs should be the same across models
checkdiff(range(firstObs,4));
checkdiff(range(lastObs,4));
firstObs = firstObs(:,:,:,1);
lastObs  = lastObs(:,:,:,1);

% LOSS0 (SPF) should be same across models
checkdiff(range(LOSS0, 5));
checkdiff(range(AVGLOSS0, 5));
LOSS    = cat(5, LOSS0(:,:,:,:,1), LOSS1);
AVGLOSS = cat(5, AVGLOSS0(:,:,:,:,1), AVGLOSS1);

%% tabulate
panellabel = {'A', 'B', 'C'};

for dd = 1 : length(DATALABELS)
    datalabel = DATALABELS{dd};

    tabname    = sprintf('tableDRPS-%s.tex', datalabel);
    tabcaption = tabname;


    % set up tab
    if ~isempty(wrap)
        tabdir = wrap.dir;
        latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
        %      latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
    else
        tabdir = fullfile(localtemp, 'foo');
    end


    % tabulate
    Ncols = 1 + length(samples);
    fid = fopen(fullfile(tabdir, tabname), 'wt');
    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.5', 1, length(samples)));
    fprintf(fid, '\\toprule\n');
   
    for hh = 1 : maxNbar
        switch hh
            case 1
                calyearlabel = 'Next year';
            case 2
                calyearlabel = 'Two years ahead';
            case 3
                calyearlabel = 'Three years ahead';
            otherwise
                error('cal year %d not implemented', hh)
        end

        if hh == 1
            fprintf(fid, '\\multicolumn{%d}{c}{\\textbf{PANEL %s: %s}} ', Ncols, panellabel{hh}, calyearlabel);
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\midrule\n');
            fprintf(fid, ' & \\multicolumn{%d}{c}{Samples} ', length(samples));
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\cmidrule(lr){%d-%d}\n',2,Ncols);
            fprintf(fid, 'Model ');
            for ss = 1 : length(samples)
                fprintf(fid, ' & \\multicolumn{1}{l}{%s} ', sprintf('%02d-%02d', year(firstObs(hh,ss,dd)), year(lastObs(hh,ss,dd))));
            end
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\midrule\n');

            for mm = 1 : length(models)
                fprintf(fid, '%s ', models(mm).pretty);
                for ss   = 1 : length(samples)
                    this = relDRPS(hh,ss,dd,mm);
                    that = dmTstat(hh,ss,dd,mm);
                    if isnan(this)
                        fprintf(fid, '& \\multicolumn{1}{c}{--} ');
                    else
                        fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
                    end
                end
                fprintf(fid, '\\\\\n');
            end
        else
            fprintf(fid, '\\midrule\n');
            fprintf(fid, '\\multicolumn{%d}{c}{\\textbf{PANEL %s: %s}} ', Ncols, panellabel{hh}, calyearlabel);
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\midrule\n');
            fprintf(fid, ' & & & \\multicolumn{%d}{c}{Samples} ', 2);
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\cmidrule(lr){%d-%d}\n',4,Ncols);
            fprintf(fid, 'Model ');
            for ss = 1 : 2
                %                 fprintf(fid, ' & \\multicolumn{1}{l}{%s} ', '--');
                fprintf(fid, ' & ');
            end
            for ss = 3 : length(samples)
                fprintf(fid, ' & \\multicolumn{1}{l}{%s} ', sprintf('%02d-%02d', year(firstObs(hh,ss,dd)), year(lastObs(hh,ss,dd))));
            end
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\midrule\n');

            for mm = 1 : length(models)
                fprintf(fid, '%s ', models(mm).pretty);
                for ss = 1 : 2
                    fprintf(fid, '& \\multicolumn{1}{c}{--} ');
                end
                for ss   = 3 : length(samples)
                    this = relDRPS(hh,ss,dd,mm);
                    that = dmTstat(hh,ss,dd,mm);
                    if isnan(this)
                        fprintf(fid, '& \\multicolumn{1}{c}{--} ');
                    else
                        fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
                    end
                end
                fprintf(fid, '\\\\\n');
            end

        end % hh
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, '\n');
    %
    %
    fprintf(fid, '\\begin{small}\n');
    fprintf(fid, 'Note: \n');



    fprintf(fid, 'DRPS of model-based calendar-year forecasts for %s relative to SPF histograms (and using SPF bins) with SPF in denominator.\n', datalabel);
    fprintf(fid, 'Quarterly forecast observations.\n');
    fprintf(fid, 'Samples listed at top of each panel refer to first and last year of observations with available SPF forecast and subsequently realized value.\n');
    fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $y\\cdot 4$ quarterly lags, where $y$ denotes the value of the annual forecast horizons.\n');
    fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
    % fprintf(fid, '``%s'''' and ``%s'''' refer to subsamples extending from %s until %s and %s, respectively (using data for realized values as far as available in %s).', ...
    %     sam(1).pretty, sam(2).pretty, datestr(sam(1).start, 'yyyyQQ'), datestr(sam(1).stop, 'yyyyQQ'), datestr(sam(2).stop, 'yyyyQQ'), datestr(dates(end), 'yyyyQQ'));
    fprintf(fid, '\\end{small}\n');
    %
    fclose(fid);
    type(fullfile(tabdir, tabname))
end % dd

%% plot DRPS (longest sample)

ss = 1;
if isempty(samples(ss).start)
    samStart = Tstart;
else
    samStart = find(dates >= samples(ss).start, 1, 'first');
end

if isempty(samples(ss).stop)
    samStop = T;
else
    samStop = find(dates <= samples(ss).stop, 1, 'last');
end
thesedates = dates(samStart:samStop);

Nmodels = size(LOSS,5);
modellist = {'SPF', models.pretty};
linestyles = {'-', '-.', ':'};
colorlist  = {Colors4Plots(8), Colors4Plots(7), Colors4Plots(3)};
linewidths = {5, 3, 3};

for dd = 1 : length(DATALABELS)
    datalabel = DATALABELS{dd};

    hanni = NaN(Nmodels, 1);
    % DRPS contributions
    for hh = 1 : maxNbar
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        for mm = 1 : Nmodels
            hanni(mm) = plot(thesedates, LOSS(samStart:samStop, hh, ss, dd, mm), ...
                'linestyle', linestyles{mm}, 'color', colorlist{mm}, 'LineWidth', linewidths{mm});
        end        
        ylim([0 max(ylim)])
        xtickdates(thesedates)
        orient landscape
        wrapthisfigure(thisfig, sprintf('DRPS-%s-y%d', datalabel, hh), wrap, [], [], [], [], true)
        legend(hanni, modellist, 'location', 'northwest', 'fontsize', 22)
        wrapthisfigure(thisfig, sprintf('DRPS-%s-y%d-WITHLEGEND', datalabel, hh), wrap)
    end
    % growing averages
    for hh = 1 : maxNbar
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        for mm = 1 : Nmodels
            hanni(mm) = plot(thesedates, AVGLOSS(samStart:samStop, hh, ss, dd, mm), ...
                'linestyle', linestyles{mm}, 'color', colorlist{mm}, 'LineWidth', linewidths{mm});
        end        
        %         ylim([0 max(ylim)])
        xtickdates(thesedates)
        orient landscape
        wrapthisfigure(thisfig, sprintf('avgDRPS-%s-y%d', datalabel, hh), wrap, [], [], [], [], true)
        legend(hanni, modellist, 'location', 'best', 'fontsize', 22)
        wrapthisfigure(thisfig, sprintf('avgDRPS-%s-y%d-WITHLEGEND', datalabel, hh), wrap)
    end

     % growing averages as difference to SPF
    for hh = 1 : maxNbar
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        for mm = 2 : Nmodels
            diffAVGDRPS = AVGLOSS(samStart:samStop, hh, ss, dd, mm) - AVGLOSS(samStart:samStop, hh, ss, dd, 1);
            hanni(mm)  = plot(thesedates, diffAVGDRPS, ...
                'linestyle', linestyles{mm}, 'color', colorlist{mm}, 'LineWidth', linewidths{mm});
        end        
        YLIM = ylim;
        ylim([min(YLIM(1), 0), max(YLIM(2), 0)])
        xtickdates(thesedates)
        plotOrigin
        orient landscape
        wrapthisfigure(thisfig, sprintf('diffavgDRPS-%s-y%d', datalabel, hh), wrap, [], [], [], [], true)
        legend(hanni(2:end), modellist(2:end), 'location', 'best', 'fontsize', 22)
        wrapthisfigure(thisfig, sprintf('diffavgDRPS-%s-y%d-WITHLEGEND', datalabel, hh), wrap)
    end
end

%% output time series
ss = 1;
for dd = 1 : length(DATALABELS)
    datalabel = DATALABELS{dd};
    for hh = 1 : maxNbar
        % LOSS
        outx    = squeeze(LOSS(:,hh,ss,dd,:));
        csvname = sprintf('DRPS-%s-y%d', datalabel, hh);
        writedatatable(wrap, csvname, dates, outx, modellist, 'yyyyQQ');
        % AVGLOSS
        outx    = squeeze(AVGLOSS(:,hh,ss,dd,:));
        csvname = sprintf('avgDRPS-%s-y%d', datalabel, hh);
        writedatatable(wrap, csvname, dates, outx, modellist, 'yyyyQQ');
    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures


