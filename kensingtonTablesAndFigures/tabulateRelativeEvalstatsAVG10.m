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

DATALABELS = {'RGDP', 'CPI', 'TBILL'};
fontsize   = 18;

Nhorizons  = 16;
Nvariables = length(DATALABELS);

doPlotsRelStats = false;
doPlotsCumLoss  = false;

doMSE = false; % to report relative MSE vs RMSE

%% define set of models to compare

PAIRS = [];

PAIRS = {[1 2]}; 

m = 1;
models(m).type   = 'STATEtrendgapSV';
models(m).pretty = 'SV';
models(m).Ndraws = 3e3;

m = m + 1;
models(m).type   = 'STATEtrendgaplongSV';
models(m).pretty = 'SV-AVG10';
models(m).Ndraws = 3e3;



%% set modeldir
resultdir = localresultsMCMC;

%% define set of eval windows
s = 1;
sam(s).start  = [];
sam(s).stop   = [];
sam(s).label  = 'fullsample';
sam(s).pretty = 'full sample';

s = 2;
sam(s).start  = [];
sam(s).stop   = datenum(2016,12,1);
sam(s).label  = 'preCOVID';
sam(s).pretty = 'pre COVID (2016)';

s = 3;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = [];
sam(s).label  = 'fullsampleSince1992';
sam(s).pretty = 'full sample (since 1992)';

s = 4;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = datenum(2016,12,1);
sam(s).label  = 'preCOVIDsince1992';
sam(s).pretty = 'pre COVID (1992-2016)';


%% define model pairs to consider

if isempty(PAIRS)
    PAIRS = arrayfun(@(x) [1 x], 2:length(models), 'UniformOutput', false);
end

%% latexwrapper
wrap   = [];
titlename = 'tabulateRelativeEvalStats';
if doPlotsCumLoss
    titlename = strcat(titlename, '-withCumLossPlots');
end
if doPlotsRelStats
    titlename = strcat(titlename, '-withRelStatPlots');
end
initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end



%% loop over model pairs

for pp = 1 : length(PAIRS)

    close all
    thispair = PAIRS{pp};
    m0       = thispair(1);
    m1       = thispair(2);

    modeltype0   = models(m0).type;
    modelpretty0 = models(m0).pretty;
    Ndraws0      = models(m0).Ndraws;

    modeltype1   = models(m1).type;
    modelpretty1 = models(m1).pretty;
    Ndraws1      = models(m1).Ndraws;


    for s = 1  : length(sam)

        %% allocate memory for forecast stats
        [MAE0, MAE1, RMSE0, RMSE1, CRPS0, CRPS1]  = deal(NaN(Nhorizons, Nvariables));
        [relRMSEtstat, relMAEtstat, relCRPStstat] = deal(NaN(Nhorizons, Nvariables));

        MAXHORIZON = NaN(length(DATALABELS), 1);

        %% collect stats for every variable

        for d = 1 : length(DATALABELS)

            %% prepare things


            datalabel = DATALABELS{d};

            modellabel0  = strcat(datalabel, modeltype0);
            modellabel1 = strcat(datalabel, modeltype1);


            %% load data

            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel0, Ndraws0);
            mat0 = matfile(fullfile(resultdir, matfilename));

            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel1, Ndraws1);
            mat1 = matfile(fullfile(resultdir, matfilename));


            % YdensityDraws  = mat.YdensityDraws;
            % YdensityDraws1 = mat1.YdensityDraws;

            Yfuture       = mat0.Yfuture;
            Nz            = mat0.Nz;
            Zdata         = mat0.Zdata;
            Znanny        = mat0.Znanny;
            Ylabel        = mat0.Ylabel;
            Yfcstlabel    = Ylabel(2:end);
            dates         = mat0.dates;
            T             = mat0.T;
            Tstart        = mat0.Tstart;

            MAXHORIZON(d)  = mat0.Nhorizons;
            theseNhorizons = mat0.Nhorizons;

            % checks
            if ~isequal(dates, mat1.dates)
                error('mismatch between input files for model1 and model2 (dates)')
            end
            if ~isequal(Tstart, mat1.Tstart)
                error('mismatch between input files for model1 and model2 (Tstart)')
            end
            if ~isequaln(Zdata, mat1.Zdata)
                warning('mismatch between input files for model1 and model2 (Zdata, %s)', datalabel)
            end
            if ~isequaln(Yfuture, mat1.Yfuture)
                warning('mismatch between input files for model1 and model2 (Yfuture)')
            end
            if Nhorizons < mat1.Nhorizons
                warning('Nhorizons does not match assumption')
            end
            if mat0.Nhorizons ~= mat1.Nhorizons
                warning('mismatch between input files (Nhorizons)')
            end

            fprintf('Processing %s ... \n', datalabel)


            if isempty(sam(s).start)
                samStart = Tstart;
            else
                samStart = find(dates >= sam(s).start, 1, 'first');
            end

            if isempty(sam(s).stop)
                samStop = T;
            else
                samStop = find(dates <= sam(s).stop, 1, 'last');
            end

            sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));

            comparisonNote = sprintf('Evaluation window from %s through %s (and as far as realized values are available).', ...
                datestr(dates(samStart), 'yyyyqq'), datestr(dates(samStop), 'yyyyqq'));


            %% RMSE
            if doMSE
                thistab = 'MSE';
            else
                thistab = 'RMSE';
            end
            tabname = sprintf('%s-%s.tex', thistab, datalabel);

            loss0 = (mat0.fcstYhatRBerror).^2;
            loss1 = (mat1.fcstYhatRBerror).^2;
            loss1 = loss1(:,1:theseNhorizons);

            % establish common sample (needed for reporting loss levels; dmtest would catch it)
            nanny        = (isnan(loss0) | isnan(loss1)) | ~sammy; % note: automatic array expansion w.r.t. sammy
            loss0(nanny) = NaN;
            loss1(nanny) = NaN;

            for hh = 1 : theseNhorizons
                if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                    [~, relRMSEtstat(hh,d)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
                else
                    relRMSEtstat(hh,d) = 0;
                    %                     warning('virtually identical means (h=%d, %s)', hh, datalabel)
                end
            end

            if doMSE
                RMSE0(1:theseNhorizons,d) = mean(loss0,1, 'omitnan');
                RMSE1(1:theseNhorizons,d) = mean(loss1,1, 'omitnan');
            else
                RMSE0(1:theseNhorizons,d) = sqrt(mean(loss0,1, 'omitnan'));
                RMSE1(1:theseNhorizons,d) = sqrt(mean(loss1,1, 'omitnan'));
            end

            %% plot cumulative losses
            if doPlotsCumLoss && s == 1

                cumlossdiff        = cumsum(loss1 - loss0, 1, 'omitnan');
                nanny              = isnan(loss1 - loss0);
                avgcumlossdiff     = cumlossdiff ./ cumsum(~nanny);
           
                cumlossdiff(nanny)    = NaN;
                avgcumlossdiff(nanny) = NaN;

                grps = {0:4,5:8,9:12,13:15};

                for gg = 1 : length(grps)
                    thisgroup = grps{gg} + 1;
                    thisgroup = thisgroup(thisgroup <= theseNhorizons);
                    if ~isempty(thisgroup)
                        thisfig = figure;
                        subplot(2,1,1)
                        h = plot(dates(sammy), cumlossdiff(sammy,thisgroup), 'linewidth', 1);
                        set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
                        %                         legend(h, Yfcstlabel(thisgroup), 'Location','best','AutoUpdate','off')
                        xtickdates(dates(sammy))
                        plotOrigin('k:')
                        box off
                        title(sprintf('cumulative difference'))
                        
                        subplot(2,1,2)
                        h = plot(dates(sammy), avgcumlossdiff(sammy,thisgroup), 'linewidth', 1);
                        set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
                        legend(h, Yfcstlabel(thisgroup), 'Location','best','AutoUpdate','off')
                        xtickdates(dates(sammy))
                        plotOrigin('k:')
                        box off
                        title('cumulative difference scaled by number of obs')

                        sgtitle(sprintf('{\\bf %s MSE:}\n  %s less %s',  datalabel, modeltype1, modeltype0))
                        wrapthisfigure(thisfig, sprintf('cumlossMSE-%s-%s-vs-%s-group%d', datalabel, ...
                            modeltype0, modeltype1, gg), wrap, [], [], [], [], false);
                    end
                end
            end

            %% CRPS

            loss0 = mat0.fcstYcrps;
            loss1 = mat1.fcstYcrps;
            loss1 = loss1(:,1:theseNhorizons);

            nanny        = (isnan(loss0) | isnan(loss1)) | ~sammy; % note: automatic array expansion w.r.t. sammy
            loss0(nanny) = NaN;
            loss1(nanny) = NaN;

            for hh = 1 : theseNhorizons
                [~, relCRPStstat(hh,d)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
            end

            CRPS0(1:theseNhorizons,d) = mean(loss0,1, 'omitnan');
            CRPS1(1:theseNhorizons,d) = mean(loss1,1, 'omitnan');


            if doPlotsCumLoss && s == 1
                cumlossdiff        = cumsum(loss1 - loss0, 1, 'omitnan');
                nanny              = isnan(loss1 - loss0);
                avgcumlossdiff     = cumlossdiff ./ cumsum(~nanny);
                cumlossdiff(nanny) = NaN;

                grps = {0:4,5:8,9:12,13:15};

                for gg = 1 : length(grps)
                    thisgroup = grps{gg} + 1;
                    thisgroup = thisgroup(thisgroup <= theseNhorizons);
                    if ~isempty(thisgroup)
                       thisfig = figure;
                        subplot(2,1,1)
                        h = plot(dates(sammy), cumlossdiff(sammy,thisgroup), 'linewidth', 1);
                        set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
                        %                         legend(h, Yfcstlabel(thisgroup), 'Location','best','AutoUpdate','off')
                        xtickdates(dates(sammy))
                        plotOrigin('k:')
                        box off
                        title(sprintf('cumulative difference in CPRS'))
                        
                        subplot(2,1,2)
                        h = plot(dates(sammy), avgcumlossdiff(sammy,thisgroup), 'linewidth', 1);
                        set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
                        legend(h, Yfcstlabel(thisgroup), 'Location','best','AutoUpdate','off')
                        xtickdates(dates(sammy))
                        plotOrigin('k:')
                        box off
                        title('cumulative difference scaled by number of obs')

                        sgtitle(sprintf('{\\bf %s CRPS:}\n  %s less %s',  datalabel, modeltype1, modeltype0))
                        wrapthisfigure(thisfig, sprintf('cumlossCRPS-%s-%s-vs-%s-group%d', datalabel, ...
                            modeltype0, modeltype1, gg), wrap, [], [], [], [], false);
                    end
                end
            end

        end

        %% tabulate relative RMSE and CRPS in joint table
        if doMSE
            statname1   = 'MSE';
        else
            statname1   = 'RMSE';
        end
        deltaloss1  = RMSE1 ./ RMSE0;
        deltaTstat1 = relRMSEtstat;

        statname2   = 'CRPS';
        deltaloss2  = CRPS1 ./ CRPS0;
        deltaTstat2 = relCRPStstat;

        tabname    = sprintf('relative-%sand%s-%s-vs-%s-%s.tex', statname1, statname2, ...
            modeltype0, modeltype1, sam(s).label);
        tabcaption = tabname;

        tabrelstats2(tabname, wrap, ...
            statname1, deltaloss1, deltaTstat1, ...
            statname2, deltaloss2, deltaTstat2, ...
            modeltype0, modeltype1, modelpretty0, modelpretty1, ...
            DATALABELS, tabcaption, comparisonNote);

        %% plot relative stats
        if doPlotsRelStats
            crit = norminv(.95,0,1);
            for nn = 1 : length(DATALABELS)

                datalabel = DATALABELS{nn};

                xndx   = 0:MAXHORIZON(nn)-1;
                ndx    = 1:MAXHORIZON(nn);
                these1 = deltaloss1(ndx,nn);
                those1 = abs(deltaTstat1(ndx,nn));
                these2 = deltaloss2(ndx,nn);
                those2 = abs(deltaTstat2(ndx,nn));

                thisfig = figure;
                hold on
                h1 = plot(xndx,these1, 'b-', 'linewidth', 2);
                plot(xndx(those1 > crit),these1(those1 > crit), 'bd', 'linewidth', 2)
                h2 = plot(xndx,these2, 'r-.', 'linewidth', 2);
                plot(xndx(those2 > crit),these2(those2 > crit), 'rd', 'linewidth', 2)
                ylim([.8 1.2])
                xlim(xndx([1 end]))
                plothorzline(1, [], 'k:', 'linewidth', 2)
                legend([h1 h2], statname1, statname2)
                title(sprintf('%s (%s):\n %s / %s', DATALABELS{nn}, sam(s).pretty, ...
                    modelpretty1, modelpretty0))
                wrapthisfigure(thisfig, sprintf('relativeStats-%s-%s-%s-vs-%s', datalabel, sam(s).label, ...
                    modeltype0, modeltype1), wrap, [], [], [], [], false);
            end
        end

    end % samples
end % PAIRS

%% finish / clean up
finishwrap
finishscript
dockAllFigures



function tabrelstats2(tabname, wrap, ...
    statname1, deltaloss1, deltaTstat1, ...
    statname2, deltaloss2, deltaTstat2, ...
    label0, label1, prettylabel0, prettylabel1, DATALABELS, tabcaption, comparisonNote) %#ok<INUSL>


%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
else
    tabdir = fullfile(localtemp, 'foo');
end

[Nhorizons, Nvar] = size(deltaloss1);

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, 2 * Nvar));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nvar, statname1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nvar, statname2);
fprintf(fid, '\\\\ ');
fprintf(fid, '\\cmidrule(lr){%d-%d}', 2, 1+Nvar);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 1+Nvar+1,1+2*Nvar);
fprintf(fid, '\n');
fprintf(fid, 'Horizon $h$');
for n = 1 : Nvar
    fprintf(fid, '& \\multicolumn{1}{c}{%s} ', DATALABELS{n});
end
for n = 1 : Nvar
    fprintf(fid, '& \\multicolumn{1}{c}{%s} ', DATALABELS{n});
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print absolute stat
    for nn = 1 : Nvar
        if isnan(deltaloss1(hh,nn))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss1(hh,nn),  Zstar(deltaTstat1(hh,nn)));
        end
    end
    % print relative stat
    for nn = 1 : Nvar
        if isnan(deltaloss2(hh,nn))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss2(hh,nn),  Zstar(deltaTstat2(hh,nn)));
        end
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


fprintf(fid, 'Note: \n');
fprintf(fid, 'Relative %s and %s of %s model (with %s in denominator).\n', ...
    statname1, statname2, prettylabel1, prettylabel0);

fprintf(fid, '%s \n', comparisonNote);
fprintf(fid, 'Variable mnemonics: RGDP denotes real growth, CPI inflation, and TBILL the interest rate on Treasury bills. (All growth rates are expressed as annualized percentage points of quarterly rates of change.)\n');

fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');

fclose(fid);
type(fullfile(tabdir, tabname))

end

