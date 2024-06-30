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
%#ok<*DEFNU>
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc


%% parameters

DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'}; % , 'TBILLcensored'};
fontsize   = 18;


Nvariables = length(DATALABELS);

Ycode      = 'Y';

doEXY1Q4 = true;


%% set modeldir
resultdir = '~/jam/lager/kensingtonStateSpaceDraws/';

%% define models
NGAP               = {'BOP'};
SAMSTART           = {'1968Q4'};
MODELTYPESpretty   = {'MDS','VAR'};

if ~doEXY1Q4
    MODELTYPES         = {'MDStrendHScycleSVt2blockNoiseHS','VAR0trendHScycleSVt2blockNoiseHS'};
    SPFquarterly       = {'-y1q4'};
    prettySPFquarterly = {'(incl Q4 of y1)'};
else
    MODELTYPES         = {'MDStrendHScycleSVt2block','VAR0trendHScycleSVt2block'};
    SPFquarterly       = {'-EXy1q4'};
    prettySPFquarterly = {'(excl Q4 of y1)'};
end

for nn = 1 : length(NGAP)
    for mm = 1 : length(MODELTYPES)
        for ss = 1 : length(SAMSTART)
            for qq = 1 : length(SPFquarterly)
                samStartLabel = SAMSTART{ss};
                models(nn,mm,ss,qq).type   = sprintf('%s%s-Ngap%s-samStart%s', MODELTYPES{mm}, SPFquarterly{qq}, NGAP{nn}, samStartLabel);
                models(nn,mm,ss,qq).Ngap   = NGAP{nn};
                models(nn,mm,ss,qq).pretty = MODELTYPESpretty{mm};
                models(nn,mm,ss,qq).Ndraws = 3e3;
            end
        end
    end
end



%% define set of eval windows
thisSam = 0;

thisSam = thisSam + 1;
sam(thisSam).start  = datenum(1990,1,1);
sam(thisSam).stop   = datenum(2023,10,1);
sam(thisSam).label  = 'fullsampleSince1990';
sam(thisSam).pretty = 'full sample (since 1990)';

thisSam = thisSam + 1;
sam(thisSam).start  = datenum(1990,1,1);
sam(thisSam).stop   = datenum(2019,10,1);
sam(thisSam).label  = 'since1990preCOVID';
sam(thisSam).pretty = '1990 -- 2019';


%% initwrap
titlename = 'prettytableCoverageRates';
if doEXY1Q4
    titelname = strcat(titlename, '-EXY1Q4');
end
initwrap

%% collect some data from first model
% note: assumes everything pulled here is equal across models and datalabels (no checks performed below)
datalabel    = DATALABELS{1};
mm           = 1;
modeltype0   = models(mm).type;
modelpretty0 = models(mm).pretty;
Ndraws0      = models(mm).Ndraws;

modellabel0  = strcat(datalabel, '-', modeltype0);

matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
fullmatfilename = fullfile(resultdir, matfilename);
mat0 = matfile(fullmatfilename);

quantileP     = mat0.quantileP;
Nquantiles    = length(quantileP);
T             = mat0.T;
dates         = mat0.dates;
Tstart        = mat0.Tstart;

%% allocate memory
Ndata        = length(DATALABELS);
maxNhorizons = 17;
lowerNdx     = [3 2];
upperNdx     = [7 8];
bandProb     = quantileP(upperNdx) - quantileP(lowerNdx);
Nbands       = length(bandProb);

[BANDmean, BANDtstat, BANDpval] = deal(NaN(maxNhorizons,Nbands,Ndata,numel(models),length(sam)));

%% loop over samples
for thisSam = 1 : length(sam)

    %% prepare sample choice
    if isempty(sam(thisSam).start)
        samStart = Tstart;
    else
        samStart = find(dates >= sam(thisSam).start, 1, 'first');
    end
    samStart = max(samStart, Tstart);

    if isempty(sam(thisSam).stop)
        samStop = T;
    else
        samStop = find(dates <= sam(thisSam).stop, 1, 'last');
    end
    samStop = min(samStop, T);

    sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));

    sam(thisSam).comparisonNote = sprintf('Evaluation window from %s through %s (and as far as realized values are available).', ...
        datestr(dates(samStart), 'yyyyqq'), datestr(dates(samStop), 'yyyyqq'));

    %% loop over models
    for mm = 1 : numel(models)

        modeltype0   = models(mm).type;
        modelpretty0 = models(mm).pretty;
        Ndraws0      = models(mm).Ndraws;


        %% loop over datalabel

        for dd = 1 : length(DATALABELS)

            datalabel    = DATALABELS{dd};

            %% load data
            modellabel0  = strcat(datalabel, '-', modeltype0);

            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
            fullmatfilename = fullfile(resultdir, matfilename);
            if ~exist(fullmatfilename, 'file')
                warning('%s not found', matfilename);
                continue
            end
            mat0 = matfile(fullmatfilename);

            Ylabel        = mat0.Ylabel;
            Nhorizons     = mat0.Nhorizons;


            fprintf('Processing %s ... \n', datalabel)



            %% collect coverage rates for uncertainty bands
            yQuantiles = mat0.(sprintf('fcst%squantiles', Ycode));
            yFuture    = mat0.(sprintf('%sfuture', Ycode));

            % cut sample for yFuture and yQuantiles
            yFuture    = yFuture(sammy,:);
            yQuantiles = yQuantiles(sammy,:,:);

            % apply NaN for Yfuture dates outside the intended horizon
            for h = 2 : Nhorizons
                thisH = h - 1; % first horizon is nowcast
                yFuture(end-thisH+1:end,h) = NaN;
            end

            [bandMean, bandTstat, bandPval] = deal(NaN(Nhorizons,Nbands));
            for ii = 1 : Nbands
                [bandMean(:,ii), bandTstat(:,ii), bandPval(:,ii)] = coverageBands(yFuture, yQuantiles(:,:,lowerNdx(ii)), yQuantiles(:,:,upperNdx(ii)), bandProb(ii));
            end

            BANDmean(1:Nhorizons,:,dd,mm,thisSam)  = bandMean;
            BANDtstat(1:Nhorizons,:,dd,mm,thisSam) = bandTstat;
            BANDpval(1:Nhorizons,:,dd,mm,thisSam)  = bandPval;

        end % datalabel d
    end % models m
end % samples

%% tabulate values and output in TeX -- separate tables per model
Nstats  = 2;
for thisSam = 1 : length(sam)
    for mm = 1 : numel(models)

        tableNote  = sam(thisSam).comparisonNote; % sprintf('Actual vs predicted coverage rates (in percent) generated by %s model.\n%s', MODELTYPESpretty{mm}, sam(thisSam).comparisonNote);

        tabname    = sprintf('table-CoverageRateBands-%s-%s.tex', models(mm).type, sam(thisSam).label);
        tabcaption = sprintf('Coverage rates of uncertainty bands (%s, %s)', models(mm).pretty, sam(thisSam).pretty);
        if ~isempty(wrap)
            tabdir = wrap.dir;
            latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
        else
            tabdir = fullfile(localtemp, 'foo');
        end
        fid = fopen(fullfile(tabdir, tabname), 'wt');
        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, Nstats * Ndata));
        fprintf(fid, '\\toprule\n');
        % header: dates row
        for dd = 1 : Ndata
            fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, DATALABELS{dd});
        end
        fprintf(fid, '\\\\\n');
        offset = 1;
        for dd = 1 : Ndata
            fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats);
            offset = offset + Nstats;
        end
        fprintf(fid, '\n');
        % header stats row
        fprintf(fid, '$h$ ');
        for dd = 1 : Ndata
            fprintf(fid, ' & \\multicolumn{1}{c}{68\\%%} & \\multicolumn{1}{c}{90\\%%}');
        end
        fprintf(fid, '\\\\\n');
        fprintf(fid, '\\midrule\n');
        for hh = 1 : maxNhorizons
            fprintf(fid, '%d ', hh-1);
            for dd = 1 : Ndata
                for ss = 1 : Nstats
                    if ~isnan(BANDmean(hh,ss,dd,mm,thisSam))
                        fprintf(fid, '& %6.2f%s ', BANDmean(hh,ss,dd,mm,thisSam), Zstar(BANDtstat(hh,ss,dd,mm,thisSam)));
                    else
                        fprintf(fid, '& \\multicolumn{1}{c}{--} ');
                    end
                end
            end
            fprintf(fid, '\\\\\n');
        end
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        fprintf(fid, 'Note: \n');
        fprintf(fid, 'Coverage rates for uncertainty bands with nominal levels of 68\\%% and 90\\%% for out-of-sample forecasts at quarterly forecast horizons, $h$.\n');
        fprintf(fid, '%s \n', tableNote);
        fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
        fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
        fclose(fid);
        type(fullfile(tabdir, tabname))

    end % models
end % samples

%% tabulate values and output in TeX -- one table with separate panels per model
Nstats  = 2;
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; % for panel labels

for thisSam = 1 : length(sam)

    tableNote  = sprintf('%s', sam(thisSam).comparisonNote);

    if ~doEXY1Q4
        tabname    = sprintf('table-CoverageRateBands-%s.tex', sam(thisSam).label);
    else
        tabname    = sprintf('table-CoverageRateBands-EXY1Q4-%s.tex', sam(thisSam).label);
    end
    tabcaption = sprintf('Coverage rates of uncertainty bands (%s)', sam(thisSam).pretty);
    if ~isempty(wrap)
        tabdir = wrap.dir;
        latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
    else
        tabdir = fullfile(localtemp, 'foo');
    end
    fid = fopen(fullfile(tabdir, tabname), 'wt');
    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, Nstats * Ndata));
    Ncols = 1 + Ndata * Nstats;
    fprintf(fid, '\\toprule\n');
    % header: datalabel row
    for dd = 1 : Ndata
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, DATALABELS{dd});
    end
    fprintf(fid, '\\\\\n');
    offset = 1;
    for dd = 1 : Ndata
        fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats);
        offset = offset + Nstats;
    end
    fprintf(fid, '\n');
    % header stats row
    fprintf(fid, '$h$ ');
    for dd = 1 : Ndata
        fprintf(fid, ' & \\multicolumn{1}{c}{68\\%%} & \\multicolumn{1}{c}{90\\%%}');
    end
    fprintf(fid, '\\\\\n');
    for mm = 1 : numel(models)
        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\multicolumn{%d}{c}{PANEL %s: %s Model} ', Ncols, alphabet(mm), MODELTYPESpretty{mm});
        fprintf(fid, '\\\\\n');
        fprintf(fid, '\\midrule\n');
        for hh = 1 : maxNhorizons
            fprintf(fid, '%d ', hh-1);
            for dd = 1 : Ndata
                for ss = 1 : Nstats
                    if ~isnan(BANDmean(hh,ss,dd,mm,thisSam))
                        fprintf(fid, '& %6.2f%s ', BANDmean(hh,ss,dd,mm,thisSam), Zstar(BANDtstat(hh,ss,dd,mm,thisSam)));
                    else
                        fprintf(fid, '& \\multicolumn{1}{c}{--} ');
                    end
                end
            end
            fprintf(fid, '\\\\\n');
        end
    end % models
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, 'Note: \n');
    fprintf(fid, 'Coverage rates for uncertainty bands with nominal levels of 68\\%% and 90\\%% for out-of-sample forecasts at quarterly forecast horizons, $h$.\n');
    fprintf(fid, '%s \n', tableNote);
    fprintf(fid, 'Reflecting the availability of annual SPF forecasts, forecasts for inflation in CPI and GDP prices are evaluated only up to $h=12$, and $h=8$, respectively.\n');
    fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
    fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
    fclose(fid);
    type(fullfile(tabdir, tabname))

end % samples

%% finish
finishwrap
finishscript

%% functions
function [coverMean, coverTstat, coverPval] = coverageBands(yFuture, lowerBand, upperBand, bandProb)
% function [tstat, pval] = coverageBands(yFuture, lowerBand, upperBand, bandProb)
%
% yFuture     : T x Nhorizons
% yQuantiles  : T x Nhorizons x Nbands

% tstat       : Nhorizons x Nbands
% pval        : Nhorizons x Nbands

Nhorizons  = size(yFuture,2);
Nbands     = length(bandProb);

quantNdx   = double((yFuture > lowerBand) & (yFuture < upperBand)) * 100;
quantNdx(isnan(yFuture) | isnan(lowerBand) | isnan(upperBand)) = NaN;

% mean coverage rate
coverMean   = permute(mean(quantNdx,1, 'omitnan'), [2 3 1]); % permute is more precise than squeeze

% test for significance with Newey-West regressions
[coverTstat, coverPval] = deal(NaN(size(coverMean)));
for q = 1 : Nbands

    for n = 1 : Nhorizons

        coverrate = quantNdx(:,n,q) - bandProb(q);
        nanny     = ~isnan(coverrate);
        reggae    = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + n));
        coverTstat(n,q) = reggae.tstat;
        coverPval(n,q)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);

    end
end

end %  function coverageBands

function tabulateCoverageRate(tabname, wrap, ...
    quantileP, coverStat, coverTstat, ...
    bandProb, bandStat, bandTstat, ...
    tabcaption, tableNote, quantHeader, bandHeader)


%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
else
    tabdir = fullfile(localtemp, 'foo');
end

[Nhorizons, Nquantiles] = size(coverStat);
if ~isequal(Nquantiles, length(quantileP))
    error('quantileP and coverStat do not match')
end

Nbands = length(bandProb);

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('.5', 1, Nquantiles+Nbands));
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Horizon & \\multicolumn{%d}{c}{%s} ', Nquantiles, quantHeader);
fprintf(fid, '& \\multicolumn{%d}{c}{%s} ', Nbands, bandHeader);
fprintf(fid, '\\\\ ');
% fprintf(fid, '\\cmidrule(lr){%d-%d}', 1, 2);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 2, 1+Nquantiles);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 1+Nquantiles+1, 1+Nquantiles+Nbands);
fprintf(fid, '\n');
fprintf(fid, '$h$');
fprintf(fid, '& %5.2f ', quantileP);
fprintf(fid, '& %5.2f ', bandProb);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print coverage rates
    for nn = 1 : Nquantiles
        if isnan(coverStat(hh,nn))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', coverStat(hh,nn),  Zstar(coverTstat(hh,nn)));
        end
    end
    % print coverage rates for bands
    for nn = 1 : Nbands
        if isnan(bandStat(hh,nn))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', bandStat(hh,nn),  Zstar(bandTstat(hh,nn)));
        end
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


fprintf(fid, 'Note: \n');
fprintf(fid, 'Coverage rates for uncertainty bands with nominal levels of 68\\%% and 90\\%% for out-of-sample forecasts at quarterly forecast horizons, $h$.\n');
fprintf(fid, '%s \n', tableNote);
fprintf(fid, 'Reflecting the availability of annual SPF forecasts, forecasts for inflation in CPI and GDP prices are evaluated only up to $h=12$, and $h=8$, respectively.\n');
fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');

fclose(fid);
type(fullfile(tabdir, tabname))

end % function tabulateCoverageRate