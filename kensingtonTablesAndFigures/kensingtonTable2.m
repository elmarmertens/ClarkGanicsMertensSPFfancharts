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

datalabel  = 'RGDP';
fontsize   = 18;
Nhorizons  = 16;

%% define set of four models to compare
m = 0;

m = m + 1;
models0(m).type    = 'STATEtrendgapSV';
models0(m).pretty  = 'MDS (SV)';
models0(m).Ndraws  = 3e3;
models0(m).ETlabel = 'none';

models1(m).type    = 'VARprior10STATEtrendgapSV';
models1(m).pretty  = 'SV';
models1(m).Ndraws  = 3e3;
models1(m).ETlabel = 'none';

m = m + 1;
models0(m).type    = 'STATEconst';
models0(m).pretty  = 'CONST';
models0(m).Ndraws  = 3e3;
models0(m).ETlabel = 'none';

models1(m).type    = 'VARprior10STATEconst';
models1(m).pretty  = 'CONST';
models1(m).Ndraws  = 3e3;
models1(m).ETlabel = 'none';


Nmodels = length(models1);
if Nmodels ~= length(models0)
    error houston
end

%% set modeldir
resultdir = localresultsMCMC;

%% define set of eval windows
ss = 1;
sam(ss).start  = datenum(1992,1,1);
sam(ss).stop   = [];
sam(ss).label  = 'since1992';
sam(ss).pretty = '92--22';

ss = 2;
sam(ss).start  = datenum(1992,1,1);
sam(ss).stop   = datenum(2016,12,1);
sam(ss).label  = 'since1992preCOVID';
sam(ss).pretty = '92--16';

Nsam = length(sam);


%% latexwrapper
wrap   = [];
titlename = strcat('kensingtonTable2-', datalabel);

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end


%% allocate memory
[relRMSE, relCRPS, relRMSEtstat, relCRPStstat] = deal(NaN(Nhorizons, Nmodels, Nsam));

%% collect data
for mm = 1 : Nmodels

    % baseline model
    modeltype0   = models0(mm).type;
    modelpretty0 = models0(mm).pretty;
    Ndraws0      = models0(mm).Ndraws;
    ETlabel0     = models0(mm).ETlabel;
    % add ET label to modelpretty
    switch ETlabel0
        case 'none'
            modelpretty0 = models0(mm).pretty;
            isET0 = false;
        case 'binsOnly'
            modelpretty0 = sprintf('ET (bins) with %s', models0(mm).pretty);
            isET0 = true;
        case 'binsAndMeans'
            modelpretty0 = sprintf('ET (bins and means) with %s', models0(mm).pretty);
            isET0 = true;
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel0)
    end
    % add ETlabel to modeltype
    modelETtype0 = strcat('ET', ETlabel0, '-', modeltype0);

    modellabel0  = strcat(datalabel, modeltype0);
    switch ETlabel0
        case 'none'
            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel0, Ndraws0);
        case {'binsOnly', 'binsAndMeans'}
            matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel0, Ndraws0, modellabel0);
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel0)
    end
    mat0        = matfile(fullfile(resultdir, matfilename));


    Yfuture       = mat0.Yfuture;
    Zdata         = mat0.Zdata;
    Nz            = size(Zdata, 2);
    Znanny        = mat0.Znanny;
    Ylabel        = mat0.Ylabel;
    Yfcstlabel    = Ylabel(2:end);
    dates         = mat0.dates;
    T             = mat0.T;
    if isET0
        Tstart = mat0.firstT;
    else
        Tstart = mat0.Tstart;
    end

    if Nhorizons ~= mat0.Nhorizons
        error('Nhorizons mismatch')
    end

    % alt model
    modeltype1   = models1(mm).type;
    modelpretty1 = models1(mm).pretty;
    Ndraws1      = models1(mm).Ndraws;

    ETlabel1     = models1(mm).ETlabel;
    switch ETlabel1
        case 'none'
            modelpretty1 = models1(mm).pretty;
            isET1        = false;
        case 'binsOnly'
            modelpretty1 = sprintf('ET (bins) with %s', models1(mm).pretty);
            isET1        = true;
        case 'binsAndMeans'
            modelpretty1 = sprintf('ET (bins and means) with %s', models1(mm).pretty);
            isET1        = true;
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel0)
    end
    % add ETlabel to modeltype
    modelETtype1 = strcat('ET', ETlabel1, '-', modeltype1);

    % load alt model
    modellabel1  = strcat(datalabel, modeltype1);
    % load model 0
    switch ETlabel1
        case 'none'
            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel1, Ndraws1);
        case {'binsOnly', 'binsAndMeans'}
            matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel1, Ndraws1, modellabel1);
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel1)
    end
    mat1 = matfile(fullfile(resultdir, matfilename));

    % checks
    if ~isequal(dates, mat1.dates)
        error('mismatch between input files for model1 and model2 (dates)')
    end
    if ~isequaln(Zdata, mat1.Zdata)
        warning('mismatch between input files for model1 and model2 (Zdata, %s)', datalabel)
    end
    if ~isequaln(Yfuture, mat1.Yfuture)
        error('mismatch between input files for model1 and model2 (Yfuture)')
    end
    if Nhorizons < mat1.Nhorizons
        error('Nhorizons does not match assumption')
    end
    %% eval windows
    for ss = 1  : Nsam

        if isempty(sam(ss).start)
            samStart     = Tstart;
        else
            samStart = find(dates >= sam(ss).start, 1, 'first');
        end
        sam(ss).start = dates(samStart); % store the chosen value

        if isempty(sam(ss).stop)
            samStop = T;
        else
            samStop = find(dates <= sam(ss).stop, 1, 'last');
        end
        sam(ss).stop = dates(samStop); % store the chosen value

        sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));


        %% RMSE
        loss0 = (mat0.fcstYhaterror).^2;
        loss1 = (mat1.fcstYhaterror).^2;
        % establish common sample (needed for reporting loss levels; dmtest would catch it)
        if ~isequal(isnan(loss0), isnan(loss1))
            error('data mismatch')
        end
        nanny        = (isnan(loss0) | isnan(loss1)) | ~sammy; % note: automatic array expansion w.r.t. sammy
        loss0(nanny) = NaN;
        loss1(nanny) = NaN;

        for hh = 1 : Nhorizons
            if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                [~, relRMSEtstat(hh,mm,ss)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
            else
                relRMSEtstat(hh,mm,ss) = 0;
            end
        end

        rmse0 = sqrt(mean(loss0,1, 'omitnan'));
        rmse1 = sqrt(mean(loss1,1, 'omitnan'));

        relRMSE(:,mm,ss) = rmse1 ./ rmse0;


        %% CRPS

        loss0 = mat0.fcstYcrps;
        loss1 = mat1.fcstYcrps;

        % establish common sample (needed for reporting loss levels; dmtest would catch it)
        if ~isequal(isnan(loss0), isnan(loss1))
            error('data mismatch')
        end
        nanny        = (isnan(loss0) | isnan(loss1)) | ~sammy; % note: automatic array expansion w.r.t. sammy
        loss0(nanny) = NaN;
        loss1(nanny) = NaN;

        for hh = 1 : Nhorizons
            [~, relCRPStstat(hh,mm,ss)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
        end

       crps0 = mean(loss0,1, 'omitnan');
       crps1 = mean(loss1,1, 'omitnan');
             
       relCRPS(:,mm,ss) = crps1 ./ crps0;

    end % samples
end % models

%% ignore DM when 1.00
ndx = round(relRMSE, 2) == 1;
relRMSEtstat(ndx) = 0;

%% table with samples across panels

tabname    = sprintf('table2-%s.tex', datalabel);
tabcaption = tabname;

% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
    %      latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
else
    tabdir = fullfile(localtemp, 'foo');
end

if Nsam ~= 2
    error houston
end


% tabulate
Ncols = 1 + Nsam * Nmodels;
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
% fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('.4', 1, Nsam * Nmodels));
fprintf(fid, '\\toprule\n');
% fprintf(fid, '& & & \\multicolumn{6}{c}{Relative to %s} ', models(1).pretty);
% fprintf(fid, '\\\\\n');
% fprintf(fid, '\\cmidrule(lr){%d-%d}', 4, Ncols);
for nn = 1 : Nmodels
    fprintf(fid, '& \\multicolumn{%d}{c}{%s} ', Nsam, models1(nn).pretty);
end
fprintf(fid, '\\\\\n');
for nn = 1 : Nmodels
    that = 2 + (nn-1) * Nsam;
    fprintf(fid, '\\cmidrule(lr){%d-%d}', that, that + (Nsam - 1));
end
fprintf(fid, '$h$');
for nn = 1 : Nmodels
    for ss = 1 : Nsam
        fprintf(fid, '& \\multicolumn{1}{c}{%s} ', sam(ss).pretty);
    end
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{%d}{c}{PANEL A: RMSE}\n', Ncols);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print ratios
    for mm   = 1 : Nmodels
        for ss   = 1 : Nsam
            this = relRMSE(hh,mm,ss);
            that = relRMSEtstat(hh,mm,ss);
            if isnan(this)
                fprintf(fid, '& \\multicolumn{1}{c}{--} ');
            else
                fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
            end
        end
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{%d}{c}{PANEL B: CRPS}\n', Ncols);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print ratios
    for mm   = 1 : Nmodels
        for ss   = 1 : Nsam
            this = relCRPS(hh,mm,ss);
            that = relCRPStstat(hh,mm,ss);
            if isnan(this)
                fprintf(fid, '& \\multicolumn{1}{c}{--} ');
            else
                fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
            end
        end
    end
    fprintf(fid, '\\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
% fprintf(fid, '\\end{small}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


fprintf(fid, '\\begin{small}\n');
fprintf(fid, 'Note: \n');
fprintf(fid, 'Comparison of forecasts with and without MDS assumption for quarterly GDP growth $h$ steps ahead.\n'); 
fprintf(fid, 'Relative RMSE and CRPS of non-MDS forecasts with corresponding MDS statistics in the denominator.');

fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
fprintf(fid, '``%s''''  and ``%s'''' refer to subsamples beginning in %s and extending until %s, and %s, respectively, (using data for realized values as far as available in %s).', ...
    sam(1).pretty, sam(2).pretty, ...
    datestr(sam(1).start, 'yyyyQQ'), ...
    datestr(sam(1).stop, 'yyyyQQ'), datestr(sam(2).stop, 'yyyyQQ'), datestr(dates(end), 'yyyyQQ'));
fprintf(fid, '\\end{small}\n');

fclose(fid);
type(fullfile(tabdir, tabname))


%% finish / clean up
finishwrap
finishscript
dockAllFigures




