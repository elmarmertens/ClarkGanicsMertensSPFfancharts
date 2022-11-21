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

%% define set of four models to compare
m = 0;

m = m + 1;
models(m).type    = 'STATEtrendgapSV';
models(m).pretty  = 'SV';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'none';

m = m + 1;
models(m).type    = 'STATEtrendgapSV';
models(m).pretty  = 'SV w/ET';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsOnly';

m = m + 1;
models(m).type    = 'STATEconst';
models(m).pretty  = 'CONST';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'none';

m = m + 1;
models(m).type    = 'STATEconst';
models(m).pretty  = 'CONST w/ET';
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsOnly';

Nmodels = length(models);

%% set modeldir
resultdir = localresultsMCMC;

%% define set of eval windows
s = 1;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = [];
sam(s).label  = 'since1992';
sam(s).pretty = '92--22';

s = 2;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = datenum(2016,12,1);
sam(s).label  = 'since1992preCOVID';
sam(s).pretty = '92--16';

Nsam = length(sam);


%% latexwrapper
wrap   = [];
titlename = strcat('kensingtonTable1-', datalabel);

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end



%% collect data (baseline model)
m0           = 1;

modeltype0   = models(m0).type;
modelpretty0 = models(m0).pretty;
Ndraws0      = models(m0).Ndraws;
ETlabel0     = models(m0).ETlabel;
% add ET label to modelpretty
switch ETlabel0
    case 'none'
        modelpretty0 = models(m0).pretty;
        isET0 = false;
    case 'binsOnly'
        modelpretty0 = sprintf('ET (bins) with %s', models(m0).pretty);
        isET0 = true;
    case 'binsAndMeans'
        modelpretty0 = sprintf('ET (bins and means) with %s', models(m0).pretty);
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

Nhorizons = mat0.Nhorizons;



%% allocate memory
[MAE0, MAE1, RMSE0, RMSE1, CRPS0, CRPS1]  = deal(NaN(Nhorizons, Nmodels, Nsam));
[relRMSEtstat, relMAEtstat, relCRPStstat] = deal(NaN(Nhorizons, Nmodels, Nsam));

%% eval windows
for s = 1  : Nsam

    if isempty(sam(s).start)
        samStart     = Tstart;
    else
        samStart = find(dates >= sam(s).start, 1, 'first');
    end
    sam(s).start = dates(samStart); % store the chosen value

    if isempty(sam(s).stop)
        samStop = T;
    else
        samStop = find(dates <= sam(s).stop, 1, 'last');
    end
    sam(s).stop = dates(samStop); % store the chosen value

    sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));

    %% RMSE and CRPS of baseline
    % RMSE
    loss0         = (mat0.fcstYhaterror).^2; % use MCMC mean Yhat, not YhatRB since ET has not RB
    nanny         = isnan(loss0)  | (~sammy);
    loss0(nanny)  = NaN;
    RMSE0(:,m0,s) = sqrt(mean(loss0,1, 'omitnan'));

    % CRPS
    loss0         = mat0.fcstYcrps;
    nanny         = isnan(loss0)  | (~sammy); % TODO: make sure that isnan identical across models
    loss0(nanny)  = NaN;
    CRPS0(:,m0,s) = mean(loss0,1, 'omitnan');

    %% collect data from altmodels
    for m1 = 2 : Nmodels

        modeltype1   = models(m1).type;
        modelpretty1 = models(m1).pretty;
        Ndraws1      = models(m1).Ndraws;

        ETlabel1     = models(m1).ETlabel;
        switch ETlabel1
            case 'none'
                modelpretty1 = models(m1).pretty;
                isET1        = false;
            case 'binsOnly'
                modelpretty1 = sprintf('ET (bins) with %s', models(m1).pretty);
                isET1        = true;
            case 'binsAndMeans'
                modelpretty1 = sprintf('ET (bins and means) with %s', models(m1).pretty);
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
                [~, relRMSEtstat(hh,m1,s)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
            else
                relRMSEtstat(hh,m1,s) = 0;
            end
        end

        RMSE0(:,m1,s) = sqrt(mean(loss0,1, 'omitnan'));
        RMSE1(:,m1,s) = sqrt(mean(loss1,1, 'omitnan'));

        checkdiff(RMSE0(:,m1,s), RMSE0(:,m0,s));


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
            [~, relCRPStstat(hh,m1,s)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
        end

        CRPS0(:,m1,s) = mean(loss0,1, 'omitnan');
        CRPS1(:,m1,s) = mean(loss1,1, 'omitnan');

        checkdiff(CRPS0(:,m1,s), CRPS0(:,m0,s));

    end % samples
end % models


%% table with samples across panels

statname1      = 'RMSE';
relstat1       =  RMSE1 ./ RMSE0;
relstat1       = relstat1(:,2:end,:);
relstat1tstat  = relRMSEtstat(:,2:end,:);

statname2      = 'CRPS';
relstat2       = CRPS1 ./ CRPS0;
relstat2       = relstat2(:,2:end,:);
relstat2tstat  = relCRPStstat(:,2:end,:);

%% ignore DM when 1.00
ndx = round(relstat1, 2) == 1;
relstat1tstat(ndx) = 0;
ndx = round(relstat2, 2) == 1;
relstat2tstat(ndx) = 0;

%% tables
Nmodels1    = Nmodels - 1;

tabname    = sprintf('table1-%s.tex', datalabel);
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
Ncols = 1 + 2 * Nmodels;
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
% fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('.4', 1, 2 * Nmodels));
fprintf(fid, '\\toprule\n');
fprintf(fid, '& & & \\multicolumn{6}{c}{Relative to %s (in denominator)} ', models(1).pretty);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\cmidrule(lr){%d-%d}', 4, Ncols);
for nn = 1 : Nmodels
    fprintf(fid, '& \\multicolumn{2}{c}{%s} ', models(nn).pretty);
end
fprintf(fid, '\\\\\n');
for nn = 1 : Nmodels
    that = 2 + (nn-1) * 2;
    fprintf(fid, '\\cmidrule(lr){%d-%d}', that, that + 1);
end
fprintf(fid, '$h$');
for nn = 1 : Nmodels
    fprintf(fid, '& \\multicolumn{1}{c}{%s} ', sam(1).pretty);
    fprintf(fid, '& \\multicolumn{1}{c}{%s} ', sam(2).pretty);
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{%d}{c}{PANEL A: %s}\n', Ncols, statname1);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print baseline levels
    for ss = 1 : Nsam
        this = RMSE0(hh,1,ss);
        if isnan(this)
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', this);
        end
    end
    % print ratios
    for mm   = 1 : Nmodels1
        for ss   = 1 : Nsam
            this = relstat1(hh,mm,ss);
            that = relstat1tstat(hh,mm,ss);
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

fprintf(fid, '\\multicolumn{%d}{c}{PANEL B: %s}\n', Ncols, statname2);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    % print baseline levels
    for ss = 1 : Nsam
        this = CRPS0(hh,1,ss);
        if isnan(this)
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', this);
        end
    end
    % print ratios
    for mm   = 1 : Nmodels1
        for ss   = 1 : Nsam
            this = relstat2(hh,mm,ss);
            that = relstat2tstat(hh,mm,ss);
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
fprintf(fid, '\\vspace{-.5\\baselineskip}\n');


fprintf(fid, '\\begin{small}\n');
fprintf(fid, 'Note: \n');
fprintf(fid, 'Forecasts for quarterly GDP growth $h$ steps ahead over ');
fprintf(fid, 'subsamples extending from %s until %s and %s, respectively (using data for realized values as far as available in %s).\n', ...
    datestr(sam(1).start, 'yyyyQQ'), datestr(sam(1).stop, 'yyyyQQ'), datestr(sam(2).stop, 'yyyyQQ'), datestr(dates(end), 'yyyyQQ'));
% fprintf(fid, 'Ratios of %s and %s use %s model in denominator.\n', ...
%     statname1, statname2, modelpretty0);
fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
% fprintf(fid, '``%s'''' and ``%s'''' refer to subsamples extending from %s until %s and %s, respectively (using data for realized values as far as available in %s).', ...
%     sam(1).pretty, sam(2).pretty, datestr(sam(1).start, 'yyyyQQ'), datestr(sam(1).stop, 'yyyyQQ'), datestr(sam(2).stop, 'yyyyQQ'), datestr(dates(end), 'yyyyQQ'));
fprintf(fid, '\\end{small}\n');

fclose(fid);
type(fullfile(tabdir, tabname))


%% finish / clean up
finishwrap
finishscript
dockAllFigures




