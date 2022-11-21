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

DATALABELS  = {'PGDP', 'UNRATE'};
Nvar        = length(DATALABELS);
Nhorizons   = 16;

fontsize    = 18;

%% define baseline and alt models

modeltype0   = 'STATEtrendgapSV';
modelpretty0 = 'SV';

modeltype1   = 'STATEtrendgapSV';
modelpretty1 = 'SV w/ET';
ETlabel1     = 'binsOnly';

% modeltype2   = 'VARprior5STATEtrendgapSV';
% modelpretty2 = 'non-MDS';
% ETlabel2     = 'none';

resultdir    = localresultsMCMC;
Ndraws       = 3e3;

Nmodels      = 1; % that is: models relative to SV

%% define set of eval windows
s = 1;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = [];
sam(s).label  = 'since1992';

s = 2;
sam(s).start  = datenum(1992,1,1);
sam(s).stop   = datenum(2016,12,1);
sam(s).label  = 'since1992preCOVID';

Nsam = length(sam);


%% latexwrapper
wrap   = [];
titlename = 'kensingtonTable3';

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% allocate memory
[relRMSE, relCRPS, relRMSEtstat, relCRPStstat] = deal(NaN(Nhorizons, Nsam, Nmodels, Nvar));

samlabel = cell(Nsam, Nvar);

for d = 1 : length(DATALABELS)
    datalabel = DATALABELS{d};

    %% collect data (baseline model)
    modellabel0   = strcat(datalabel, modeltype0);
    matfilename   = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel0, Ndraws);
    mat0          = matfile(fullfile(resultdir, matfilename));

    Yfuture       = mat0.Yfuture;
    Zdata         = mat0.Zdata;
    Nz            = size(Zdata, 2);
    Znanny        = mat0.Znanny;
    dates         = mat0.dates;
    T             = mat0.T;
    Tstart        = mat0.Tstart;

    theseNhorizons  = mat0.Nhorizons;

    modellabel1   = strcat(datalabel, modeltype1);
    switch ETlabel1
        case 'none'
            matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel1, Ndraws);
        case {'binsOnly', 'binsAndMeans'}
            matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel1, Ndraws, modellabel1);
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel1)
    end
    mat1          = matfile(fullfile(resultdir, matfilename));

    %     modellabel2   = strcat(datalabel, modeltype2);
    %     switch ETlabel2
    %         case 'none'
    %             matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel2, Ndraws);
    %         case {'binsOnly', 'binsAndMeans'}
    %             matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel2, Ndraws, modellabel2);
    %         otherwise
    %             error('ETlabel <<%s>> not recognized', ETlabel1)
    %     end
    %     mat2          = matfile(fullfile(resultdir, matfilename));

    % todo: check data across mat files

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

        %% RMSE

        loss0         = (mat0.fcstYhaterror).^2;
        loss1         = (mat1.fcstYhaterror).^2;
        %         loss2         = (mat2.fcstYhaterror).^2;

        nanny         = isnan(loss0)  | isnan(loss1)  | (~sammy);
        
        loss0(nanny)  = NaN;
        loss1(nanny)  = NaN;
        %         loss2(nanny)  = NaN;

        this0         = sqrt(mean(loss0,1, 'omitnan'));
        this1         = sqrt(mean(loss1,1, 'omitnan'));
        %         this2         = sqrt(mean(loss2,1, 'omitnan'));

        relRMSE(1:theseNhorizons,s,1,d) = this1 ./ this0;
        %         relRMSE(1:theseNhorizons,s,2,d) = this2 ./ this0;

        for hh = 1 : theseNhorizons
            if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                [~, relRMSEtstat(hh,s,1,d)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
            else
                relRMSEtstat(hh,s,1,d) = 0;
            end
        end

        %% samlabel
        nanny          = all(nanny, 2);
        firstyear      = year(min(dates(~nanny)));
        lastyear       = year(sam(s).stop); % max(dates) will return last obs for which realized value is available
        if firstyear > 2000
            firstyear = firstyear - 2000;
        else
            firstyear = firstyear - 1900;
        end
        if lastyear > 2000
            lastyear = lastyear - 2000;
        else
            lastyear = lastyear - 1900;
        end
        samlabel{s, d} = sprintf('%02d--%02d', firstyear, lastyear);

        %% CRSP

        loss0         = mat0.fcstYcrps;
        loss1         = mat1.fcstYcrps;
        %         loss2         = mat2.fcstYcrps;

        nanny         = isnan(loss0)  | isnan(loss1) | (~sammy);

        loss0(nanny)  = NaN;
        loss1(nanny)  = NaN;
        %         loss2(nanny)  = NaN;

        this0         = mean(loss0,1, 'omitnan');
        this1         = mean(loss1,1, 'omitnan');
        %         this2         = mean(loss2,1, 'omitnan');

        relCRPS(1:theseNhorizons,s,1,d) = this1 ./ this0;
        %         relCRPS(1:theseNhorizons,s,2,d) = this2 ./ this0;

        for hh = 1 : theseNhorizons
            if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                [~, relCRPStstat(hh,s,1,d)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
            else
                relCRPStstat(hh,s,1,d) = 0;
            end
        end
    end % sam
end % datalabel

%% tabulate
tabname    = sprintf('table3compact-%s-%s.tex', DATALABELS{1}, DATALABELS{2});
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
Ncols = 1 + (Nsam * Nmodels) * Nvar * 2;
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
% fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('.4', 1, Ncols - 1));
fprintf(fid, '\\toprule\n');

% data header
for d = 1 : Nvar
    fprintf(fid, '& \\multicolumn{%d}{c}{%s} ', Nsam * Nmodels * 2, DATALABELS{d});
end
fprintf(fid, '\\\\\n');
that = 1;
fprintf(fid, '\\cmidrule(lr){%d-%d}', that + 1, that + Nsam * Nmodels * 2);
that = 1 + Nsam * Nmodels * 2;
fprintf(fid, '\\cmidrule(lr){%d-%d}', that + 1, that + Nsam * Nmodels * 2);

% model header
fprintf(fid, '& \\multicolumn{2}{c}{%s} ', 'RMSE');
fprintf(fid, '& \\multicolumn{2}{c}{%s} ', 'CRPS');
fprintf(fid, '& \\multicolumn{2}{c}{%s} ', 'RMSE');
fprintf(fid, '& \\multicolumn{2}{c}{%s} ', 'CRPS');
fprintf(fid, '\\\\\n');
for nn = 1 : Nmodels * Nsam * 2
    that = 2 + (nn-1) * 2;
    fprintf(fid, '\\cmidrule(lr){%d-%d}', that, that + 1);
end
fprintf(fid, '$h$');
for d = 1 : Nvar
    for nn = 1 : Nmodels * 2
     fprintf(fid, '& \\multicolumn{1}{c}{%s} ', samlabel{1,d});
     fprintf(fid, '& \\multicolumn{1}{c}{%s} ', samlabel{2,d});
    end
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

% fprintf(fid, '\\multicolumn{%d}{c}{PANEL A: RMSE}\n', Ncols);
% fprintf(fid, '\\\\\n');
% fprintf(fid, '\\midrule\n');
for hh = 1 : Nhorizons
    fprintf(fid, '%d ', hh-1);
    for dd   = 1 : Nvar
        mm = 1;
        % RMSE
        for ss   = 1 : Nsam
            this = relRMSE(hh,ss,mm,dd);
            that = relRMSEtstat(hh,ss,mm,dd);
            if isnan(this)
                fprintf(fid, '& \\multicolumn{1}{c}{--} ');
            else
                fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
            end
        end
        % CRPS
        for ss   = 1 : Nsam
            this = relCRPS(hh,ss,mm,dd);
            that = relCRPStstat(hh,ss,mm,dd);
            if isnan(this)
                fprintf(fid, '& \\multicolumn{1}{c}{--} ');
            else
                fprintf(fid, '& %6.2f%s ', this,  Zstar(that));
            end
        end
    end
    fprintf(fid, '\\\\\n');
end
% fprintf(fid, '\\midrule\n');

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
% fprintf(fid, '\\end{small}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


fprintf(fid, '\\begin{small}\n');
fprintf(fid, 'Note: \n');
fprintf(fid, 'RMSE and CRPS of %s relative to %s (in denominator).\n', modelpretty1, modelpretty0);

fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');
fprintf(fid, 'Evaluation windows extend from %s and %s until %s and %s, respectively.\n', ...
     datestr(datenum(1992,1,1), 'yyyyQQ'), ...
     datestr(datenum(2009,4,1), 'yyyyQQ'), ...
     datestr(datenum(2022,4,1), 'yyyyQQ'), ...
     datestr(datenum(2016,10,1), 'yyyyQQ') ...
     );
% fprintf(fid, 'and as indicated by the column label ``%s,'''' ``%s,'''' ``%s,'''' and ``%s.''''\n', ...
%       samlabel{1,1}, samlabel{1,2}, samlabel{1,2}, samlabel{2,2});
fprintf(fid, '(using data for realized values as far as available by %s).\n', datestr(dates(end), 'yyyyQQ'));
fprintf(fid, '\\end{small}\n');

fclose(fid);
type(fullfile(tabdir, tabname))


%% finish / clean up
finishwrap
finishscript
dockAllFigures




