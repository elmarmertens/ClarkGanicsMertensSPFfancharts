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

%% clear workspace
clear variables

close all
fclose all;
clc

%% parameters

DATALABELS = {'RGDP', 'PGDP'}; %{'UNRATE'};  % 
DATALABELS_CONC = [DATALABELS{:}];
DATALABELS_pretty_all = {'GDP growth','Inflation in the GDP price index','Unemployment rate'};
DATALABELS_pretty = DATALABELS_pretty_all(1:2); % adjust this depending on which table you want to produce - make sure it matches DATALABELS

fontsize   = 18;

Nhorizons  = 16;
Nvariables = length(DATALABELS);

doMSE = false; % to report relative MSE vs RMSE

doYAVG = false;

if doYAVG
    datatype = 'YAVG';
else
    datatype = [];
end
    
%% define set of models to compare

PAIRS = {[1  8  10 17],  ...     % SV:    raw      vs GBMVS      and CONST: raw      vs GBMVS
         [2  8  11 17],  ...     % SV:    binsOnly vs GBMVS      and CONST: binsOnly vs GBMVS
         [1  8  2  8],   ...     % SV:    raw      vs GBMVS      and SV:    binsOnly vs GBMVS
         [10 17 11 17],  ...     % CONST: raw      vs GBMVS      and CONST: binsOnly vs GBMVS
         [1  9  10 18]};         % SV:    raw      vs GBVSSPFAnn and CONST: raw      vs GBVSSPFAnn

% PAIRS = {[2 6 11 15], ...        % SV: binsOnly vs GBM   and CONST: binsOnly vs GBM
%          [2 7 11 16], ...        % SV: binsOnly vs GBMV  and CONST: binsOnly vs GBMV
%          [2 8 11 17], ...        % SV: binsOnly vs GBMVS and CONST: binsOnly vs GBMVS
%          [2 4 11 13], ...        % SV: binsOnly vs NM    and CONST: binsOnly vs NM
%          [2 5 11 14]};           % SV: binsOnly vs NMV   and CONST: binsOnly vs NMV

m = 0;

thisType   = 'STATEtrendgapSV';
thisPretty = 'SV';

% 1
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'none';
models(m).header  = 'SV';

% 2
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsOnly';
models(m).header  = 'ET (bins)';

% 3
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsAndMeans';

% 4
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'NMeansOnly';

% 5
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'NMeansAndVars';

% 6
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansOnly';

% 7
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansAndVars';

% 8
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansVarsAndSkew';
models(m).header  = 'ET (moments)';

% 9
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBVarsSkewAndSPFAnnPoints';
models(m).header  = 'ET (moments)';


thisType   = 'STATEconst';
thisPretty = 'CONST';

% 10
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'none';
models(m).header  = 'CONST';

% 11
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsOnly';
models(m).header  = 'ET (bins)';

% 12
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'binsAndMeans';

% 13
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'NMeansOnly';

%14
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'NMeansAndVars';

% 15
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansOnly';

% 16
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansAndVars';

% 17
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBMeansVarsAndSkew';
models(m).header  = 'ET (moments)';

% 18
m = m + 1;
models(m).type    = thisType;
models(m).pretty  = thisPretty;
models(m).Ndraws  = 3e3;
models(m).ETlabel = 'GBVarsSkewAndSPFAnnPoints';
models(m).header  = 'ET (moments)';

%% set modeldir
resultdir = localresultsMCMC;
ETresultdir = localresultsET;

%% define set of eval windows

s = 1;

if strcmpi(DATALABELS,'UNRATE')
    sam(s).start  = datenum(2009,4,1);
    sam(s).stop   = [];
    sam(s).label  = 'since2009';
    sam(s).pretty = '09--22';
else
    sam(s).start  = datenum(1992,1,1);
    sam(s).stop   = [];
    sam(s).label  = 'since1992';
    sam(s).pretty = '92--22';
end

samstartpretty   = string(datestr(sam(s).start,'yyyyQQ'));
samendlongpretty = string(datestr(datenum(2022,4,1),'yyyyQQ'));

s = 2;
if strcmpi(DATALABELS,'UNRATE')
    sam(s).start  = datenum(2009,4,1);
    sam(s).label  = 'preCOVIDsince2009';
    sam(s).pretty = '09--16';
else
    sam(s).start  = datenum(1992,1,1);
    sam(s).label  = 'preCOVIDsince1992';
    sam(s).pretty = '92--16';
end
sam(s).stop   = datenum(2016,12,1);

samendshortpretty = string(datestr(sam(s).stop,'yyyyQQ'));

Nsam = length(sam);

%% latexwrapper
wrap   = [];
titlename = strcat('tabulateRelativeEvalStatsETSPFCDFforPaper',DATALABELS_CONC);
if doYAVG
    titlename = strcat(titlename, '-YAVG');
end

initwrap
if isempty(wrap) && ~isdesktop
    initwrap
end

%% loop over model pairs

for pp = 5%1 : length(PAIRS)

    close all
    thispair = PAIRS{pp};
    m0       = thispair(1);
    m1       = thispair(2);
    m02      = thispair(3);
    m12      = thispair(4);

    modeltype0   = models(m0).type;
    Ndraws0      = models(m0).Ndraws;
    ETlabel0     = models(m0).ETlabel;

    modeltype1   = models(m1).type;
    Ndraws1      = models(m1).Ndraws;
    ETlabel1     = models(m1).ETlabel;

    modeltype02   = models(m02).type;
    Ndraws02      = models(m02).Ndraws;
    ETlabel02     = models(m02).ETlabel;

    modeltype12   = models(m12).type;
    Ndraws12      = models(m12).Ndraws;
    ETlabel12     = models(m12).ETlabel;

    % add ET label to modelpretty
    switch ETlabel0
        case 'none'
            modelpretty0 = [];
        case 'binsOnly'
            modelpretty0 = [];
        case 'binsAndMeans'
            modelpretty0 = [];
        case 'NMeansOnly'
            modelpretty0 = sprintf('mean of a normal distribution fitted to the SPF probability bins');
        case 'NMeansAndVars'
            modelpretty0 = sprintf('mean and variance of a normal distribution fitted to the SPF probability bins');
        case 'GBMeansOnly'
            modelpretty0 = sprintf('mean of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansAndVars'
            modelpretty0 = sprintf('mean and variance of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansVarsAndSkew'
            modelpretty0 = sprintf('mean, variance, and skewness of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBVarsSkewAndSPFAnnPoints'
            modelpretty0 = sprintf('variance and skewness of a generalized beta distribution fitted to the SPF probability bins, keeping the annual SPF point forecasts');
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel0)
    end

    switch ETlabel02
        case 'none'
            modelpretty02 = [];
        case 'binsOnly'
            modelpretty02 = [];
        case 'binsAndMeans'
            modelpretty02 = [];
        case 'NMeansOnly'
            modelpretty02 = sprintf('mean of a normal distribution fitted to the SPF probability bins');
        case 'NMeansAndVars'
            modelpretty02 = sprintf('mean and variance of a normal distribution fitted to the SPF probability bins');
        case 'GBMeansOnly'
            modelpretty02 = sprintf('mean of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansAndVars'
            modelpretty02 = sprintf('mean and variance of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansVarsAndSkew'
            modelpretty02 = sprintf('mean, variance, and skewness of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBVarsSkewAndSPFAnnPoints'
            modelpretty02 = sprintf('variance and skewness of a generalized beta distribution fitted to the SPF probability bins, keeping the annual SPF point forecasts');
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel02)
    end


    switch ETlabel1
        case 'none'
            modelpretty1 = [];
        case 'binsOnly'
            modelpretty1 = [];
        case 'binsAndMeans'
            modelpretty1 = [];
        case 'NMeansOnly'
            modelpretty1 = sprintf('mean of a normal distribution fitted to the SPF probability bins');
        case 'NMeansAndVars'
            modelpretty1 = sprintf('mean and variance of a normal distribution fitted to the SPF probability bins');
        case 'GBMeansOnly'
            modelpretty1 = sprintf('mean of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansAndVars'
            modelpretty1 = sprintf('mean and variance of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansVarsAndSkew'
            modelpretty1 = sprintf('mean, variance, and skewness of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBVarsSkewAndSPFAnnPoints'
            modelpretty1 = sprintf('variance and skewness of a generalized beta distribution fitted to the SPF probability bins, keeping the annual SPF point forecasts');
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel1)
    end

    switch ETlabel12
        case 'none'
            modelpretty12 = [];
        case 'binsOnly'
            modelpretty12 = [];
        case 'binsAndMeans'
            modelpretty12 = [];
        case 'NMeansOnly'
            modelpretty12 = sprintf('mean of a normal distribution fitted to the SPF probability bins');
        case 'NMeansAndVars'
            modelpretty12 = sprintf('mean and variance of a normal distribution fitted to the SPF probability bins');
        case 'GBMeansOnly'
            modelpretty12 = sprintf('mean of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansAndVars'
            modelpretty12 = sprintf('mean and variance of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBMeansVarsAndSkew'
            modelpretty12 = sprintf('mean, variance, and skewness of a generalized beta distribution fitted to the SPF probability bins');
        case 'GBVarsSkewAndSPFAnnPoints'
            modelpretty12 = sprintf('variance and skewness of a generalized beta distribution fitted to the SPF probability bins, keeping the annual SPF point forecasts');
        otherwise
            error('ETlabel <<%s>> not recognized', ETlabel12)
    end

    % add ETlabel to modeltype
    modelETtype0 = strcat('ET', ETlabel0, '-', modeltype0);
    modelETtype1 = strcat('ET', ETlabel1, '-', modeltype1);
    % drop "STATE" moniker (meaningful only when VARSTATE is around)
    modelETtype0 = replace(modelETtype0, 'STATE', '');
    modelETtype1 = replace(modelETtype1, 'STATE', '');
    % add ETlabel to modeltype
    modelETtype02 = strcat('ET', ETlabel02, '-', modeltype02);
    modelETtype12 = strcat('ET', ETlabel12, '-', modeltype12);
    % drop "STATE" moniker (meaningful only when VARSTATE is around)
    modelETtype02 = replace(modelETtype02, 'STATE', '');
    modelETtype12 = replace(modelETtype12, 'STATE', '');

    %% allocate memory for forecast stats
    [MAE0, MAE1, RMSE0, RMSE1, CRPS0, CRPS1]        = deal(NaN(Nhorizons, Nvariables, Nsam));
    [relRMSEtstat, relMAEtstat, relCRPStstat]       = deal(NaN(Nhorizons, Nvariables, Nsam));
    [MAE02, MAE12, RMSE02, RMSE12, CRPS02, CRPS12]  = deal(NaN(Nhorizons, Nvariables, Nsam));
    [relRMSEtstat2, relMAEtstat2, relCRPStstat2]    = deal(NaN(Nhorizons, Nvariables, Nsam));

    %% loop over samples
    for s = 1  : length(sam)

        MAXHORIZON = NaN(length(DATALABELS), 1);

        %% collect stats for every variable

        for d = 1 : length(DATALABELS)

            %% prepare things


            datalabel = DATALABELS{d};

            modellabel0  = strcat(datalabel, modeltype0);
            modellabel1  = strcat(datalabel, modeltype1);

            modellabel02  = strcat(datalabel, modeltype02);
            modellabel12  = strcat(datalabel, modeltype12);


            %% load data

            % load model 0
            switch ETlabel0
                case 'none'
                    isET0 = false;
                    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel0, Ndraws0);
                    mat0 = matfile(fullfile(resultdir, matfilename));
                case {'binsOnly', 'binsAndMeans', 'NMeansOnly', 'NMeansAndVars', 'GBMeansOnly', 'GBMeansAndVars','GBMeansVarsAndSkew','GBVarsSkewAndSPFAnnPoints'}
                    isET0 = true;
                    matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel0, Ndraws0, modellabel0);
                    mat0 = matfile(fullfile(ETresultdir, matfilename));
                otherwise
                    error('ETlabel <<%s>> not recognized', ETlabel0)
            end

            % load model 02
            switch ETlabel02
                case 'none'
                    isET02 = false;
                    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel02, Ndraws02);
                    mat02 = matfile(fullfile(resultdir, matfilename));
                case {'binsOnly', 'binsAndMeans', 'NMeansOnly', 'NMeansAndVars', 'GBMeansOnly', 'GBMeansAndVars','GBMeansVarsAndSkew','GBVarsSkewAndSPFAnnPoints'}
                    isET02 = true;
                    matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel02, Ndraws02, modellabel02);
                    mat02 = matfile(fullfile(ETresultdir, matfilename));
                otherwise
                    error('ETlabel <<%s>> not recognized', ETlabel02)
            end
            

            % load model 1
            switch ETlabel1
                case 'none'
                    isET1 = false;
                    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel1, Ndraws1);
                    mat1 = matfile(fullfile(resultdir, matfilename));
                case {'binsOnly', 'binsAndMeans', 'NMeansOnly', 'NMeansAndVars', 'GBMeansOnly', 'GBMeansAndVars','GBMeansVarsAndSkew','GBVarsSkewAndSPFAnnPoints'}
                    isET1 = true;
                    matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel1, Ndraws1, modellabel1);
                    mat1 = matfile(fullfile(ETresultdir, matfilename));
                otherwise
                    error('ETlabel <<%s>> not recognized', ETlabel1)
            end

            % load model 12
            switch ETlabel12
                case 'none'
                    isET12 = false;
                    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel12, Ndraws12);
                    mat12 = matfile(fullfile(resultdir, matfilename));
                case {'binsOnly', 'binsAndMeans', 'NMeansOnly', 'NMeansAndVars', 'GBMeansOnly', 'GBMeansAndVars','GBMeansVarsAndSkew','GBVarsSkewAndSPFAnnPoints'}
                    isET12 = true;
                    matfilename = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s.mat', ETlabel12, Ndraws12, modellabel12);
                    mat12 = matfile(fullfile(ETresultdir, matfilename));
                otherwise
                    error('ETlabel <<%s>> not recognized', ETlabel12)
            end

            Yfuture        = mat0.Yfuture;
            Zdata          = mat0.Zdata;
            dates          = mat0.dates;
            MAXHORIZON(d)  = mat0.Nhorizons;
            theseNhorizons = mat0.Nhorizons;

            T              = length(dates);
            Nz             = size(Zdata, 2);
            Znanny         = mat0.Znanny;
            Ylabel         = mat0.Ylabel;
            Yfcstlabel     = Ylabel(2:end);

            if isET0
                Tstart0 = mat0.firstT;
            else
                Tstart0 = mat0.Tstart;
            end
            if isET02
                Tstart02 = mat02.firstT;
            else
                Tstart02 = mat02.Tstart;
            end
            if isET1
                Tstart1 = mat1.firstT;
            else
                Tstart1 = mat1.Tstart;
            end
            if isET12
                Tstart12 = mat12.firstT;
            else
                Tstart12 = mat12.Tstart;
            end


            Tstart        = max([Tstart0, Tstart02, Tstart1, Tstart12]);

            % checks
            if ~isequal(dates, mat1.dates)
                error('mismatch between input files for model1 and model2 (dates)')
            end
            if ~isequal(mat02.dates, mat12.dates)
                error('mismatch between input files for model12 and model22 (dates)')
            end

            if ~isequaln(Zdata, mat1.Zdata)
                warning('mismatch between input files for model1 and model2 (Zdata, %s)', datalabel)
            end
            if ~isequaln(mat02.Zdata, mat12.Zdata)
                warning('mismatch between input files for model12 and model22 (Zdata, %s)', datalabel)
            end

            if ~isequaln(Yfuture, mat1.Yfuture)
                error('mismatch between input files for model1 and model2 (Yfuture)')
            end
             if ~isequaln(mat02.Yfuture, mat12.Yfuture)
                error('mismatch between input files for model12 and model22 (Yfuture)')
            end
            if Nhorizons < mat1.Nhorizons
                error('Nhorizons does not match assumption')
            end
            if mat02.Nhorizons < mat12.Nhorizons
                error('Nhorizons does not match assumption')
            end
            if mat0.Nhorizons ~= mat1.Nhorizons
                error('mismatch between input files (Nhorizons)')
            end
            if mat02.Nhorizons ~= mat12.Nhorizons
                error('mismatch between input files (Nhorizons)')
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

            %comparisonNote = sprintf('Evaluation window from %s (2009Q2 in the case of UNRATE) through %s (and as far as realized values are available).', ...
            %    datestr(dates(samStart), 'yyyyqq'), datestr(dates(samStop), 'yyyyqq'));

            %% RMSE
            if doMSE
                thistab = 'MSE';
            else
                thistab = 'RMSE';
            end
            tabname = sprintf('%s-%s.tex', thistab, datalabel);

            if doYAVG
                loss0 = (mat0.fcstYAVGhaterror).^2;
                loss1 = (mat1.fcstYAVGhaterror).^2;
                loss02 = (mat02.fcstYAVGhaterror).^2;
                loss12 = (mat12.fcstYAVGhaterror).^2;
            else
                loss0 = (mat0.fcstYhaterror).^2;
                loss1 = (mat1.fcstYhaterror).^2;
                loss02 = (mat02.fcstYhaterror).^2;
                loss12 = (mat12.fcstYhaterror).^2;
            end

            % establish common sample (needed for reporting loss levels; dmtest would catch it)
            nanny        = (isnan(loss0) | isnan(loss1) | isnan(loss02) | isnan(loss12)) | ~sammy; % note: automatic array expansion w.r.t. sammy
            loss0(nanny) = NaN;
            loss02(nanny) = NaN;
            loss1(nanny) = NaN;
            loss12(nanny) = NaN;

            for hh = 1 : theseNhorizons
                if max(abs(loss0(:,hh) - loss1(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                    [~, relRMSEtstat(hh,d,s)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
                else
                    relRMSEtstat(hh,d,s) = 0;
                    %                     warning('virtually identical means (h=%d, %s)', hh, datalabel)
                end
                if max(abs(loss02(:,hh) - loss12(:,hh))) > 1e-6 % for some horizons, means are (virtually) equal
                    [~, relRMSEtstat2(hh,d,s)] = dmtest(loss02(:,hh), loss12(:,hh), hh + 1);
                else
                    relRMSEtstat2(hh,d,s) = 0;
                    %                     warning('virtually identical means (h=%d, %s)', hh, datalabel)
                end
            end

            if doMSE
                RMSE0(1:theseNhorizons,d,s)  = mean(loss0, 1, 'omitnan');
                RMSE1(1:theseNhorizons,d,s)  = mean(loss1, 1, 'omitnan');
                RMSE02(1:theseNhorizons,d,s) = mean(loss02, 1, 'omitnan');
                RMSE12(1:theseNhorizons,d,s) = mean(loss12, 1, 'omitnan');
            else
                RMSE0(1:theseNhorizons,d,s)  = sqrt(mean(loss0, 1, 'omitnan'));
                RMSE1(1:theseNhorizons,d,s)  = sqrt(mean(loss1, 1, 'omitnan'));
                RMSE02(1:theseNhorizons,d,s) = sqrt(mean(loss02, 1, 'omitnan'));
                RMSE12(1:theseNhorizons,d,s) = sqrt(mean(loss12, 1, 'omitnan'));
            end

            
            %% CRPS

            if doYAVG
                loss0  = mat0.fcstYAVGcrps;
                loss1  = mat1.fcstYAVGcrps;
                loss02 = mat02.fcstYAVGcrps;
                loss12 = mat12.fcstYAVGcrps;
            else
                loss0  = mat0.fcstYcrps;
                loss1  = mat1.fcstYcrps;
                loss02 = mat02.fcstYcrps;
                loss12 = mat12.fcstYcrps;
            end

            nanny        = (isnan(loss0) | isnan(loss1) | isnan(loss02) | isnan(loss12)) | ~sammy; % note: automatic array expansion w.r.t. sammy
            loss0(nanny) = NaN;
            loss1(nanny) = NaN;
            loss02(nanny) = NaN;
            loss12(nanny) = NaN;

            for hh = 1 : theseNhorizons
                [~, relCRPStstat(hh,d,s)] = dmtest(loss0(:,hh), loss1(:,hh), hh + 1);
                [~, relCRPStstat2(hh,d,s)] = dmtest(loss02(:,hh), loss12(:,hh), hh + 1);
            end

            CRPS0(1:theseNhorizons,d,s)  = mean(loss0,1, 'omitnan');
            CRPS02(1:theseNhorizons,d,s) = mean(loss02,1, 'omitnan');
            CRPS1(1:theseNhorizons,d,s)  = mean(loss1,1, 'omitnan');
            CRPS12(1:theseNhorizons,d,s) = mean(loss12,1, 'omitnan');

        end % datalabels
    end % samples

        %% tabulate relative RMSE and CRPS in joint table
        if doMSE
            statname1   = 'MSE';
        else
            statname1   = 'RMSE';
        end
        deltaloss1  = RMSE1 ./ RMSE0;
        deltaloss12  = RMSE12 ./ RMSE02;
        deltaTstat1 = relRMSEtstat;
        deltaTstat12 = relRMSEtstat2;
        % TODO: are the following lines still needed ?
        % weed out "significant" relRMSE of 1.00
         ndx              = round(deltaloss1, 2) == 1;
         deltaTstat1(ndx) = 0;
         ndx              = round(deltaloss12, 2) == 1;
         deltaTstat12(ndx) = 0;

        statname2   = 'CRPS';
        deltaloss2  = CRPS1 ./ CRPS0;
        deltaloss22  = CRPS12 ./ CRPS02;
        deltaTstat2 = relCRPStstat;
        deltaTstat22 = relCRPStstat2;
        % TODO: are the following lines still needed ?
        % weed out "significant" relCRPS of 1.00
         ndx              = round(deltaloss2, 2) == 1;
         deltaTstat2(ndx) = 0;
         ndx              = round(deltaloss22, 2) == 1;
         deltaTstat22(ndx) = 0;

        
        
            tabname    = sprintf('relative%s-%sand%s-%s-%s-vs-%s-AND-%s-vs-%s-%s.tex', datatype, statname1, statname2, DATALABELS_CONC, ...
                modelETtype0, modelETtype1,modelETtype02, modelETtype12, [sam(:).label]);
            tabcaption = tabname;
            tabrelstats(tabname, wrap, ...
                statname1, deltaloss1, deltaTstat1, ...
                statname2, deltaloss2, deltaTstat2, ...
                deltaloss12, deltaTstat12, ...
                deltaloss22, deltaTstat22, ...
                modelETtype0, modelETtype1, models(m0).pretty, models(m02).pretty, ...
                DATALABELS, tabcaption, ...
                sam(1).pretty,sam(2).pretty, samstartpretty, samendlongpretty, samendshortpretty, ...
                modelpretty1 ,DATALABELS_pretty);
        
end % PAIRS


%% finish / clean up
finishwrap
finishscript
dockAllFigures



function tabrelstats(tabname, wrap, ...
    statname1, deltaloss1, deltaTstat1, ...
    statname2, deltaloss2, deltaTstat2, ...
    deltaloss12, deltaTstat12,...
    deltaloss22, deltaTstat22,...
    modelETtype0, ~, prettylabel0, prettylabel02, DATALABELS, tabcaption, ...
    sam1, sam2, samstartpretty, samendlongpretty, samendshortpretty, ...
    modelpretty1, DATALABELS_pretty) 


%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
else
    tabdir = fullfile(localtemp, 'foo');
end

Nhorizons_vec = sum(~isnan(deltaloss1(:,:,1)));

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, 1 + 8));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 4, statname1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 4, statname2);
fprintf(fid, '\\\\ ');
fprintf(fid, '\n');
fprintf(fid, '\\cmidrule(lr){%d-%d}', 2, 1+4);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 1+4+1, 1+8);
fprintf(fid, '\n');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 2, prettylabel0);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 2, prettylabel02);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 2, prettylabel0);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 2, prettylabel02);
fprintf(fid, '\\\\\n');
fprintf(fid, '\\cmidrule(lr){%d-%d}', 2, 3);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 4, 5);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 6, 7);
fprintf(fid, '\\cmidrule(lr){%d-%d}', 8, 9);
fprintf(fid, '\n');
fprintf(fid, '$h$');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam2);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam2);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam2);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam1);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, sam2);

fprintf(fid, '\\\\\n');


panID = {'A','B'};

for nn = 1 : length(DATALABELS)
fprintf(fid, '\\midrule\n');
if length(DATALABELS) > 1
    fprintf(fid, ' & \\multicolumn{%d}{c}{%s} \\\\', 8, append('Panel ',panID{nn},': ',DATALABELS_pretty{nn}));
    fprintf(fid, '\\midrule\n');
end
for hh = 1 : Nhorizons_vec(nn)
    fprintf(fid, '%d ', hh-1);
    % print relative RMSE of model 1 in sam1, sam2
    for s = 1 : 2
        if isnan(deltaloss1(hh,nn,s))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss1(hh,nn,s),  Zstar(deltaTstat1(hh,nn,s)));
        end
    end
    % print relative RMSE of model 2 in sam1, sam2
    for s = 1 : 2
        if isnan(deltaloss12(hh,nn,s))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss12(hh,nn,s),  Zstar(deltaTstat12(hh,nn,s)));
        end
    end
    % print relative CRPS of model 1
    for s = 1 : 2
        if isnan(deltaloss2(hh,nn,s))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss2(hh,nn,s),  Zstar(deltaTstat2(hh,nn,s)));
        end
    end
     % print relative CRPS of model 2
    for s = 1 : 2
        if isnan(deltaloss22(hh,nn,s))
            fprintf(fid, '& \\multicolumn{1}{c}{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaloss22(hh,nn,s),  Zstar(deltaTstat22(hh,nn,s)));
        end
    end
    fprintf(fid, '\\\\\n');
end

end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


fprintf(fid, 'Note: \n');
if ~contains(modelETtype0,'none')
    fprintf(fid, 'Relative %s and %s of forecasts obtained from tilting to moments imputed from the SPF probability bins vs. tilting to the bins directly (with the latter in the denominator). \n',statname1, statname2);
else
    fprintf(fid, 'Relative %s and %s of forecasts obtained from tilting to moments imputed from the SPF probability bins vs. CGM model output (with the latter in the denominator). \n',statname1, statname2);
end
if strcmpi(modelpretty1(6:7),'of') % singular
    fprintf(fid, 'The imputed moment is %s. \n',modelpretty1);
else % plural
    fprintf(fid, 'The imputed moments are %s. \n',modelpretty1);
end

fprintf(fid, 'Forecasts are for quarterly outcomes and evaluation windows extend from %s until %s and %s, respectively (using data for realized values as far as available in %s). \n', samstartpretty, samendlongpretty, samendshortpretty,samendlongpretty);
fprintf(fid, 'Significance assessed by Diebold-Mariano tests using Newey-West standard errors with $h + 1$ lags.\n');
fprintf(fid, '%s\n', '$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\%, and 10\% level, respectively.');

% fprintf(fid, append('Forecasts for quarterly '));
% if length(DATALABELS_pretty) == 1
%     fprintf(fid, DATALABELS_pretty{1});
% elseif length(DATALABELS_pretty) == 2
%     fprintf(fid, append(DATALABELS_pretty{1}, ' and ', DATALABELS_pretty{2}));
% elseif length(DATALABELS_pretty) > 2
%     for i = 1 : length(DATALABELS_pretty) - 1
%         fprintf(fid, append(DATALABELS_pretty{i}, ', '));
%     end
%     fprintf(fid, append(' and ',DATALABELS_pretty{i}));
% end

%fprintf(fid, append([' $h$ steps ahead ' ...
%    'over subsamples extending from 1992Q1 until 2022Q2 and 2016Q4, respectively (using data for realized values as far as available in 2022Q2). \n']));
%fprintf(fid, 'Relative %s and %s of ``%s'''' vs. ``%s'''' (in denominator) and ``%s'''' vs. ``%s'''' (in denominator).\n', ...
%    statname1, statname2, modelpretty1, modelpretty0, modelpretty12, modelpretty02);

%fprintf(fid, '%s \n', comparisonNote);
%fprintf(fid, 'Variable mnemonics: RGDP denotes real growth, PGDP changes in the GDP deflator and UNRATE the unemployment rate. (All growth rates are expressed as annualized percentage points of quarterly rates of change.)\n');



fclose(fid);
type(fullfile(tabdir, tabname))

end
