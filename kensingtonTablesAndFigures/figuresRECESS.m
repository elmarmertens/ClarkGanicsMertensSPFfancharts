%% calculate recession probabilities and compare against SPF measures

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

fontsize = 18;

doSince1992 = true;

%% some parameters


modeltype   = 'STATEtrendgapSV';
modelpretty = 'SV';
% modeltype   = 'STATEconst';
% modelpretty = 'CONST';
Ndraws      = 3e3;

ETlabel     = 'binsOnly';

resultdir = localresultsMCMC;

fontsize  = 18;

datalabel = 'RGDP';

doTitle = false;

%% set up wrapper
titlename = sprintf('recess-%s', modeltype);
wrap = [];
initwrap


%% load data
modellabel    = strcat(datalabel, modeltype);
matfilename   = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
mat           = matfile(fullfile(resultdir, matfilename));
matfilename   = sprintf('kensington-EVAL-ET%s-Ndraws%d-%s%s',ETlabel,Ndraws,datalabel,modeltype);
matET         = matfile(fullfile(resultdir, matfilename));

Yfuture        = mat.Yfuture;

Nz            = mat.Nz;
Zdata         = mat.Zdata;
Znanny        = mat.Znanny;

Nhorizons     = mat.Nhorizons;
dates         = mat.dates;
T             = mat.T;
Tstart        = mat.Tstart;
datesQ        = mat.datesQ;
doNIPA        = mat.doNIPA;

fcstYrecess   = mat.fcstYrecess;
fcstYrecessET = matET.fcstYrecess;


if doSince1992
    Tstart = find(dates == datenum(1992,1,1));
end

%% load SPF RECESS
spfDATA   = importdata('../kensingtonDataUSSPF/Mean_RECESS_Level.xlsx');
spfdates  = datenum(spfDATA.data(:,1), (spfDATA.data(:,2)-1)*3+1, 1);
spfRECESS = spfDATA.data(:,3:end);

ndx       = ismember(spfdates, dates);
spfRECESS = spfRECESS(ndx,:);
spfdates  = spfdates(ndx);

H = 5;
if size(spfRECESS, 2) ~= H
    error('dimension mismatch')
end

%% plot RECESS prob
for hh = 1 : H
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hSPF   = plot(spfdates, spfRECESS(:,hh), '-', 'color', Colors4Plots(8), 'LineWidth', 5);
    hModel = plot(dates, fcstYrecess(:,hh),   '-.', 'color', Colors4Plots(7), 'LineWidth', 3);
    hET    = plot(dates, fcstYrecessET(:,hh), ':',  'color', Colors4Plots(3), 'LineWidth', 3);
    % xtickdates(dates(Tstart:T))
    nbershades(dates(Tstart:T))
    ylim([0 max(ylim)])
    wrapthisfigure(thisfig, sprintf('RECESSh%d-%s', ...
        hh - 1, modeltype), wrap);
    legend([hModel hET hSPF], modelpretty, sprintf('ET (%s)', modelpretty), 'SPF')
    if doTitle
        title(sprintf('h = %d', hh - 1)) 
    end
    wrapthisfigure(thisfig, sprintf('RECESSh%d-%s-WITHLEGEND', ...
        hh - 1, modeltype), wrap);
end

%% scale recess probs to (0,1)
spfRECESS     = spfRECESS ./ 100;
fcstYrecess   = fcstYrecess ./ 100;
fcstYrecessET = fcstYrecessET ./ 100;

%% compute DRPS
spfDRPS = NaN(T,H);

isRECESS = Yfuture(:,1:H) < 0;

spfDRPS  = (spfRECESS - isRECESS).^2 + ((1 - spfRECESS) - 1).^2;
% check simplified formula
checkdiff(spfDRPS, (spfRECESS - isRECESS).^2 + spfRECESS.^2);

modelDRPS    = (fcstYrecess - isRECESS).^2 + fcstYrecess.^2;
ETmodelDRPS  = (fcstYrecessET - isRECESS).^2 + fcstYrecessET.^2;

%% collected losses into matrix
LOSS = cat(3, spfDRPS, modelDRPS, ETmodelDRPS);

% cut to eval window
LOSS       = LOSS(Tstart:T,:,:);
thesedates = dates(Tstart:T);

Nmodels = size(LOSS,3);

%% plot DRPS
linestyles = {'-', '-.', ':'};
colorlist  = {Colors4Plots(8), Colors4Plots(7), Colors4Plots(3)};
linewidths = {5, 3, 3};

for hh = 1 : H
    hanni = NaN(Nmodels,1);
    thisfig = figure;
    hold on
    set(gca, 'FontSize', fontsize)
    for mm = 1 : Nmodels
        hanni(mm) = plot(thesedates, LOSS(:,hh,mm), 'color', colorlist{mm}, 'LineStyle', linestyles{mm}, 'LineWidth', linewidths{mm});
    end
    nbershades(thesedates)
    if doTitle
        title(sprintf('h = %d (DRPS)', hh - 1))
    end
    legend(hanni, 'SPF', modelpretty, sprintf('ET (%s)', modelpretty), 'location', 'best')
    wrapthisfigure(thisfig, sprintf('DRPS-RECESSh%d-%s', ...
        hh - 1, modeltype), wrap);
end

%% plot cum-avg DRPS

avgLOSS    = cumsum(LOSS, 1) ./ (1:size(LOSS,1))';

linestyles = {'-', '-.', ':'};
colorlist  = {Colors4Plots(8), Colors4Plots(7), Colors4Plots(3)};
linewidths = {5, 3, 3};

for hh = 1 : H
    hanni = NaN(Nmodels,1);
    thisfig = figure;
    hold on
    set(gca, 'FontSize', fontsize)
    for mm = 1 : Nmodels
        hanni(mm) = plot(thesedates, avgLOSS(:,hh,mm), 'color', colorlist{mm}, 'LineStyle', linestyles{mm}, 'LineWidth', linewidths{mm});
    end
    nbershades(thesedates)
    if doTitle
        title(sprintf('h = %d (avgDRPS)', hh - 1))
    end
    wrapthisfigure(thisfig, sprintf('avgDRPS-RECESSh%d-%s', ...
        hh - 1, modeltype), wrap);
    legend(hanni, 'SPF', modelpretty, sprintf('ET (%s)', modelpretty), 'location', 'best')
    wrapthisfigure(thisfig, sprintf('avgDRPS-RECESSh%d-%s-WITHLEGEND', ...
        hh - 1, modeltype), wrap);
end


%% finish / clean up
finishwrap
finishscript
dockAllFigures
