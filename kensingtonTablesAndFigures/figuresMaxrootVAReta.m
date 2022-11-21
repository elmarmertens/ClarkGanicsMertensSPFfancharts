%% compute and plot bias in SPF-consistent expectations of VAReta

% todo:
% - load predictive means directly (once stored with MCMC)
% - in MCMC: compute predictive means analytically (instead of MC)
% - compute bias directly in MCMC to generate also uncertainty around bias

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


%% some parameters


modeltype   = 'VARprior10STATEtrendgapSV';
modelpretty = 'non-MDS';
% modeltype   = 'VARprior5STATEtrendgapSV';
% modelpretty = 'non-MDS (tighter prior SV)';
% modeltype   = 'VARprior5STATEconst';
% modelpretty = 'non-MDS (tighter prior; CONST)';
% Ndraws      = 3e3;
% modeltype   = 'VARprior50STATEtrendgapSV';
% modelpretty = 'non-MDS (loose prior; SV)';
Ndraws      = 3e3;

resultdir = localresultsMCMC;
resultdir = '~/jam/lager/KENSINGTON/kensingtonresults3';

fontsize = 18;

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};

%% initwrap
wrap = [];
titlename = sprintf('SPFgeigenvalues-%s', modeltype);
initwrap

%% loop over datalabels

for d = 1 : length(DATALABELS)

    datalabel   = DATALABELS{d};
    modellabel  = strcat(datalabel, modeltype);

    %% prepare latexwrapper
    %     wrap = [];
    %     titlename = sprintf('SPFbias-%s-%s', datalabel, modeltype);
    %     initwrap

    fprintf('Processing %s ... \n', datalabel)

    %% load data
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    %     if doMedian
    %         Yhat  = mat.fcstYmedian;
    %     else
    %         Yhat  = mat.fcstYhatRB;
    %     end
    %     Ydraws        = mat.Ydraws;

    %     Yfuture       = mat.Yfuture;
    %     Nz            = mat.Nz;
    %     Zdata         = mat.Zdata;
    %     Znanny        = mat.Znanny;

    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    Tstart        = mat.Tstart;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    horizonLabels = arrayfun(@(x) sprintf('h=%d', x), 0:Nhorizons-1, 'uniformoutput', false);


    %% plot G-eigenvalues
    Gpriorlambda = mat.Gpriorlambda;
    Glambda      = mat.Glambda;

    thisfig = figure;
    hold on
    set(gca, 'FontSize', fontsize)
    hprior = plot(dates, median(Gpriorlambda,1), 'k--', 'linewidth', 2);
    plot(dates, prctile(Gpriorlambda, [25 75], 1), 'k--')

    hposterior = plot(dates, median(Glambda,1), 'b-', 'linewidth', 2);
    plot(dates, prctile(Glambda, [25 75], 1), 'b-')

    % ylim([min([ylim, 0]) max(ylim)])
    ylim([0 1])
    yticks(0 : .2 : 1)
    xtickdates(dates([Tstart T]))
    wrapthisfigure(thisfig, sprintf('Geigenvalues-%s', modellabel), wrap)
    legend([hprior hposterior], 'Prior', 'Posterior')
    wrapthisfigure(thisfig, sprintf('Geigenvalues-%s-WITHLEGEND', modellabel), wrap)
    
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
