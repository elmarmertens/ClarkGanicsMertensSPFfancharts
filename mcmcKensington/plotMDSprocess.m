%% MDS vs VAR0 model -- single run for specific thisT
% with constants added to each gap equation

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
addpath ../matlabtoolbox/empbsbox/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*PFBNS>
%#ok<*UNRCH>
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc
rng('default')

%% some parameters

rndStream  = getDefaultStream;

quicky     = false;  % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'};

% select parameters
thisDate       = datenum(2019,10,1);
doSamStartSPF  = false;

NGAP1    = 'BOP';
NGAP2    = 'BOP';

LABEL1   = 'MDS';
LABEL2   = 'VAR';

thinSTEPS = 1;



%% some settings

fontsize = 18;

color1   = colors4plots("darkblue");
color2   = colors4plots("darkred");

resultsdir     = '~/jam/fooMDSvVAR';

startDateLabel = '1968Q4';
thisDateLabel = datestr(thisDate, 'yyyyqq');

%% loop over SV choices
for SVCHOICE = {'mean', 'median', 'min', 'max'}
    SVchoice = SVCHOICE{1};

    titlename = sprintf('plotMDSprocess-%sSV', SVchoice);
    initwrap

    %% allocate memory to collect draws for all datalabels
    Hstar          = 20; % horizons to plot
    MCMCdraws      = 3e3;
    NGAP           = 17;
    MDSprocess     = zeros(Hstar, MCMCdraws, length(DATALABELS));
    MDSkalmanTilde = NaN(NGAP, MCMCdraws, length(DATALABELS));
    MDSkalmanStar  = NaN(1, MCMCdraws, length(DATALABELS));

    %% loop over datalabels
    for d = 1 : length(DATALABELS)

        close all
        %% load dataset for variable
        datalabel = DATALABELS{d};

        datadir = fullfile('..', 'matdataKensington');
        modellabel = sprintf('%s-MDSvsVAR0trendHScycleSV2blockNoiseHS-Ngap1%s-Ngap2%s-samStart%s-%s-thin%d', ...
            datalabel, NGAP1, NGAP2, startDateLabel, thisDateLabel, thinSTEPS);
        matfilename = fullfile(resultsdir, strcat(modellabel, '.mat'));
        mat = matfile(matfilename);

        % MDSkalmanStar = mat.MDSkalmanStar;
        % MDSprocess    = mat.MDSprocess;
        % Hstar         = mat.Hstar;

        Ngap1         = mat.Ngap1;
        % MCMCdraws     = mat.MCMCdraws;
        SVtdraws1     = mat.SVtdraws1;
        sqrtsigmaDraws1 = mat.sqrtsigmaDraws1;
        mustarvardraws1 = mat.mustarvardraws1;

        %% set modellabel for outputs
        modellabel = sprintf('%s-MDStrendHScycleSV2blockNoiseHS-Ngap%s-samStart%s-%s-thin%d-SV%s', ...
            datalabel, NGAP1, startDateLabel, thisDateLabel, thinSTEPS, SVchoice);

        switch SVchoice
            case {'min', 'max'}
                SVfun = str2func(sprintf('@(x,dim) %s(x,[],dim)', SVchoice));
            otherwise
                SVfun = str2func(SVchoice);
        end

        %% plot SVt paths
        dates    = mat.datesT;
        midSV    = median(SVtdraws1, 3);
        tailSV   = prctile(SVtdraws1, normcdf11, 3);
        midSVfun = SVfun(midSV,2);
        thisfig = figure;
        for nn = 1 : 2
            subplot(2,1,nn)
            hold on
            set(gca, 'fontsize', fontsize)
            plot(dates,midSV(nn,:), '-', 'color', colors4plots('blue'), 'LineWidth', 2)
            plot(dates,squeeze(tailSV(nn,:,:)), '-', 'color', colors4plots('blue'), 'LineWidth', 1)
            hfun = yline(midSVfun(nn), ':', 'color', colors4plots('lightblue'), 'LineWidth', 2);
            plotOrigin
            xtickdates(dates)
            legend(hfun, sprintf('%s along median path', SVchoice), 'location', 'best')
            title(sprintf('SV block %d', nn))
        end
        sgtitle(sprintf('%s: SV-t paths', datalabel))
        wrapthisfigure(thisfig, sprintf('SVt-%s', modellabel), wrap)

        %% MDS: construct univariate process for y(t)
        Ngap      = Ngap1; % "H" in the notes
        ndxSV11   = 1:5;
        ndxSV22   = ndxSV11(end)+1:Ngap;

        % prepare state space matrices
        PSItilde       = diag(ones(Ngap-1,1),1);
        AA             = blkdiag(PSItilde, 1);
        BB             = zeros(Ngap+1);
        cc             = zeros(1,Ngap+1); cc(1) = 1; cc(end) = 1;


        for mm = 1 : MCMCdraws

            thisSV              = SVfun(SVtdraws1(:,:,mm),2);

            volTrend         = sqrt(mustarvardraws1(mm));

            sqrtSigmaTilde            = sqrtsigmaDraws1(:,:,mm);
            sqrtSigmaTilde(:,ndxSV11) = sqrtSigmaTilde(:,ndxSV11) .* thisSV(1);
            sqrtSigmaTilde(:,ndxSV22) = sqrtSigmaTilde(:,ndxSV22) .* thisSV(2);

            BB(1:Ngap,1:Ngap) = sqrtSigmaTilde;
            BB(end,end)       = volTrend;
            % Compute steady-state gain
            [~, K]                      = abckalman(AA,BB,cc);
            MDSkalmanTilde(1:Ngap,mm,d) = K(1:Ngap);
            MDSkalmanStar(1,mm,d)       = K(end);
        end

        MDSprocess(1:Ngap,:,d)   = MDSkalmanTilde(1:Ngap,:,d);
        MDSprocess(:,:,d)        = MDSprocess(:,:,d) + MDSkalmanStar(:,:,d);
        MDSprocess(1,:,d)        = 1;

        %% plot MDS process
        midstar   = median(MDSkalmanStar(:,:,d),2);
        midstar68 = prctile(MDSkalmanStar(:,:,d), normcdf11,2);
        mid       = median(MDSprocess(:,:,d),2);
        tail68    = prctile(MDSprocess(:,:,d), normcdf11, 2);
        thisfig   = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        hMA   = plot(0:Hstar-1, mid, '-', 'color', colors4plots('blue'), 'LineWidth', 2);
        hMA68 = plot(0:Hstar-1, tail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
        hStar = yline(midstar, ':', 'color', colors4plots('lightblue'), 'LineWidth', 3);
        yline(midstar68, ':', 'color', colors4plots('lightblue'), 'LineWidth', 1);
        plotOrigin
        xlim([0 Hstar-1])
        wrapthisfigure(thisfig, sprintf('MDSprocess-%s', modellabel), wrap, [], [], [], [], true);
        legend([hMA,hStar], 'IMA process', 'Endpoint shift', 'location', 'best')
        wrapthisfigure(thisfig, sprintf('MDSprocess-%s-WITHLEGEND', modellabel), wrap)

    end

    %% compare IRF for pairs of variables
    for d1 = [1 3]
        for d2 = 3 : 4
            if d1 == d2
                continue
            end
            mid1 = median(MDSprocess(:,:,d1),2);
            tail1 = prctile(MDSprocess(:,:,d1), normcdf11, 2);
            mid2 = median(MDSprocess(:,:,d2),2);
            tail2 = prctile(MDSprocess(:,:,d2), normcdf11, 2);
            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)
            hMA1 = plot(0:Hstar-1, mid1, '-', 'color', colors4plots('blue'), 'LineWidth', 2);
            plot(0:Hstar-1, tail1, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
            hMA2 = plot(0:Hstar-1, mid2, '--', 'color', colors4plots('red'), 'LineWidth', 2);
            plot(0:Hstar-1, tail2, '--', 'color', colors4plots('red'), 'LineWidth', 1);
            plotOrigin
            xlim([0 Hstar-1])
            wrapthisfigure(thisfig, sprintf('MDSprocess-%s-vs-%s', DATALABELS{d1}, DATALABELS{d2}), wrap, [], [], [], [], true);
            legend([hMA1,hMA2], DATALABELS{d1}, DATALABELS{d2}, 'location', 'best')
            wrapthisfigure(thisfig, sprintf('MDSprocess-%s-vs-%s-WITHLEGEND', DATALABELS{d1}, DATALABELS{d2}), wrap);
        end
    end

    dockAllFigures
    finishwrap
end % SVchoices

%% finish / clean up
finishscript
