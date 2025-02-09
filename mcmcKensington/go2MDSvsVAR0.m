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
thisDate       = datenum(2024,1,1);
doSamStartSPF  = false;

NGAP1    = 'BOP';
NGAP2    = 'BOP';

LABEL1   = 'MDS';
LABEL2   = 'VAR';

thinSTEPS = 1;

doStore = true;

%% general parameters
samStart = []; % 52 is 1981Q3

%% batch parameters
% SED-PARAMETERS-HERE
thisDate   = datenum(2019,10,1);

%% some settings

fontsize = 24;

color1   = colors4plots("darkblue");
color2   = colors4plots("darkred");

if quicky
    MCMCdraws    = 1e2;
    Nfedraws     = 10;
else
    MCMCdraws    = 3e3;
    Nfedraws     = 100;
end

%% load SEP spreadsheet (for fan charts)
SEPfcstTable    = readtable('../rawdataKensington/SEPforecasterrors_absvalues.xlsx','Sheet','SEP forecasts');
SEPfcstDates    = datenum(SEPfcstTable.Var1, 'yyyy:QQ');
SEPfcstVarnames = SEPfcstTable.Properties.VariableNames;

%% loop over datalabels
for d = 1 : length(DATALABELS)
    
    close all
    datalabel = DATALABELS{d};
    
    
    
    %% load dataset for variable
    
    datadir = fullfile('..', 'matdataKensington');
    matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
    load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', 'Nbar', ...
        'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czhat', 'Czbar')
    
    ndxFixedHorizon     = 1 + 1 + (0:4);
    
    
    % select start of sample and cut sample
    if isempty(samStart)
        if doSamStartSPF
            samStart   = find(sum(~Znanny,2) > 1,1); % first availability of SPF
        else
            samStart   = find(any(~Znanny,2),1); % for CPI/TBILL: leaves early sample with only Yrealized ...
        end
    end
    startDateLabel = datestr(dates(samStart), 'yyyyqq');
    dates          = dates(samStart:end); datesQ = datesQ(samStart:end);
    T              = length(dates);
    Zdata      = Zdata(samStart:end,:);
    Znanny     = Znanny(samStart:end,:);
    Cz         = Cz(:,:,samStart:end);
    ndxZbar    = cat(1, false(Nz-Nbar,1), true(Nbar,1));
    
    thisT    = find(dates <= thisDate, 1, 'last');
    thisDate = dates(thisT);
    thisDateLabel = datestr(thisDate, 'yyyyqq');
    
    % process NGAP1
    [Ngap1, NGAP1] = switchNGAP(NGAP1, Ny, datalabel, thisDate);
    
    
    % process NGAP2
    [Ngap2, NGAP2] = switchNGAP(NGAP2, Ny, datalabel, thisDate);
    
    modellabel = sprintf('%s-MDSvsVAR0trendHScycleSV2blockNoiseHS-Ngap1%s-Ngap2%s-samStart%s-%s-thin%d', ...
        datalabel, NGAP1, NGAP2, startDateLabel, thisDateLabel, thinSTEPS);
    
    
    wrap = [];
    titlename = modellabel;
    if ~quicky
        initwrap
    end
    
    fprintf('Processing %s ... \n', modellabel)
    
    %% load SEP forecast data (to produce fan charts later)
    switch datalabel
        case 'RGDP'
            varNdx = contains(SEPfcstVarnames, 'GDP');
        case 'UNRATE'
            varNdx = contains(SEPfcstVarnames, 'UNRATE');
        case 'CPI'
            varNdx = contains(SEPfcstVarnames, 'PCEINFL');
            warning('Matching SPF for CPI against SEP for PCE')
        otherwise
            warning('no SEP for <<%s>>', datalabel)
            varNdx = [];
    end
    if isempty(varNdx)
        SEPfcst = [];
    else
        SEPfcst = SEPfcstTable{:,varNdx};
    end
    
    %% call MCMC sampler Ngap1
    modelname1 = LABEL1; % sprintf('Ngap%d', Ngap1);
    fprintf('\n');
    OK = false;
    while ~OK
        try
            fprintf('STARTED mcmc1  ...\n')
            tic
            warning('off', 'backtrace');
            [YdensityDraws1, Ydraws1, ETAdraws1, RESIDSVdraws1, ...
                SVdraws1, SVtdraws1, hvcvdraws1, hrhodraws1, tdofdraws1, tScaleDraws1, sqrtsigmaDraws1, ~, mustarvardraws1, ...
                NOISEdraws1, NOISEvarstate1, NOISEsig1, ...
                ] = ...
                mcmcsamplerMDStrendHScycleSVt2blockNoiseHS(Zdata, Znanny, Cz, ...,
                Ngap1, ...
                ndxZbar, dates, ...
                thisT, MCMCdraws, Nfedraws, rndStream, true, thinSTEPS);
            toc
            warning('on', 'backtrace');
            fprintf('DONE mcmc1.\n')
            OK = true;
        catch mcmcME
            warning('mcmc1 aborted at thisT=%d, message: %s', thisT, getReport(mcmcME, 'basic'))
        end
    end
    
    
    Nhorizons1 = size(YdensityDraws1,1);
    Nstates1   = Ngap1 + 1;
    
    %% collect PSRF and INFF for MCMC1
    % Ydraws1
    thesedraws          = Ydraws1;
    % NaN out the known values
    for ii = 1 : 6
        ndx = Znanny(1:thisT,ii);
        thesedraws(ii,~ndx,:) = NaN;
    end
    thesedraws = transpose(reshape(thesedraws, Nstates1 * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: Ydraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: Ydraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    % SV
    Nsv = 2;
    thesedraws = transpose(reshape(SVdraws1, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    % SVt
    thesedraws = transpose(reshape(SVtdraws1, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVtdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVtdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    % tdofdraws
    thesedraws = transpose(tdofdraws1);
    plotPSRF(thesedraws, sprintf('%s: tdofdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: tdofdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
    
    %% MDS: construct univariate process for y(t)
    Ngap      = Ngap1; % "H" in the notes
    ndxSV11   = 1:5;
    ndxSV22   = ndxSV11(end)+1:Ngap;
    
    % prepare state space matrices
    PSItilde       = diag(ones(Ngap-1,1),1);
    AA             = blkdiag(PSItilde, 1);
    BB             = zeros(Ngap+1);
    cc             = zeros(1,Ngap+1); cc(1) = 1; cc(end) = 1;
    
    
    MDSkalmanTilde = NaN(Ngap, MCMCdraws);
    MDSkalmanStar  = NaN(1, MCMCdraws);
    % maxlogESV      = 10; % max vol of about 150
    for mm = 1 : MCMCdraws
        
        ESV              = median(SVtdraws1(:,:,mm),2);
        
        volTrend         = sqrt(mustarvardraws1(mm));
        
        sqrtSigmaTilde            = sqrtsigmaDraws1(:,:,mm);
        sqrtSigmaTilde(:,ndxSV11) = sqrtSigmaTilde(:,ndxSV11) .* ESV(1);
        sqrtSigmaTilde(:,ndxSV22) = sqrtSigmaTilde(:,ndxSV22) .* ESV(2);
        
        BB(1:Ngap,1:Ngap) = sqrtSigmaTilde;
        BB(end,end)       = volTrend;
        % Compute steady-state gain
        [~, K]                 = abckalman(AA,BB,cc);
        MDSkalmanTilde(:,mm)   = K(1:Ngap);
        MDSkalmanStar(mm)      = K(end);
    end
    
    Hstar                    = 20; % horizons to plot
    MDSprocess               = zeros(Hstar, MCMCdraws);
    MDSprocess(1:Ngap,:)     = MDSkalmanTilde;
    MDSprocess               = MDSprocess + MDSkalmanStar;
    MDSprocess(1,:)          = 1;
    
    %% plot MDS process
    midstar   = median(MDSkalmanStar);
    midstar68 = prctile(MDSkalmanStar, normcdf11);
    mid       = median(MDSprocess,2);
    tail68    = prctile(MDSprocess, normcdf11, 2);
    thisfig   = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hMA   = plot(0:Hstar-1, mid, '-', 'color', colors4plots('blue'), 'LineWidth', 2);
    hMA68 = plot(0:Hstar-1, tail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
    hStar = yline(midstar, ':', 'color', colors4plots('orange'), 'LineWidth', 3);
    yline(midstar68, '-', 'color', colors4plots('orange'), 'LineWidth', 1);
    plotOrigin
    xlim([0 Hstar-1])
    wrapthisfigure(thisfig, sprintf('MDSprocess-%s', modellabel), wrap)
    legend([hMA,hStar], 'IMA process', 'endpoint shift', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('MDSprocess-%s-WITHLEGEND', modellabel), wrap)
    
    %% prepare VAR prior
    Ng2    = Ngap2 * Ngap2;
    G0var2 = NaN(Ng2,1); % diagonal elements
    Gtheta = [.2^2, 0.5^2, 3];
    % Gtheta = [.1^2, 0.5^2, 3];
    ndx = 0;
    for i = 1 : Ngap2
        for j = 1 : Ngap2
            ndx  = ndx + 1;
            if (j==i)
                G0var2(ndx)=Gtheta(1);
                % G0var(ndx) = G0var(ndx) / i.^Gtheta(3);
            else
                G0var2(ndx) =Gtheta(1)*Gtheta(2); % /(abs(i-j).^Gtheta(3));
                % G0var(ndx) = G0var(ndx) / (i * abs(i-j)).^Gtheta(3);
                % G0var2(ndx) = G0var2(ndx) / abs(i-j).^Gtheta(3);
            end
            % G0var(ndx) = G0var(ndx) / (i.^Gtheta(3));
        end
    end
    
    %% call MCMC sampler 2
    modelname2 = LABEL2; % sprintf('Ngap%d', Ngap1);
    fprintf('\n');
    
    OK = false;
    while ~OK
        try
            fprintf('STARTED mcmc2 ...\n')
            tic
            warning('off', 'backtrace');
            [YdensityDraws2, Ydraws2, YhatRB2, ETAdraws2, RESIDSVdraws2, ZRESIDdraws2, ...
                Gdraws2, Glambda2, ~, YgapCONSTdraws2, ...
                SVdraws2, SVtdraws2, hvcvdraws2, hrhodraws2, tdofdraws2, tScaleDraws2, sqrtsigmaDraws2, ~, mustarvardraws2, ...
                NOISEdraws2, NOISEvarstate2, ...
                ] = ...
                mcmcsamplerVAR0trendHScycleSVt2blockNoiseHS(Zdata, Znanny, Cz, ...,
                Ngap2, ...
                G0var2, ...
                ndxZbar, dates, ...
                thisT, MCMCdraws, Nfedraws, rndStream, true, thinSTEPS);
            toc
            warning('on', 'backtrace');
            fprintf('DONE mcmc2.\n')
            OK = true;
        catch mcmcME
            warning('mcmc1 aborted at thisT=%d, message: %s', thisT, getReport(mcmcME, 'basic'))
        end
    end
    Nhorizons2 = size(YdensityDraws2,1);
    Nstates2   = Ngap2 + 1;
    
    
    %% collect PSRF and INFF for MCMC2
    % Ydraws2
    thesedraws          = Ydraws2;
    % NaN out the known values
    for ii = 1 : 6
        ndx = Znanny(1:thisT,ii);
        thesedraws(ii,~ndx,:) = NaN;
    end
    thesedraws = transpose(reshape(thesedraws, Nstates2 * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: Ydraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: Ydraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    % SV
    Nsv = 2;
    thesedraws = transpose(reshape(SVdraws2, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    % SVt
    thesedraws = transpose(reshape(SVtdraws2, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVtdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVtdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    % tdofdraws
    thesedraws = transpose(tdofdraws2);
    plotPSRF(thesedraws, sprintf('%s: tdofdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: tdofdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    % Gdraws2
    thesedraws = transpose(reshape(Gdraws2, Ngap2 * Ngap2, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: Gdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: Gdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    % YgapCONSTdraws2
    thesedraws = transpose(YgapCONSTdraws2);
    plotPSRF(thesedraws, sprintf('%s: YgapCONSTdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: YgapCONSTdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
    
    %% VAR: construct univariate process for y(t)
    Ngap      = Ngap2; % "H" in the notes
    Nstates   = Ngap * 2 + 1; % VAR(2) in gaps plus trend
    ndxSV11   = 1:5;
    ndxSV22   = ndxSV11(end)+1:Ngap;
    
    % prepare state space matrices (trend at bottom)
    PSItilde       = diag(ones(Ngap-1,1),1);
    AA             = zeros(Nstates);
    AA(Ngap+(1:Ngap),(1:Ngap)) = eye(Ngap); % kompanion part
    AA(end,end)    = 1; % trend
    BB             = zeros(Nstates, Ngap + 1);
    cc             = zeros(1,Nstates); cc(1) = 1; cc(end) = 1;
    CC             = [eye(Ngap) zeros(Ngap) ones(Ngap,1)]; % to select SPF term structure
    
    
    Hstar          = 20; % horizons to plot
    VARprocess     = zeros(Hstar, MCMCdraws);
    VARkalmanStar  = NaN(1, MCMCdraws);
    SPFresponses   = zeros(Ngap,Hstar,MCMCdraws);
    
    % maxlogESV      = 10; % max vol of about 150
    for mm = 1 : MCMCdraws
        
        PItilde             = Gdraws2(:,:,mm);
        OMEGA1              = PItilde + PSItilde;
        OMEGA2              = - PItilde * PSItilde;
        AA(1:Ngap,1:Ngap*2) = [OMEGA1, OMEGA2];
        
        % VCV
        ESV                       = median(SVtdraws2(:,:,mm),2);
        
        sqrtSigmaTilde            = sqrtsigmaDraws2(:,:,mm);
        sqrtSigmaTilde(:,ndxSV11) = sqrtSigmaTilde(:,ndxSV11) .* ESV(1);
        sqrtSigmaTilde(:,ndxSV22) = sqrtSigmaTilde(:,ndxSV22) .* ESV(2);
        
        BB(1:Ngap,1:Ngap) = sqrtSigmaTilde;
        
        % trend shocks
        volTrend        = sqrt(mustarvardraws2(mm));
        BB(end,end)     = volTrend;
        
        % Compute steady-state gain
        [~, K]                 = abckalman(AA,BB,cc);
        VARkalmanStar(mm)      = K(end);
        
        % simulate state space
        Atransition        = eye(Nstates);
        hh                = 1;
        Shat              = K;
        VARprocess(hh,mm) = cc * Shat;
        checkdiff(1, cc * Shat); % should be unit impact
        SPFresponses(:,hh,mm) = CC * Shat;
        
        for hh = 2 : Hstar
            Atransition           = AA * Atransition;
            Shat                  = Atransition * Shat;
            VARprocess(hh,mm)     = cc * Shat;
            SPFresponses(:,hh,mm) = CC * Shat;
        end
        
    end
    
    
    %% plot VAR vs MDS process
    MDSmid       = median(MDSprocess,2);
    MDStail68    = prctile(MDSprocess, normcdf11, 2);
    VARmid       = median(VARprocess, 2);
    VARtail68    = prctile(VARprocess, normcdf11, 2);
    thisfig   = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hMDS = plotCI(MDSmid, MDStail68, 0:Hstar-1, []);
    hVAR    = plot(0:Hstar-1, VARmid, '-', 'color', colors4plots('blue'), 'LineWidth', 2);
    plot(0:Hstar-1, VARtail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
    plotOrigin
    xlim([0 Hstar-1])
    legend([hMDS,hVAR], 'MDS', 'VAR', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('MDSvsVARprocess-%s', modellabel), wrap)
    
    %% plot VAR process
    midstar   = median(VARkalmanStar);
    midstar68 = prctile(VARkalmanStar, normcdf11);
    mid       = median(VARprocess,2);
    tail68    = prctile(VARprocess, normcdf11, 2);
    thisfig   = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hMA   = plot(0:Hstar-1, mid, '-', 'color', colors4plots('blue'), 'LineWidth', 2);
    hMA68 = plot(0:Hstar-1, tail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
    hStar = yline(midstar, ':', 'color', colors4plots('orange'), 'LineWidth', 3);
    yline(midstar68, '-', 'color', colors4plots('orange'), 'LineWidth', 1);
    plotOrigin
    xlim([0 Hstar-1])
    wrapthisfigure(thisfig, sprintf('VARprocess-%s', modellabel), wrap)
    legend([hMA,hStar], 'IMA process', 'endpoint shift', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('VARprocess-%s-WITHLEGEND', modellabel), wrap)
    
    %% plot VAR process with SPF responses
    midstar   = median(VARkalmanStar);
    midstar68 = prctile(VARkalmanStar, normcdf11);
    mid       = median(VARprocess,2);
    tail68    = prctile(VARprocess, normcdf11, 2);
    midSPF    = median(SPFresponses,3);
    tailSPF   = prctile(SPFresponses, normcdf11, 3);
    thisfig   = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    hMA   = plot(0:Hstar-1, mid, '-', 'color', colors4plots('black'), 'LineWidth', 5);
    % hMA68 = plot(0:Hstar-1, tail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
    % hStar = yline(midstar, ':', 'color', colors4plots('orange'), 'LineWidth', 3);
    % yline(midstar68, '-', 'color', colors4plots('orange'), 'LineWidth', 1);
    Hshow = 9;
    hSPF = NaN(Hshow,1);
    for hh = 1 : Hshow
        hSPF(hh) = plot((0:Ngap-1) + (hh-1), midSPF(:,hh), '--', 'Linewidth', 2);
    end
    plotOrigin
    xlim([0 Hstar-1])
    legend(hSPF, arrayfun(@(i) sprintf('%d', i), 0:Hshow -1, 'uniformoutput', false));
    % legend([hMA,hStar], 'IMA process', 'endpoint shift', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('VARprocessSPF-%s', modellabel), wrap)
    
    
    % plot cumulative responses
    if ~strcmpi(datalabel, 'UNRATE')
        thisfig   = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        thisVARpath = cumsum(mid);
        hMA   = plot(0:Hstar-1, thisVARpath, '-', 'color', colors4plots('black'), 'LineWidth', 5);
        % hMA68 = plot(0:Hstar-1, tail68, '-', 'color', colors4plots('blue'), 'LineWidth', 1);
        % hStar = yline(midstar, ':', 'color', colors4plots('orange'), 'LineWidth', 3);
        % yline(midstar68, '-', 'color', colors4plots('orange'), 'LineWidth', 1);
        hSPF = NaN(Hshow,1);
        for hh = 1 : Hshow
            thisSPFpath  = cumsum(midSPF(:,hh));
            if hh > 1
                thisSPFpath = thisSPFpath + thisVARpath(hh-1);
            end
            hSPF(hh) = plot((0:Ngap-1) + (hh-1), thisSPFpath, '--', 'Linewidth', 2);
        end
        plotOrigin
        xlim([0 Hstar-1])
        legend(hSPF, arrayfun(@(i) sprintf('%d', i), 0:Hshow -1, 'uniformoutput', false), 'location', 'best');
        % legend([hMA,hStar], 'IMA process', 'endpoint shift', 'location', 'best')
        wrapthisfigure(thisfig, sprintf('cumulativeVARprocessSPF-%s', modellabel), wrap)
    end
    
    %% adjust datevector etc for plotting
    datesT     = dates(1:thisT);
    
    %% plot Ydraws over time
    Nstates = max(Nstates1,Nstates2);
    for n = 1 : Nstates
        these  = squeeze(Ydraws1(min(n,Nstates1),:,:));
        mid1   = median(these, 2);
        tails1 = prctile(these, normcdf([-1 1]) * 100, 2);
        
        these  = squeeze(Ydraws2(min(n,Nstates2),:,:));
        mid2   = median(these, 2);
        tails2 = prctile(these, normcdf([-1 1]) * 100, 2);
        
        thisfig = figure;
        hold on
        h0 = plot(datesT, mid1, '-', 'color', color1, 'LineWidth', 2);
        plot(datesT, tails1, '--', 'color', color1, 'LineWidth', 1);
        h2 = plot(datesT, mid2, '-', 'color', color2, 'LineWidth', 2);
        plot(datesT, tails2, '--', 'color', color2, 'LineWidth', 1);
        if any(ismember(ndxFixedHorizon, n))
            thisZ = Zdata(1:thisT, n);
            thisZ(Znanny(1:thisT,n)) = NaN;
            hSPF = plot(datesT, thisZ, 'rx', 'LineWidth', .5);
        else
            hSPF = [];
        end
        
        xtickdates(datesT)
        if n == 1
            title(sprintf('data (lag 1)'))
        elseif ~isempty(hSPF)
            title(sprintf('forecast h=%d', n-2))
            legend([h0 h2 hSPF], modelname1, modelname2, 'SPF', 'location', 'best')
        else
            title(sprintf('forecast h=%d', n-2))
            legend([h0 h2], modelname1, modelname2, 'location', 'best')
        end
        wrapthisfigure(thisfig, sprintf('Ydraws%d-%s', n, modellabel), wrap)
    end
    clear Nstates
    
    
    %% plot SV
    Nsv = 2;
    % SVtdraws1 = SVtdraws1 ./ SVtdraws1(:,1,:);
    % SVtdraws2 = SVtdraws2 ./ SVtdraws2(:,1,:);
    % SVdraws1  = SVdraws1 ./ SVdraws1(:,1,:);
    % SVdraws2  = SVdraws2 ./ SVdraws2(:,1,:);
    
    midt1   = median(SVtdraws1,3);
    % tailst1 = prctile(SVtdraws1, normcdf([-1 1]) * 100, 3);
    midt2   = median(SVtdraws2,3);
    % tailst2 = prctile(SVtdraws2, normcdf([-1 1]) * 100, 3);
    
    mid1   = median(SVdraws1,3);
    % tails1 = prctile(SVdraws1, normcdf([-1 1]) * 100, 3);
    mid2   = median(SVdraws2,3);
    % tails2 = prctile(SVdraws2, normcdf([-1 1]) * 100, 3);
    
    for ii = 1 : Nsv
        thisfig = figure;
        subplot(2,1,1)
        hold on
        set(gca, 'fontsize', 18)
        h1 = plot(datesT,midt1(ii,:), '-', 'color', color1, 'LineWidth', 2);
        % plot(datesT,squeeze(tailst1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
        h2 = plot(datesT,midt2(ii,:), '-', 'color', color2, 'LineWidth', 2);
        % plot(datesT,squeeze(tailst2(ii,:,:)), '-', 'color', color2, 'LineWidth', 1)
        xtickdates(datesT)
        legend([h1, h2], modelname1, modelname2, 'location', 'best')
        title('SV-t')
        subplot(2,1,2)
        hold on
        set(gca, 'fontsize', 18)
        h1 = plot(datesT,mid1(ii,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
        % plot(datesT,squeeze(tails1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
        h2 = plot(datesT,mid2(ii,:), '-', 'color', color2, 'LineWidth', 2); %#ok<NASGU>
        % plot(datesT,squeeze(tails2(ii,:,:)), '-', 'color', color2, 'LineWidth', 1)
        xtickdates(datesT)
        % legend([h1, h2], modelname1, modelname2, 'location', 'best')
        title('SV')
        wrapthisfigure(thisfig, sprintf('SVblock%d-%s', ii, modellabel), wrap)
    end
    
    %% compare histograms of tdofdraws
    for ii = 1 : Nsv
        thisfig = figure;
        hold on
        set(gca, 'fontsize', 18)
        histogram(tdofdraws1(ii,:), 'Normalization', 'pdf', 'FaceAlpha', .5);
        histogram(tdofdraws2(ii,:), 'Normalization', 'pdf', 'FaceAlpha', .5);
        legend(modelname1, modelname2, 'location', 'best')
        title(sprintf('tdof Block %d', ii))
        wrapthisfigure(thisfig, sprintf('tdof-block%d-%s', ii, modellabel), wrap)
    end
    
    %% report VAR for ETA
    Gprior = chol(diag(G0var2), 'lower') * randn(Ng2, MCMCdraws);
    Gprior = reshape(Gprior, Ngap2, Ngap2, MCMCdraws);
    Gprior = permute(Gprior, [2 1 3]);
    GpriorlambdaDraws = NaN(MCMCdraws,1);
    for mm = 1 : MCMCdraws
        GpriorlambdaDraws(mm) = max(abs(eig(Gprior(:,:,mm))));
    end
    thisfig = figure;
    [~, h2, hprior] = plotpriorposteriordraws(Glambda2, GpriorlambdaDraws, 0 : .01 : 1);
    h2.Color=color2;
    legend([h2 hprior], 'VAR', 'prior')
    wrapthisfigure(thisfig, sprintf('Geigenvalues2-%s', modellabel), wrap)
    
    
    
    %% Compute CG regressions model2
    Ngap      = Ngap2;
    Nstates   = Ngap + 1;
    H         = Nstates - 2 - 1;
    Hpools    = {0:3};
    ndxSV11   = 1:5;
    ndxSV22   = ndxSV11(end)+1:Ngap2; % OK, let's hard code this. Could use cell arrays for easier extension to N blocks, but leave that for later
    
    CGbeta2   = NaN(H,MCMCdraws);
    CGpooled2 = NaN(length(Hpools),MCMCdraws);
    maxlogESV = 10; % max vol of about 150
    for m = 1 : MCMCdraws
        
        
        % thisrho           = hrhodraws2(:,m);
        % logESV            = diag(hvcvdraws2(:,:,m)) ./ (1 - thisrho.^2) ...
        %     + log(tdofdraws2(:,m)) - log((tdofdraws2(:,m) - 2));
        % logESV(logESV > maxlogESV) = maxlogESV;
        % ESV               = exp(.5 .* logESV);
        
        ESV              = median(SVtdraws2(:,:,mm),2);
        
        sigTrend                  = mustarvardraws2(m);
        
        sqrtSigmaTilde            = sqrtsigmaDraws2(:,:,m);
        sqrtSigmaTilde(:,ndxSV11) = sqrtSigmaTilde(:,ndxSV11) .* ESV(1);
        sqrtSigmaTilde(:,ndxSV22) = sqrtSigmaTilde(:,ndxSV22) .* ESV(2);
        SigmaTilde                = sqrtSigmaTilde * sqrtSigmaTilde';
        
        G           = Gdraws2(:,:,m);
        VCV         = dlyapdoubling(G, SigmaTilde);
        
        CGpooled2(:,m) = CGregressionpooledTC(G,VCV,sigTrend,Hpools);
        % CG regression
        CGbeta2(:,m) = CGregressionTC(G,VCV,sigTrend,H);
        
    end
    clear Ngap Nstates
    
    %% plot densities of CG slopes
    CGbeta2 = transpose(CGbeta2);
    thisfig = figure;
    for h = 1 : H
        subplot(ceil(H/2),2,h)
        [~, h2] = plotpriorposteriordraws(CGbeta2(:,h));
        h2.Color=color2;
        xline(0, '-', 'LineWidth',1)
        title(sprintf('h = %d', h))
    end
    sgtitle(sprintf('%s \n CG slopes', datalabel))
    wrapthisfigure(thisfig, sprintf('CGslopes%d-%s', H, modellabel), wrap)
    
    %% plot pooled CG slopes
    for ii = 1 : length(Hpools)
        thisfig = figure;
        [~,h2] = plotpriorposteriordraws(CGpooled2(ii,:));
        h2.Color=color2;
        xline(0, '-', 'LineWidth',1)
        thisPool = Hpools{ii};
        title(sprintf('%s -- Pooled CG slope (h=%d:%d)', datalabel, thisPool(1), thisPool(end)))
        wrapthisfigure(thisfig, sprintf('CGpooled%d-%s', ii, modellabel), wrap)
    end
    
    %% plot medians of noise levels
    Zbarlabel = Zlabel(ndxZbar);
    ndxQ4 = quarter(datesT) == 4;
    mid = median(NOISEdraws1,3);
    mid2 = median(NOISEdraws2,3);
    for nn = 1 : Nbar
        firstObs = find(~isnan(mid(nn,:)), 1);
        if ~isempty(firstObs)
            thisfig = figure;
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2);
            h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
            
            h2 = plot(datesT, mid2(nn,:), ':', 'color', color2, 'LineWidth', 2);
            h2q4 = plot(datesT(ndxQ4), mid2(nn,ndxQ4), 'd', 'color', color2, 'LineWidth', 2);
            xtickdates(datesT(firstObs:end))
            plotOrigin;
            legend([h1, h1q4, h2, h2q4], modelname1, 'Q4', modelname2, 'Q4', 'location', 'best')
            title(sprintf('NOISE levels %s', Zbarlabel{nn}))
            wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s', nn, modellabel), wrap)
        end
    end
    
    %% plot medians of noise variances
    mid = median(NOISEvarstate1,3);
    mid2 = median(NOISEvarstate2,3);
    for nn = 1 : Nbar
        firstObs = find(~isnan(mid(nn,:)), 1);
        if ~isempty(firstObs)
            thisfig = figure;
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2);
            h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
            
            h2 = plot(datesT, mid2(nn,:), ':', 'color', color2, 'LineWidth', 2);
            h2q4 = plot(datesT(ndxQ4), mid2(nn,ndxQ4), 'd', 'color', color2, 'LineWidth', 2);
            xtickdates(datesT(firstObs:end))
            plotOrigin;
            legend([h1, h1q4, h2, h2q4], modelname1, 'Q4', modelname2, 'Q4', 'location', 'best')
            title(sprintf('NOISE variances %s', Zbarlabel{nn}))
            wrapthisfigure(thisfig, sprintf('NOISEvariances-y%d-%s', nn, modellabel), wrap)
        end
    end
    
    %% plot predictive density of Y
    colorALT = colors4plots("lightblue");
    
    Zdata(Znanny) = NaN;
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    
    mid1   = mean(YdensityDraws1,2);
    tails1 = prctile(YdensityDraws1, [5 normcdf([-1 1]) * 100 95], 2);
    h1 = plotCI(mid1, tails1, 0:Nhorizons1-1,[], '-', 'color', 'black', 'LineWidth', 4);

    h2 = plot(0:Nhorizons2-1,mean(YdensityDraws2,2), '-', 'color', colorALT, 'LineWidth', 4);
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, normcdf([-1 1]) * 100, 2), '-', 'color', colorALT, 'LineWidth', 2)
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, [5 95], 2), '-.', 'color', colorALT, 'LineWidth', 2)
    
    % redo mid1 as dashes
    plot(0:Nhorizons1-1, mid1, 'k--', 'linewidth', 3);

    xmax = max(Nhorizons1,Nhorizons2) - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s', modellabel), wrap)
    hl = legend([h1 h2], modelname1, modelname2,'location', 'best');
    wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGEND', modellabel), wrap)
    title(sprintf('%s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGENDTITLE', modellabel), wrap)
    YtermLIM = ylim; % for Yterm plot
    
    %% plot term structures Ydraws (at end of sample)
    
    thisH = max(Nstates1,Nstates2)-2;
    
    Zdata(Znanny) = NaN;
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % add plot of term structure estimates
    h1 = plot(0:Nstates1-2,mean(Ydraws1(2:Nstates1,end,:),3), '-', 'color', color1, 'Linewidth', 2);
    theseTtails = squeeze(prctile(Ydraws1(2:Nstates1,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates1-2,theseTtails, '--', 'color', color1, 'Linewidth', 2);
    
    h2 = plot(0:Nstates2-2,mean(Ydraws2(2:Nstates2,end,:),3), '-', 'color', color2, 'Linewidth', 2);
    theseTtails = squeeze(prctile(Ydraws2(2:Nstates2,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates2-2,theseTtails, '--', 'color', color2, 'Linewidth', 2);
    % add SPF plots
    SPFcolor = 'black';
    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end
    ndxB            = 2 + 4 + 1;
    horizonsB       = (4 - datesQ(thisT)) + MAweights;
    ndxC            = 2 + 4 + 2;
    horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
    ndxD            = 2 + 4 + 3;
    horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
    hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', SPFcolor, 'linewidth', 2); %#ok<NASGU>
    if ~Znanny(thisT,ndxB)
        hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', SPFcolor, 'linewidth', 2);
    end
    if ndxC <= Nz && ~Znanny(thisT,ndxC)
        hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', SPFcolor, 'linewidth', 2);
    end
    if ndxD <= Nz && ~Znanny(thisT,ndxD)
        hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', SPFcolor, 'linewidth', 2);
    end
    
    xticks(0 : 2 : thisH)
    xlim([0 thisH])
    legend([h1 h2], modelname1, modelname2,'location', 'best');
    title(sprintf('Term structure: %s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
    % ylim(YtermLIM)
    wrapthisfigure(thisfig, sprintf('Yterm-%s-WITHLEGENDTITLE', modellabel), wrap)
    
    
    %% for model1: plot predictive density together with term structure of Ydraws (at end of sample)
    colorTerm = colors4plots('black');
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % plot density
    hDense = plot(0:Nhorizons1-1,mean(YdensityDraws1,2), '-', 'color', color1, 'LineWidth', 4);
    plot(0:Nhorizons1-1,prctile(YdensityDraws1, normcdf([-1 1]) * 100, 2), '-', 'color', color1, 'LineWidth', 2)
    plot(0:Nhorizons1-1,prctile(YdensityDraws1, [5 95], 2), '-.', 'color', color1, 'LineWidth', 2)
    % plot term structure estimates
    hTerm = plot(0:Nstates1-2,mean(Ydraws1(2:Nstates1,end,:),3), '--', 'color', colorTerm, 'Linewidth', 4);
    theseTtails = squeeze(prctile(Ydraws1(2:Nstates1,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates1-2,theseTtails, '--', 'color', colorTerm, 'Linewidth', 1);
    % axis settings
    xmax = max(Nhorizons1,Nhorizons2) - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s', modelname1, modellabel), wrap)
    legend([hDense hTerm], 'Predictive Density', 'SPF-consistent term structure','location', 'best');
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGEND', modelname1, modellabel), wrap)
    title(sprintf('%s: %s per %s', modelname1, datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGENDTITLE', modelname1, modellabel), wrap)
    
    %% for model1: plot predictive density plus SPF and SEP (at end of sample)
    colorSPF = colors4plots('orange');
    colorSEP = colors4plots('green');
    colorCGM = 'black';
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % plot density
    hDense = plot(0:Nhorizons1-1,mean(YdensityDraws1,2), '-', 'color', colorCGM, 'LineWidth', 4);
    plot(0:Nhorizons1-1,prctile(YdensityDraws1, normcdf([-1 1]) * 100, 2), '-', 'color', colorCGM, 'LineWidth', 2)
    plot(0:Nhorizons1-1,prctile(YdensityDraws1, [5 95], 2), '-.', 'color', colorCGM, 'LineWidth', 2)
    % add SPF
    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end
    ndxB            = 2 + 4 + 1;
    horizonsB       = (4 - datesQ(thisT)) + MAweights;
    ndxC            = 2 + 4 + 2;
    horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
    ndxD            = 2 + 4 + 3;
    horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
    hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', colorSPF, 'linewidth', 2);
    if ~Znanny(thisT,ndxB)
        hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
    if ndxC <= Nz && ~Znanny(thisT,ndxC)
        hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
    if ndxD <= Nz && ~Znanny(thisT,ndxD)
        hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
 
    ndxCalendarYears = (4 - datesQ(thisT)) + 1:4:12;
    xline(ndxCalendarYears, ':')
    % axis settings
    xmax = max(Nhorizons1,Nhorizons2) - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s', modelname1, modellabel), wrap)
    legend([hDense hSPF], 'Predictive Density', 'SPF', 'location', 'best');
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGEND', modelname1, modellabel), wrap)
    title(sprintf('%s: %s per %s', modelname1, datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGENDTITLE', modelname1, modellabel), wrap)
    if ~isempty(SEPfcst)
        % add SEP
        thisSEPndx = SEPfcstDates == dates(thisT);
        horizonsY{1} = (4 - datesQ(thisT)) + (1:4) - 4;
        horizonsY{2} = (4 - datesQ(thisT)) + (1:4);
        horizonsY{3} = (4 - datesQ(thisT)) + (1:4) + 4;
        horizonsY{4} = (4 - datesQ(thisT)) + (1:4) + 8;
        for yy = 1 : 4
            hSEP = plot(horizonsY{yy}, repmat(SEPfcst(thisSEPndx,yy), 1, 4), '-d', 'color', colorSEP, 'linewidth', 2);
        end
        legend([hDense hSPF hSEP(1)], 'Predictive Density', 'SPF', 'SEP', 'location', 'best');
        title(sprintf('%s: %s per %s', modelname1, datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('YpredictivedensitySPFSEP-%s-%s-WITHLEGENDTITLE', modelname1, modellabel), wrap)
    end

    %% for model2: plot predictive density together with term structure of Ydraws (at end of sample)
    colorTerm = colors4plots('black');
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % plot density
    hDense = plot(0:Nhorizons2-1,mean(YdensityDraws2,2), '-', 'color', color2, 'LineWidth', 4);
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, normcdf([-1 1]) * 100, 2), '-', 'color', color2, 'LineWidth', 2)
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, [5 95], 2), '-.', 'color', color2, 'LineWidth', 2)
    % plot term structure estimates
    hTerm = plot(0:Nstates2-2,mean(Ydraws2(2:Nstates2,end,:),3), '--', 'color', colorTerm, 'Linewidth', 4);
    theseTtails = squeeze(prctile(Ydraws2(2:Nstates2,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates2-2,theseTtails, '--', 'color', colorTerm, 'Linewidth', 1);
    % axis settings
    xmax = max(Nhorizons1,Nhorizons2) - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s', modelname2, modellabel), wrap)
    legend([hDense hTerm], 'Predictive Density', 'SPF-consistent term structure','location', 'best');
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGEND', modelname2, modellabel), wrap)
    title(sprintf('%s: %s per %s', modelname2, datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGENDTITLE', modelname2, modellabel), wrap)
    
    %% for model2: plot predictive density plus SPF and SEP (at end of sample)
    colorSPF = colors4plots('orange');
    colorSEP = colors4plots('green');
    colorCGM = 'black';
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % plot density
    hDense = plot(0:Nhorizons2-1,mean(YdensityDraws2,2), '-', 'color', colorCGM, 'LineWidth', 4);
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, normcdf([-1 1]) * 100, 2), '-', 'color', colorCGM, 'LineWidth', 2)
    plot(0:Nhorizons2-1,prctile(YdensityDraws2, [5 95], 2), '-.', 'color', colorCGM, 'LineWidth', 2)
    % add SPF
    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end
    ndxB            = 2 + 4 + 1;
    horizonsB       = (4 - datesQ(thisT)) + MAweights;
    ndxC            = 2 + 4 + 2;
    horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
    ndxD            = 2 + 4 + 3;
    horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
    SPFcolor =  [0 0 0];
    hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', colorSPF, 'linewidth', 2);
    if ~Znanny(thisT,ndxB)
        hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
    if ndxC <= Nz && ~Znanny(thisT,ndxC)
        hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
    if ndxD <= Nz && ~Znanny(thisT,ndxD)
        hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', colorSPF, 'linewidth', 2);
    end
    
    ndxCalendarYears = (4 - datesQ(thisT)) + 1:4:12;
    xline(ndxCalendarYears, ':')
    % axis settings
    xmax = max(Nhorizons1,Nhorizons2) - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s', modelname2, modellabel), wrap)
    legend([hDense hSPF], 'Predictive Density', 'SPF', 'location', 'best');
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGEND', modelname2, modellabel), wrap)
    title(sprintf('%s: %s per %s', modelname2, datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('YpredictivedensitySPF-%s-%s-WITHLEGENDTITLE', modelname2, modellabel), wrap)
    
    if ~isempty(SEPfcst)
        % add SEP
        thisSEPndx = SEPfcstDates == dates(thisT);
        horizonsY{1} = (4 - datesQ(thisT)) + (1:4) - 4;
        horizonsY{2} = (4 - datesQ(thisT)) + (1:4);
        horizonsY{3} = (4 - datesQ(thisT)) + (1:4) + 4;
        horizonsY{4} = (4 - datesQ(thisT)) + (1:4) + 8;
        for yy = 1 : 4
            hSEP = plot(horizonsY{yy}, repmat(SEPfcst(thisSEPndx,yy), 1, 4), '-d', 'color', colorSEP, 'linewidth', 2);
        end
        legend([hDense hSPF hSEP(1)], 'Predictive Density', 'SPF', 'SEP', 'location', 'best');
        title(sprintf('%s: %s per %s', modelname2, datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('YpredictivedensitySPFSEP-%s-%s-WITHLEGENDTITLE', modelname2, modellabel), wrap)
    end
    
    %% wrap up
    dockAllFigures
    finishwrap
    
    if ~isempty(wrap) && doStore
        close all
        filename = fullfile(wrap.dir, modellabel);
        save(filename, '-v7.3');
    end
    
end

%% finish / clean up
finishscript
