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

NGAP       = 'BOP';
MDSLABEL   = 'MDS';

thinSTEPS = 1;

doStore = true;

%% general parameters
samStart = []; % 52 is 1981Q3

%% batch parameters
% SED-PARAMETERS-HERE

%% some settings

fontsize = 18;

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

SEPrmseData     = importdata('../rawdataKensington/SEPTable2Data.xlsx');

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
    
    RTdata   = matfile(fullfile(datadir, sprintf('kensington%sdataRT', upper(datalabel))));
    
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
    [Ngap, NGAP] = switchNGAP(NGAP, Ny, datalabel, thisDate);
    
    modellabel = sprintf('%s-%strendHScycleSV2blockNoiseHS-Ngap%s-samStart%s-%s-thin%d', ...
        datalabel, MDSLABEL, NGAP, startDateLabel, thisDateLabel, thinSTEPS);
    
    
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
    
    %% read SEP uncertainty data
    switch datalabel
        case 'RGDP'
            sheetname = 'RealGDP';
        case 'UNRATE'
            sheetname = 'UnemploymentRate';
        case 'CPI'
            sheetname = 'TotalConsumerPrices';
            warning('Matching CPI uncertainty with PCE from SEP')
        otherwise
            warning('no SEP for <<%s>>', datalabel)
            sheetname = [];
    end
    
    if isempty(sheetname)
        SEPbandsDates = [];
        SEPrmse       = [];
    else
        if ispc
            SEPbandsDates = datenum(SEPrmseData.textdata.(sheetname)(2:end,1));
            SEPrmse       = SEPrmseData.data.(sheetname);
        else
            SEPbandsDates = x2mdate(SEPrmseData.data.(sheetname)(:,1)); %#ok<XMDATE>
            SEPrmse       = SEPrmseData.data.(sheetname)(:,2:end);
        end
    end
    
    %% call MCMC sampler Ngap1
    modelname = MDSLABEL; % sprintf('Ngap%d', Ngap1);
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
                Ngap, ...
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
    
    
    Nhorizons = size(YdensityDraws1,1);
    Nstates   = Ngap + 1;
    
    %% collect PSRF and INFF for MCMC1
    % Ydraws1
    thesedraws          = Ydraws1;
    % NaN out the known values
    for ii = 1 : 6
        ndx = Znanny(1:thisT,ii);
        thesedraws(ii,~ndx,:) = NaN;
    end
    thesedraws = transpose(reshape(thesedraws, Nstates * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: Ydraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: Ydraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    % SV
    Nsv = 2;
    thesedraws = transpose(reshape(SVdraws1, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    % SVt
    thesedraws = transpose(reshape(SVtdraws1, Nsv * thisT, MCMCdraws));
    plotPSRF(thesedraws, sprintf('%s: SVtdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: SVtdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    % tdofdraws
    thesedraws = transpose(tdofdraws1);
    plotPSRF(thesedraws, sprintf('%s: tdofdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    plotINFF(thesedraws, sprintf('%s: tdofdraws1', MDSLABEL), strcat(modelname, '-', modellabel), wrap);
    
    %% MDS: construct univariate process for y(t)
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
        
        % collect Sigma
        % thisrho          = hrhodraws1(:,mm);
        % logESV           = diag(hvcvdraws1(:,:,mm)) ./ (1 - thisrho.^2) ...
        %     + log(tdofdraws1(:,mm)) - log((tdofdraws1(:,mm) - 2));
        % logESV(logESV > maxlogESV) = maxlogESV;
        % ESV              = exp(.5 .* logESV);
        
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
    
    
    %% adjust datevector etc for plotting
    datesT     = dates(1:thisT);
    
    %% plot Ydraws over time
    for n = 1 : Nstates
        these  = squeeze(Ydraws1(min(n,Nstates),:,:));
        mid1   = median(these, 2);
        tails1 = prctile(these, normcdf([-1 1]) * 100, 2);
        
        thisfig = figure;
        hold on
        h0 = plot(datesT, mid1, '-', 'color', color1, 'LineWidth', 2);
        plot(datesT, tails1, '--', 'color', color1, 'LineWidth', 1);
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
        else
            title(sprintf('forecast h=%d', n-2))
        end
        wrapthisfigure(thisfig, sprintf('Ydraws%d-%s', n, modellabel), wrap)
    end
    
    
    %% plot SV
    Nsv = 2;
    
    midt1   = median(SVtdraws1,3);
    % tailst1 = prctile(SVtdraws1, normcdf([-1 1]) * 100, 3);
    
    mid1   = median(SVdraws1,3);
    % tails1 = prctile(SVdraws1, normcdf([-1 1]) * 100, 3);
    
    for ii = 1 : Nsv
        thisfig = figure;
        subplot(2,1,1)
        hold on
        set(gca, 'fontsize', 18)
        h1 = plot(datesT,midt1(ii,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
        % plot(datesT,squeeze(tailst1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
        xtickdates(datesT)
        % legend([h1, h2], modelname, modelname2, 'location', 'best')
        title('SV-t')
        subplot(2,1,2)
        hold on
        set(gca, 'fontsize', 18)
        h1 = plot(datesT,mid1(ii,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
        % plot(datesT,squeeze(tails1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
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
        % histogram(tdofdraws2(ii,:), 'Normalization', 'pdf', 'FaceAlpha', .5);
        % legend(modelname, modelname2, 'location', 'best')
        title(sprintf('tdof Block %d', ii))
        wrapthisfigure(thisfig, sprintf('tdof-block%d-%s', ii, modellabel), wrap)
    end
    
    
    %% plot medians of noise levels
    Zbarlabel = Zlabel(ndxZbar);
    ndxQ4 = quarter(datesT) == 4;
    mid = median(NOISEdraws1,3);
    % mid2 = median(NOISEdraws2,3);
    for nn = 1 : Nbar
        firstObs = find(~isnan(mid(nn,:)), 1);
        if ~isempty(firstObs)
            thisfig = figure;
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
            h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
            
            % h2 = plot(datesT, mid2(nn,:), ':', 'color', color2, 'LineWidth', 2);
            % h2q4 = plot(datesT(ndxQ4), mid2(nn,ndxQ4), 'd', 'color', color2, 'LineWidth', 2);
            xtickdates(datesT(firstObs:end))
            plotOrigin;
            % legend([h1, h1q4, h2, h2q4], modelname, 'Q4', modelname2, 'Q4', 'location', 'best')
            title(sprintf('NOISE levels %s', Zbarlabel{nn}))
            wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s', nn, modellabel), wrap)
        end
    end
    
    %% plot medians of noise variances
    mid = median(NOISEvarstate1,3);
    % mid2 = median(NOISEvarstate2,3);
    for nn = 1 : Nbar
        firstObs = find(~isnan(mid(nn,:)), 1);
        if ~isempty(firstObs)
            thisfig = figure;
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
            h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
            
            % h2 = plot(datesT, mid2(nn,:), ':', 'color', color2, 'LineWidth', 2);
            % h2q4 = plot(datesT(ndxQ4), mid2(nn,ndxQ4), 'd', 'color', color2, 'LineWidth', 2);
            xtickdates(datesT(firstObs:end))
            plotOrigin;
            % legend([h1, h1q4, h2, h2q4], modelname, 'Q4', modelname2, 'Q4', 'location', 'best')
            title(sprintf('NOISE variances %s', Zbarlabel{nn}))
            wrapthisfigure(thisfig, sprintf('NOISEvariances-y%d-%s', nn, modellabel), wrap)
        end
    end
    
    %% plot predictive density of Y
    Zdata(Znanny) = NaN;
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % note: using the mean (to match observed)
    % plot density
    h1 = plot(0:Nhorizons-1,mean(YdensityDraws1,2), '-', 'color', color1, 'LineWidth', 4); %#ok<NASGU>
    plot(0:Nhorizons-1,prctile(YdensityDraws1, normcdf([-1 1]) * 100, 2), '-', 'color', color1, 'LineWidth', 2)
    plot(0:Nhorizons-1,prctile(YdensityDraws1, [5 95], 2), '-.', 'color', color1, 'LineWidth', 2)
    
    xmax = Nhorizons - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    title(sprintf('%s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGENDTITLE', modellabel), wrap)
    YtermLIM = ylim; % for Yterm plot
    
    %% plot term structures Ydraws (at end of sample)
    
    thisH = Nstates - 2;
    
    Zdata(Znanny) = NaN;
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % add plot of term structure estimates
    h1 = plot(0:Nstates-2,mean(Ydraws1(2:Nstates,end,:),3), '-', 'color', color1, 'Linewidth', 2);
    theseTtails = squeeze(prctile(Ydraws1(2:Nstates,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates-2,theseTtails, '--', 'color', color1, 'Linewidth', 2);
    
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
    hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', SPFcolor, 'linewidth', 2); 
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
    % legend([h1 h2], modelname, modelname2,'location', 'best');
    title(sprintf('Term structure: %s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
    % ylim(YtermLIM)
    wrapthisfigure(thisfig, sprintf('Yterm-%s-WITHLEGENDTITLE', modellabel), wrap)
    
    
    %% for model1: plot predictive density together with term structure of Ydraws (at end of sample)
    colorTerm = colors4plots('black');
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % plot density
    hCGM = plot(0:Nhorizons-1,mean(YdensityDraws1,2), '-', 'color', color1, 'LineWidth', 4);
    plot(0:Nhorizons-1,prctile(YdensityDraws1, normcdf([-1 1]) * 100, 2), '-', 'color', color1, 'LineWidth', 2)
    plot(0:Nhorizons-1,prctile(YdensityDraws1, [5 95], 2), '-.', 'color', color1, 'LineWidth', 2)
    % plot term structure estimates
    hTerm = plot(0:Nstates-2,mean(Ydraws1(2:Nstates,end,:),3), '--', 'color', colorTerm, 'Linewidth', 4);
    theseTtails = squeeze(prctile(Ydraws1(2:Nstates,end,:), normcdf([-1 1]) * 100, 3));
    plot(0:Nstates-2,theseTtails, '--', 'color', colorTerm, 'Linewidth', 1);
    % axis settings
    xmax = Nhorizons - 1;
    xticks(0 : 2 : xmax)
    xlim([0 xmax])
    legend([hCGM hTerm], 'Predictive Density', 'SPF-consistent term structure','location', 'best');
    title(sprintf('%s: %s per %s', modelname, datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGENDTITLE', modelname, modellabel), wrap)
    
    
    %% construct Q4 draws for YAVG
    % YAVGdraws are constructued for each quarterly horizon, Q4 will be picked below
    switch datalabel
        case {'RGDP','CPI','PGDP'}
            Kavg = 4;
            YAVGdraws = constructYAVGthisT(YdensityDraws1, thisT, dates, datalabel, RTdata, Kavg);
        case 'UNRATE'
            YAVGdraws = YdensityDraws1;
        otherwise
            error('datalabel <<%s>> not recognized', datalabel)
    end

    %% for model1: ANNUAL FAN CHART with SPF and SEP (at end of sample)
    if ~isempty(SEPfcst)
        ndxQ4horizons = (4 - datesQ(thisT)) + (1:4:Nhorizons);
        ndxQ4horizons = ndxQ4horizons(ndxQ4horizons<=Nhorizons);
        Nannual       = length(ndxQ4horizons);
        
        tickYears     = year(dates(thisT)) + (1:Nannual) - 1;

        colorSEP = colors4plots('orange');
        colorCGM = colors4plots('blue');
        
        % collect estimates
        mid    = mean(YAVGdraws,2);
        tail68 = prctile(YAVGdraws, normcdf11, 2);
        % tail90 = prctile(YAVGdraws, [5 95], 2);

        % plot figure
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density
        hCGM = plot(1:Nannual,mid(ndxQ4horizons), '-d', 'color', colorCGM, 'LineWidth', 2, 'MarkerSize', 10);
        plot(1:Nannual,mid(ndxQ4horizons), '-', 'color', colorCGM, 'LineWidth', 4);
        plot(1:Nannual, tail68(ndxQ4horizons,:), '-', 'color', colorCGM, 'LineWidth', 2)
        % plot(1:Nannual, tail90(ndxQ4horizons,:), '-.', 'color', colorCGM, 'LineWidth', 2)
        % add SEP
        thisSEPndx = SEPfcstDates == dates(thisT);
        thisSEPfcst = SEPfcst(thisSEPndx,:);
        hSEP       = plot(1:4, thisSEPfcst, '--s', 'color', colorSEP, 'linewidth', 2, 'MarkerSize', 10);
        plot(1:4, thisSEPfcst, '--', 'color', colorSEP, 'linewidth', 4);
        thisSEPndx = SEPbandsDates == dates(thisT);
        SEPupper   = SEPfcst(thisSEPndx,:) + SEPrmse(thisSEPndx,:);
        SEPlower   = SEPfcst(thisSEPndx,:) - SEPrmse(thisSEPndx,:);
        plot(1:4, SEPupper, '--', 'color', colorSEP, 'linewidth', 2);
        plot(1:4, SEPlower, '--', 'color', colorSEP, 'linewidth', 2);
        
        % title(sprintf('%s: %s per %s', modelname, datalabel, datestr(datesT(end), 'yyyyqq')))
        % axis settings
        xticks(1 : 4)
        xticklabels(tickYears)
        xlim([1 4])
        ylim0
        wrapthisfigure(thisfig, sprintf('fanchartSEP-%s', modellabel), wrap)
        legend([hCGM hSEP(1)], 'CGM', 'SEP', 'location', 'best');
        wrapthisfigure(thisfig, sprintf('fanchartSEP-%s-WITHLEGEND', modellabel), wrap)
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
