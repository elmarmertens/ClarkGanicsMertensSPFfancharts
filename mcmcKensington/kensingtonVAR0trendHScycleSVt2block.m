%% VAR-TRENDCYCLE-SV model with horseshoe shocks to trend and multivariate SV-t-scale, and 2block SV
% replaces MDS assumption by VAR(1) in forecast updates


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

Nstreams    = max(1,getparpoolsize);
rndStreams  = initRandStreams(Nstreams, [], 0);

% beginning and end of eval window (set to empty if to use max)
samStart     = [];
samEnd       = [];

doSamStartSPF   = false; % if TRUE: set samStart to first SPF forecast origin
doTBILLcensored = true;
ELB             = .25; % note: not effective if ELB set below observed FFR

quicky     = false; % if TRUE: very short MCMC chains, no looping across variables,

doStore    = true; % to store result files in directory specified by "localstore" m-file
doStoreXXL = false; % to store extended result files (with MCMC draws)
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'};

Tstart     = 60;

quantileP    = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
Nquantiles   = length(quantileP);
percentiles  = 1 : 99;
Npercentiles = length(percentiles);

CGpools      = {0:3};
CGmaxhorizon = 5;

Kavg = 4;

NGAP = 'BOP';


fontsize = 18;
datadir = fullfile('..', 'matdataKensington');

doSPFquarterlyOnly = false; % note: Noise model not needed if true
doY1Q4             = false; % do *not* include Y1Q4 (default)

%% placeholder for batch parameters
% SED-PARAMETERS-HERE

%% process parameters: MCMCdraws etc
if ~exist('MCMCdraws', 'var') || isempty(MCMCdraws)
    if quicky
        MCMCdraws    = 1e2;
    else
        MCMCdraws    = 3e3;
    end
end

if ~exist('Nfedraws', 'var') || isempty(Nfedraws)
    if quicky
        Nfedraws     = 10;
    else
        Nfedraws     = 100;
    end
end

%% loop over variables
tic

for d =  1 : length(DATALABELS)
    
    close all
    datalabel = DATALABELS{d};
    
    NfcstDraws  = Nfedraws * MCMCdraws;
    NfcstDraws2 = 2 * NfcstDraws; % for density draws with antithetic sampling
    
    %% load data
    matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
    load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', 'YBARlabel', 'Nbar', 'YBARfuture', 'YCURRENTfuture', 'YAVGfuture', ...
        'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czhat', 'Czbar')
    
    RTdata   = matfile(fullfile(datadir, sprintf('kensington%sdataRT', upper(datalabel))));
    
    matfilename = fullfile(datadir,sprintf('kensington%sdataHISTOGRAMS.mat',DATALABELS{d}));
    if exist(matfilename, 'file')
        PROBdata = matfile(matfilename);
        if dates ~= PROBdata.dates
            error('mismatch in dates with SPF histogram data')
        end
        SPFhistograms = PROBdata.SPFhistograms;
    else
        PROBdata = [];
        SPFhistograms = [];
    end
    
    %% pick trend-cycle parameters
    Ngap      = switchNGAP(NGAP,Ny,datalabel);
    Nstates   = Ngap + 1;
    Nhorizons = max(Ny, Nstates - 1); % Ny includes lagged data, Nhorizons counts forecast steps
    Nhorizons = 17;
    
    %% select start of sample
    if isempty(samStart)
        if doSamStartSPF
            samStart   = find(sum(~Znanny,2) > 1,1); % first availability of SPF
        else
            samStart   = find(any(~Znanny,2),1); % for CPI/TBILL: leaves early sample with only Yrealized ...
        end
    end
    
    %% prepare wrapper
    
    if doTBILLcensored & strcmpi(datalabel, 'TBILL')
        modellabel = sprintf('%scensored-VAR0trendHScycleSVt2block-Ngap%s-samStart%s', datalabel, NGAP, datestr(dates(samStart), 'yyyyQQ'));
        doCENSOR  = true;
    else
        modellabel = sprintf('%s-VAR0trendHScycleSVt2block-Ngap%s-samStart%s', datalabel, NGAP, datestr(dates(samStart), 'yyyyQQ'));
        doCENSOR   = false;
    end
    
    if ~doSPFquarterlyOnly

        if doY1Q4
            modellabel = replace(modellabel, '-Ngap', '-y1q4-Ngap');
        else
            modellabel = replace(modellabel, '-Ngap', '-EXy1q4-Ngap');

            % drop Y1Q4 from Zdata
            ndxQ4                = quarter(dates) == 4;
            ndxY1                = 7; % y(-1), h=0:4 means six elements in Z before Y1
            Zdata(ndxQ4, ndxY1)  = NaN;
            Znanny(ndxQ4, ndxY1) = true;
            Cz(ndxY1,:,ndxQ4)    = 0; % should be superfluous, since mcmcsamplers align Cz with Zanny anyway

        end

    else
        
        modellabel = replace(modellabel, '-Ngap', '-SPFquarterlyOnly-Ngap');
        
        %% drop annual predictions from Zdata
        Zdata  = Zdata(:,1:Nz-Nbar);
        Znanny = Znanny(:,1:Nz-Nbar);
        Cz     = Cz(1:Nz-Nbar,:,:);
        Nz     = Nz - Nbar;
    end
    
    wrap = [];
    titlename = modellabel;
    if ~quicky
        initwrap
    end
    
    fprintf('Processing %s ... \n', modellabel)
    
    %% prepare data
    dates      = dates(samStart:end); datesQ = datesQ(samStart:end);
    T          = length(dates);
    
    Zdata      = Zdata(samStart:end,:);
    Znanny     = Znanny(samStart:end,:);
    Cz         = Cz(:,:,samStart:end);
    
    % adjustment to align outcomes with estimation sample
    Yfuture        = Yfuture(samStart:end, :);
    YBARfuture     = YBARfuture(samStart:end, :);
    YCURRENTfuture = YCURRENTfuture(samStart:end, :);
    YAVGfuture     = YAVGfuture(samStart:end, :);
    if ~isempty(SPFhistograms)
        SPFhistograms = SPFhistograms(samStart:end);
    end

    Nsv = 2; % two SV blocks

    %% plot Z elements
    %     if ~quicky
    %         for n = 1 : Nz
    %             thisZ = Zdata(:,n);
    %             thisZ(Znanny(:,n)) = NaN;
    %             thisfig = figure;
    %             hold on
    %             plot(dates,thisZ)
    %             title(sprintf('%s', Zlabel{n}))
    %             xtickdates(dates)
    %             wrapthisfigure(thisfig,sprintf('Z%d%s', n, modellabel), wrap)
    %             close(thisfig)
    %         end
    %     end
    
    %% allocate memory for QRT runs
    Ydraws          = NaN(Nstates,MCMCdraws,T);
    ETAdraws        = NaN(Nstates,MCMCdraws,T);
    SVdraws         = NaN(Nsv,MCMCdraws,T);
    SVtdraws        = NaN(Nsv,MCMCdraws,T);
    tDOFdraws       = NaN(Nsv,MCMCdraws,T);
    hvcvdraws       = NaN(Nsv,Nsv,MCMCdraws,T);
    hrhodraws       = NaN(Nsv,MCMCdraws,T);

    mustarvoldraws  = NaN(MCMCdraws,T);
    
    Gdraws          = cell(T,1); % cell to handle minStateT
    sqrtSIGMAdraws  = cell(T,1);
    G0var           = cell(T,1); % cell to handle minStateT
    
    % in keeping with the CMM template, allocate memory for T jumpoffs,
    % even though we have only T-Tstart+1 jumpoffs
    
    fcstZhat          = NaN(T,Nz);
    fcstZhaterror     = NaN(T,Nz);
    fcstZmvlogscore   = NaN(T,1);
    
    YdensityDraws     = NaN(Nhorizons,NfcstDraws2,T); % permute later
    fcstYhat          = NaN(T,Nhorizons);  % predictive mean
    fcstYhatRB        = NaN(T,Nhorizons); % predictive mean (linear RB)
    fcstYmedian       = NaN(T,Nhorizons); % predictive median
    fcstYhaterror     = NaN(T,Nhorizons);
    fcstYhatRBerror   = NaN(T,Nhorizons);
    fcstYmederror     = NaN(T,Nhorizons);
    fcstYlogscore     = NaN(T,Nhorizons);
    fcstYcrps         = NaN(T,Nhorizons);
    fcstYdrps         = NaN(T,Nhorizons);
    fcstYquantiles    = NaN(T,Nhorizons,Nquantiles);
    fcstYpercentiles  = NaN(T,Nhorizons,Npercentiles);
    
    fcstYvol        = NaN(T,Nhorizons);
    fcstYskew       = NaN(T,Nhorizons);
    
    YBARdensityDraws     = NaN(Nbar,NfcstDraws2,T); % permute later
    fcstYBARmedian       = NaN(T,Nbar);
    fcstYBARhat          = NaN(T,Nbar);
    fcstYBARvol          = NaN(T,Nbar);
    fcstYBARskew         = NaN(T,Nbar);
    fcstYBARquantiles    = NaN(T,Nbar,Nquantiles);
    fcstYBARpercentiles  = NaN(T,Nbar,Npercentiles);
    fcstYBARhaterror     = NaN(T,Nbar);
    fcstYBARmederror     = NaN(T,Nbar);
    fcstYBARlogscore     = NaN(T,Nbar);
    fcstYBARcrps         = NaN(T,Nbar);
    fcstYBARdrps         = NaN(T,Nbar);
    
    YCURRENTdensityDraws     = NaN(1,NfcstDraws2, T); % permute later
    fcstYCURRENThat          = NaN(T,1);
    fcstYCURRENTmedian       = NaN(T,1);
    fcstYCURRENTvol          = NaN(T,1);
    fcstYCURRENTskew         = NaN(T,1);
    fcstYCURRENTquantiles    = NaN(T,1,Nquantiles);
    fcstYCURRENTpercentiles  = NaN(T,1,Npercentiles);
    fcstYCURRENThaterror     = NaN(T,1);
    fcstYCURRENTmederror     = NaN(T,1);
    fcstYCURRENTlogscore     = NaN(T,1);
    fcstYCURRENTcrps         = NaN(T,1);
    fcstYCURRENTdrps         = NaN(T,1);
    
    NgapT                = NaN(T,1);
    
    CGPOOLEDdraws = NaN(length(CGpools),MCMCdraws,T);
    CGSLOPESdraws = NaN(CGmaxhorizon,MCMCdraws,T);
    
    [YFINALdraws, SVFINALdraws, ...
        ETAFINALdraws] = deal(cell(T,1)); % trick to collect only thisT==T output inside parfor
    
    logtwopi = log(2 * pi);
    
    %% if TBILL: apply censoring to future values
    if doCENSOR
        ndxELB = Yfuture < ELB;
        Yfuture(ndxELB) = ELB;
        % TODO: apply censoring to future values of YBAR, YCURRENT and YAVG
    end
    
    %% parfor loop over QRT jumpoffs
    warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
    parfor thisT = Tstart : T % @parfor
        warning('off','backtrace')
        TID = parid;
        
        thisNgap     = switchNGAP(NGAP,Ny,datalabel,dates(thisT));
        thisNstates  = thisNgap + 1;
        NgapT(thisT) = thisNgap;
        
        %% prepare prior for G (or leave empty)
        thisG0var  = NaN(thisNgap*thisNgap,1); % diagonal elements
        Gtheta = [0.2 0.5].^2;
        ndx = 0;
        for i=1:thisNgap
            for j=1:thisNgap
                ndx  = ndx + 1;
                if (j==i)
                    thisG0var(ndx)=Gtheta(1);
                else
                    thisG0var(ndx)=Gtheta(1)*Gtheta(2);
                end
            end
        end
        thisG0var = diag(thisG0var);
        
        %% call mcmc sampler, and catch any numerical errors
        OK = false;
        while ~OK
            try
                [mcmc_Ydensitydraws, mcmc_Y, mcmc_YhatRB, mcmc_ETA, ~, ~, ...
                    mcmc_G, mcmc_Glambda, ~, mcmc_VARCONST, ... 
                    mcmc_SV, mcmc_SVt, mcmc_hVCV, mcmc_RHO, mcmc_tDOF, ~, ...
                    mcmc_sqrtSIGMA, ~, mcmc_mustarSIG, ...
                    Ztp1mean, Ztp1meanerror, Zmvlogscore] = ...
	                     mcmcsamplerVAR0trendHScycleSVt2block(Zdata, Znanny, Cz, thisNgap, ...
                    thisG0var, ...
                    thisT, MCMCdraws, Nfedraws, rndStreams{TID}, false);
                OK = true;
            catch mcmcME
                warning('MCMC ERROR: at thisT=%d, message: %s', thisT, getReport(mcmcME, 'basic'))
            end
        end
        
        %% apply censoring (if desired)
        if doCENSOR
            ndxELB = mcmc_Ydensitydraws < ELB;
            mcmc_Ydensitydraws(ndxELB) = ELB;
            
            % notes:
            % - censoring applied only to TBILL
            % - calendar year TBILL draws are averages of quarterly draws; thus no censoring for calendar year draws needed
        end
        %% collect mcmc objects
        % need to go via tmp variable for parfor compatibility
        
        % ETA
        tmp                      = NaN(Nstates,MCMCdraws);
        tmp(1:thisNstates,:)     = mcmc_ETA(:,end,:);
        ETAdraws(:,:,thisT)      = tmp;
        
        % Y
        tmp                      = NaN(Nstates,MCMCdraws);
        tmp(1:thisNstates,:)     = mcmc_Y(:,end,:);
        if thisNstates < Nstates
            tmp(thisNstates+1:Nstates,:) = repmat(tmp(thisNstates,:), Nstates - thisNstates,1);
        end
        Ydraws(:,:,thisT)          = tmp;


        SVdraws(:,:,thisT)         = mcmc_SV(:,end,:);
        SVtdraws(:,:,thisT)        = mcmc_SVt(:,end,:);
        tDOFdraws(:,:,thisT)       = mcmc_tDOF;

        hvcvdraws(:,:,:,thisT)     = mcmc_hVCV;
        hrhodraws(:,:,thisT)       = mcmc_RHO;
        mustarvoldraws(:,thisT)    = sqrt(mcmc_mustarSIG);
        
        G0var{thisT}               = thisG0var;
        
        Gdraws{thisT}              = mcmc_G;
        sqrtSIGMAdraws{thisT}      = mcmc_sqrtSIGMA;
        
        % predictive densities for quarterly horizons
        YdensityDraws(:,:,thisT)    = mcmc_Ydensitydraws;
        fcstYquantiles(thisT,:,:)   = prctile(mcmc_Ydensitydraws, quantileP, 2);
        fcstYpercentiles(thisT,:,:) = prctile(mcmc_Ydensitydraws, percentiles, 2);
        fcstYhat(thisT,:)           = mean(mcmc_Ydensitydraws, 2);
        fcstYhatRB(thisT,:)       = mean(mcmc_YhatRB, 2);
        fcstYmedian(thisT,:)      = median(mcmc_Ydensitydraws, 2);
        fcstYhaterror(thisT,:)    = fcstYhat(thisT,:)    - Yfuture(thisT,:);
        fcstYhatRBerror(thisT,:)  = fcstYhatRB(thisT,:)  - Yfuture(thisT,:);
        fcstYmederror(thisT,:)    = fcstYmedian(thisT,:) - Yfuture(thisT,:);
        
        fcstYvol(thisT,:)         = std(mcmc_Ydensitydraws, 0, 2);
        fcstYskew(thisT,:)        = skewBowleyKelly(mcmc_Ydensitydraws, 2);
        
        
        thesedraws = transpose(sort(mcmc_Ydensitydraws, 2)); % just flip dimensions for more efficient computations in next few steps
        % univariate log scores via kerneldensities
        for hh = 1 : Nhorizons
            fcstYlogscore(thisT,hh) = log(ksdensity(thesedraws(:,hh), Yfuture(thisT,hh)));
        end
        % CRPS
        for hh = 1 : Nhorizons
            fcstYcrps(thisT,hh) = crpsDraws(Yfuture(thisT,hh), thesedraws(:,hh), true);
        end
        
        %% DRPS
        if ~isempty(SPFhistograms)
            theseBinEdges = SPFhistograms(thisT).binEdges;
            if ~isempty(theseBinEdges)
                fcstYdrps(thisT,:) = drpsDraws(Yfuture(thisT,:), thesedraws, theseBinEdges, [], true);
            end
        end
        
        %% CG slopes
        ndxSV11     = 1:5;
        ndxSV22     = ndxSV11(end)+1:thisNgap; % OK, let's hard code this. Could use cell arrays for easier extension to N blocks, but leave that for later

        maxlogESV      = 10; % max vol of about 150

        for mm = 1 : MCMCdraws
            % compute VCV of shocks
  
            % thisrho                    = mcmc_RHO(:,mm);
            % logESV                     = diag(mcmc_hVCV(:,:,mm)) ./ (1 - thisrho.^2) + log(mcmc_tDOF(:,mm)) - log((mcmc_tDOF(:,mm) - 2));
            % logESV(logESV > maxlogESV) = maxlogESV;
            % ESV                        = exp(.5 .* logESV);

            ESV = median(mcmc_SVt(:,:,mm),2); % median of SVt

            sigTrend                  = mcmc_mustarSIG(mm);
            sqrtSigmaTilde            = mcmc_sqrtSIGMA(:,:,mm);
            sqrtSigmaTilde(:,ndxSV11) = sqrtSigmaTilde(:,ndxSV11) .* ESV(1);
            sqrtSigmaTilde(:,ndxSV22) = sqrtSigmaTilde(:,ndxSV22) .* ESV(2);
            SigmaTilde                = sqrtSigmaTilde * sqrtSigmaTilde';

            % VCV of ETA
            thisGtilde   = mcmc_G(:,:,mm);
            thisVCVtilde = dlyapdoubling(thisGtilde, SigmaTilde); 
            % pooled CG regression
            CGPOOLEDdraws(:,mm,thisT) = CGregressionpooledTC(thisGtilde,thisVCVtilde,sigTrend,CGpools);
            CGSLOPESdraws(:,mm,thisT) = CGregressionTC(thisGtilde,thisVCVtilde,sigTrend,CGmaxhorizon);
        end % mm
        
        %% Z predictions and density score
        if thisT < T
            fcstZmvlogscore(thisT)  = Zmvlogscore;
            thisZ                   = NaN(Nz,1); % helper variable to make parfor work
            ndx                     = ~Znanny(thisT+1,:);
            thisZ(ndx)              = Ztp1mean;
            fcstZhat(thisT,:)       = thisZ;
            thisZ(ndx)              = Ztp1meanerror;
            fcstZhaterror(thisT,:)  = thisZ;
        end
        
        %% predictive densities for calendar years
        
        % prepare RTdata
        thisDate      = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI'}
                % use latest vintage, not real-time data
                RT_vintage_idx = length(RTdata.vindates); % pick last vintage
            case {'RGDP', 'PGDP'}
                RT_vintage_idx = find(datenum(RTdata.vindates) == thisDate,1,'first'); % vintage, time of SPF round
            otherwise
                RT_vintage_idx = []; %#ok<NASGU> % to avoid parfor warning
                error('datalabel <<%s>> not supported', datalabel)
        end
        % last observation available at time of SPF round
        RT_obsdate_idx = find(datenum(RTdata.obsdates)==datenum(dateshift(datetime(thisDateLabel,'InputFormat','yyyyQQQ'),'start','quarter',-1)),1,'first');
        % vector of vintage data
        RT_vec = RTdata.data(1:RT_obsdate_idx,RT_vintage_idx);
        
        % transform draws
        [theseYBARdraws, theseYCURRENTdraws] = trfQ2A(doNIPA, thisT, datesQ, RT_obsdate_idx, ...
            mcmc_Ydensitydraws, RT_vec, NfcstDraws2, Nbar);
        
        YBARdensityDraws(:,:,thisT)    = theseYBARdraws;
        fcstYBARhat(thisT,:)           = mean(theseYBARdraws,2);
        fcstYBARmedian(thisT,:)        = median(theseYBARdraws,2);
        fcstYBARvol(thisT,:)           = std(theseYBARdraws,0,2);
        fcstYBARskew(thisT,:)          = skewBowleyKelly(theseYBARdraws, 2);
        fcstYBARquantiles(thisT,:,:)   = prctile(theseYBARdraws,quantileP,2);
        fcstYBARpercentiles(thisT,:,:) = prctile(theseYBARdraws,percentiles,2);
        
        fcstYBARhaterror(thisT,:)      = fcstYBARhat(thisT,:) - YBARfuture(thisT,:);
        fcstYBARmederror(thisT,:)      = fcstYBARmedian(thisT,:) - YBARfuture(thisT,:);
        
        theseYBARdraws = transpose(sort(theseYBARdraws, 2)); % just flip dimensions for more efficient computations in next few steps
        % univariate log scores via kerneldensities
        for hh = 1 : Nbar
            fcstYBARlogscore(thisT,hh) = log(ksdensity(theseYBARdraws(:,hh), YBARfuture(thisT,hh)));
        end
        % CRPS
        for hh = 1 : Nbar
            fcstYBARcrps(thisT,hh) = crpsDraws(YBARfuture(thisT,hh), theseYBARdraws(:,hh), true);
        end
        % DRPS
        if ~isempty(SPFhistograms)
            theseBinEdges = SPFhistograms(thisT).binEdges;
            if ~isempty(theseBinEdges)
                fcstYBARdrps(thisT,:) = drpsDraws(YBARfuture(thisT,:), theseYBARdraws, theseBinEdges, [], true);
            end
        end
        
        %% YCURRENT
        if ~any(isnan(theseYCURRENTdraws), 'all') % can occur due to data gap in lagged data used to construct YCURRENT
            YCURRENTdensityDraws(:,:,thisT)    = theseYCURRENTdraws;
            fcstYCURRENThat(thisT,:)           = mean(theseYCURRENTdraws,2);
            fcstYCURRENTmedian(thisT,:)        = median(theseYCURRENTdraws,2);
            fcstYCURRENTvol(thisT,:)           = std(theseYCURRENTdraws,0,2);
            fcstYCURRENTskew(thisT,:)          = skewBowleyKelly(theseYCURRENTdraws,2);
            fcstYCURRENTquantiles(thisT,:,:)   = prctile(theseYCURRENTdraws,quantileP,2);
            fcstYCURRENTpercentiles(thisT,:,:) = prctile(theseYCURRENTdraws,percentiles,2);
            
            fcstYCURRENThaterror(thisT,:)      = fcstYCURRENThat(thisT,:) - YCURRENTfuture(thisT,:);
            fcstYCURRENTmederror(thisT,:)      = fcstYCURRENTmedian(thisT,:) - YCURRENTfuture(thisT,:);
            
            theseYCURRENTdraws = transpose(sort(theseYCURRENTdraws, 2)); % just flip dimensions for more efficient computations in next few steps
            % univariate log scores via kerneldensities
            fcstYCURRENTlogscore(thisT, 1) = log(ksdensity(theseYCURRENTdraws(:,1), YCURRENTfuture(thisT,1)));
            % CRPS
            fcstYCURRENTcrps(thisT,1) = crpsDraws(YCURRENTfuture(thisT,1), theseYCURRENTdraws(:,1), true);
            % DRPS
            if ~isempty(SPFhistograms)
                theseBinEdges = SPFhistograms(thisT).binEdges;
                if ~isempty(theseBinEdges)
                    fcstYCURRENTdrps(thisT,:) = drpsDraws(YCURRENTfuture(thisT,:), theseYCURRENTdraws, theseBinEdges, [], true);
                end
            end
        end
        
        %% collect output for thisT=T
        if  thisT == T
            YFINALdraws{thisT}    = mcmc_Y;
            ETAFINALdraws{thisT}  = mcmc_ETA;
            SVFINALdraws{thisT}   = mcmc_SV;
            
        else
            YFINALdraws{thisT}      = [];
            ETAFINALdraws{thisT}    = [];
            SVFINALdraws{thisT}     = [];
        end
        
        fprintf('%s - QRT: done with t=%d (Ngap=%d)\n', modellabel, thisT, thisNgap)
        
        warning('on','backtrace')
        
    end % parfor
    
    YFINALdraws     = YFINALdraws{T};
    ETAFINALdraws   = ETAFINALdraws{T};
    SVFINALdraws    = SVFINALdraws{T};
    
    % permute draws
    Ydraws               = permute(Ydraws, [1 3 2]);
    ETAdraws             = permute(ETAdraws, [1 3 2]);
    SVdraws              = permute(SVdraws, [1 3 2]);
    SVtdraws             = permute(SVtdraws, [1 3 2]);
    tDOFdraws            = permute(tDOFdraws, [1 3 2]);
    
    
    YdensityDraws        = permute(YdensityDraws, [3 1 2]);
    YBARdensityDraws     = permute(YBARdensityDraws, [3 1 2]);
    YCURRENTdensityDraws = permute(YCURRENTdensityDraws, [3 1 2]);
    
    
    %% YRECESS (for quarterly 0:4)
    if strcmpi(datalabel, 'RGDP')
        thesedraws  = (exp((YdensityDraws(:,1:5,:) ./ 100)) - 1) .* 100;
        % note: undoing the logs should be an abundance of caution; at the
        % zero-threshold both log- and simple changes should be identical
        fcstYrecess = mean(thesedraws < 0, 3) * 100;
        
        checkdiff(fcstYrecess, mean(YdensityDraws(:,1:5,:) < 0, 3) * 100);
    end
    
    %% YAVG: draws
    % note: the construction of YAVG could also be moved inside the parfor
    % loop (kept here for legacy reasons)
    
    % collect sequence of historical realizations to complete the current year
    ylags = NaN(T,1,Kavg-1);
    parfor thisT = Tstart : T
        thisDate      = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI'}
                % use latest vintage, not real-time data
                RT_vintage_idx = length(RTdata.vindates); % pick last vintage
            case {'RGDP', 'PGDP'}
                RT_vintage_idx = find(datenum(RTdata.vindates) == thisDate,1,'first'); % vintage, time of SPF round
            otherwise
                RT_vintage_idx = []; %#ok<NASGU> % to avoid parfor warning
                error('datalabel <<%s>> not supported', datalabel)
        end
        % last observation available at time of SPF round
        RT_obsdate_idx = find(datenum(RTdata.obsdates)==datenum(dateshift(datetime(thisDateLabel,'InputFormat','yyyyQQQ'),'start','quarter',-1)),1,'first');
        % vector of vintage data
        RT_vec = RTdata.data(1:RT_obsdate_idx,RT_vintage_idx);
        
        % collect three lags, note that RT_obsdate_idx is one lag relative to thisT
        ylags(thisT,1,:) = RT_vec(end-(Kavg-1)+1:end); % store lag three first
    end
    
    % collect draws
    ydraws    = permute(YdensityDraws, [1 3 2]);
    YAVGdensityDraws = NaN(size(ydraws));
    % h < Kavg
    for h = 1 : (Kavg - 1)
        YAVGdensityDraws(:,:,h) = (sum(ydraws(:,:,1:h), 3) + sum(ylags(:,1,end-(Kavg-1)+h:end), 3)) / Kavg;
    end
    % h >= Kavg
    for h = Kavg : size(YAVGdensityDraws, 3)
        YAVGdensityDraws(:,:,h) = sum(ydraws(:,:,h+(-(Kavg-1):0)), 3) / Kavg;
    end
    YAVGdensityDraws = permute(YAVGdensityDraws, [1 3 2]);
    
    
    %% clear helper objects
    clear ydraws
    clear ylags
    
    %% patch YAVG: compute stats
    
    fcstYAVGhat          = mean(YAVGdensityDraws, 3);
    fcstYAVGmedian       = median(YAVGdensityDraws, 3);% predictive median
    fcstYAVGvol          = std(YAVGdensityDraws, 0, 3);
    fcstYAVGskew         = skewBowleyKelly(YAVGdensityDraws, 3);
    fcstYAVGquantiles    = prctile(YAVGdensityDraws, quantileP, 3);
    fcstYAVGpercentiles  = prctile(YAVGdensityDraws, percentiles, 3);
    
    fcstYAVGhatRB  = NaN(size(fcstYhatRB));
    for h = Kavg : size(fcstYAVGhatRB, 2)
        fcstYAVGhatRB(:,h) = sum(fcstYhatRB(:,h+(-(Kavg-1):0)), 2) / Kavg;
    end
    
    fcstYAVGhaterror   = fcstYAVGhat    - YAVGfuture;
    fcstYAVGhatRBerror = fcstYAVGhatRB  - YAVGfuture;
    fcstYAVGmederror   = fcstYAVGmedian - YAVGfuture;
    
    fcstYAVGlogscore   = NaN(T,Nhorizons);
    fcstYAVGcrps       = NaN(T,Nhorizons);
    for thisT = Tstart : T
        thesedraws = permute(YAVGdensityDraws(thisT,:,:), [3 2 1]);
        for hh = Kavg : Nhorizons
            fcstYAVGcrps(thisT,hh)      = crpsDraws(YAVGfuture(thisT,hh), thesedraws(:,hh));
            fcstYAVGlogscore(thisT,hh)  = log(ksdensity(thesedraws(:,hh), YAVGfuture(thisT,hh)));
        end
    end
    clear thesedraws
    
    %% patch PITvalues
    % note: comparisons with NaN return 0 (or false); hence some corrections
    fcstYpits                        = mean(YdensityDraws < Yfuture, 3);
    fcstYpits(isnan(Yfuture))        = NaN;
    fcstYpits(1:Tstart-1,:)          = NaN;
    
    fcstYBARpits                     = mean(YBARdensityDraws < YBARfuture, 3);
    fcstYBARpits(isnan(YBARfuture))  = NaN;
    fcstYBARpits(1:Tstart-1,:)       = NaN;
    
    fcstYCURRENTpits                        = mean(YCURRENTdensityDraws  < YCURRENTfuture, 3);
    fcstYCURRENTpits(isnan(YCURRENTfuture)) = NaN;
    fcstYCURRENTpits(any(isnan(YCURRENTdensityDraws), 3)) = NaN; % to handle missing obs in construct YCURRENTdraws
    fcstYCURRENTpits(1:Tstart-1,:)          = NaN;
    
    fcstYAVGpits                     = mean(YAVGdensityDraws < YAVGfuture, 3);
    fcstYAVGpits(isnan(YAVGfuture))  = NaN;
    fcstYAVGpits(1:Tstart-1,:)       = NaN;
    
    %% G posterior
    
    % posterior
    [GlambdaDraws, GpriorlambdaDraws] = deal(NaN(MCMCdraws, T));
    parfor tt = Tstart : T
        
        Ng = NgapT(tt) * NgapT(tt);
        
        % prior
        Gprior = chol(G0var{tt}, 'lower') * randn(Ng, MCMCdraws);
        Gprior = reshape(Gprior, NgapT(tt), NgapT(tt), MCMCdraws);
        Gprior = permute(Gprior, [2 1 3]);
        
        for mm = 1 : MCMCdraws
            GpriorlambdaDraws(mm,tt) = max(abs(eig(Gprior(:,:,mm))));
        end
        
        
        theseGdraws = Gdraws{tt};
        for mm = 1 : MCMCdraws
            GlambdaDraws(mm,tt) = max(abs(eig(theseGdraws(:,:,mm))));
        end
    end
    
    GpriorlambdaMedian = median(GpriorlambdaDraws,1);
    GpriorlambdaTails  = prctile(GpriorlambdaDraws, [25 75], 1);
    
    GlambdaMedian = median(GlambdaDraws,1);
    GlambdaTails  = prctile(GlambdaDraws, [25 75], 1);
    
    thisfig = figure;
    hold on
    hprior = plot(dates, GpriorlambdaMedian, 'k--', 'linewidth', 2);
    plot(dates, GpriorlambdaTails, 'k--')
    
    hpost = plot(dates, GlambdaMedian, 'b-', 'linewidth', 2);
    plot(dates, GlambdaTails, 'b-')
    
    
    ylim([min([ylim, 0]) max(ylim)])
    xtickdates(dates([Tstart T]))
    
    legend([hprior, hpost], 'prior', 'posterior');
    wrapthisfigure(thisfig, sprintf('Geigenvalues-%s', modellabel), wrap)
    
    clear GlambdaDraws GpriorlambdaDraws
    
    %% plot vol and skew of predictive density
    thisfig = figure;
    plot(dates, fcstYvol)
    set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
    legend(Ylabel(2:end), 'location', 'best')
    xtickdates(dates(Tstart:end))
    wrapthisfigure(thisfig, sprintf('predictiveVol-%s', modellabel), wrap)
    
    thisfig = figure;
    plot(dates, fcstYskew)
    set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
    legend(Ylabel(2:end), 'location', 'best')
    xtickdates(dates(Tstart:end))
    wrapthisfigure(thisfig, sprintf('predictiveSkew-%s', modellabel), wrap)
    
    thisfig = figure;
    plot(dates, fcstYBARvol)
    set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
    legend(YBARlabel, 'location', 'best')
    xtickdates(dates(Tstart:end))
    wrapthisfigure(thisfig, sprintf('predictiveYBARVol-%s', modellabel), wrap)
    
    thisfig = figure;
    plot(dates, fcstYBARskew)
    set(gca,'linestyleorder',{'-','-.',':','--','*','+'})
    legend(YBARlabel, 'location', 'best')
    xtickdates(dates(Tstart:end))
    wrapthisfigure(thisfig, sprintf('predictiveYBARSkew-%s', modellabel), wrap)
    
    %% compare CRPS vs DRPS
    for nn = 1 : Nhorizons
        thisfig = figure;
        hanni = plot(dates, [fcstYdrps(:,nn) fcstYcrps(:,nn)]);
        xtickdates(dates(Tstart:T))
        legend(hanni, 'DRPS', 'CRPS')
        title(sprintf('h=%d', nn-1))
        wrapthisfigure(thisfig, sprintf('DRPS-CRPS-h%d-%s', nn, modellabel), wrap)
    end
    
    %% check censoring of calendar year draws
    if doCENSOR
        if any(YBARdensityDraws < ELB, 'all')
            warning('YBARdensityDraws contains draws below ELB')
        end
        if any(YCURRENTdensityDraws < ELB, 'all')
            warning('YCURRENTdensityDraws contains draws below ELB')
        end
    end
    %% store results
    diary off
    if  doStore && ~quicky
        if doStoreXXL
            save(fullfile(localstore, modellabel), ...
                '-v7.3');
        end
        clear foo thesedraws
        close all
        if doStoreXXL
            matname = sprintf('CGMmcmc-%s-Ndraws%d', modellabel, MCMCdraws);
            save(fullfile(localstore, matname), '-v7.3');
        end
        clear Y*density*
        clearvars *draws *Draws -except SV*draws tDOFdraws CG*draws Ydraws YFINALdraws MCMCdraws  Nfedraws
        matname = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, MCMCdraws);
        save(fullfile(localstore, matname), '-v7.3');
    end
    finishwrap % if latex compilation fails, edit this script to turn the compilation off
    if quicky
        dockAllFigures
        toc
        return
    else
        close all
    end
end

toc

%% finish / clean up
finishscript
dockAllFigures
