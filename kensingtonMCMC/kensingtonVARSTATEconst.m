%% VAR-STATE-CONST model
% replaces MDS assumption by VAR(1) in forecast updates

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*PFBNS>
%#ok<*UNRCH>

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
samStart = [];
samEnd   = [];
% samEnd   = find(dates < datenum(2020,1,1), 1, 'last'); % stop before COVID

quicky     = false; % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
doStore    = true; % to store result files in directory specified by "localstore" m-file
doStoreXXL = false; % to store extended result files (with MCMC draws)
DATALABELS = {'RGDP', 'PGDP', 'UNRATE', 'TBILL', 'CPI'};

Tstart     = 60;

quantileP    = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
Nquantiles   = length(quantileP);
percentiles  = 1 : 99;
Npercentiles = length(percentiles);

Kavg = 4;

VARpriorshrinkage = 5; % will be divided by 1000

fontsize = 18;

datadir = fullfile('..', 'kensingtonDataMatfiles');

tic

%% loop over variables

for d =  1 : length(DATALABELS)

    close all
    datalabel = DATALABELS{d};

    if quicky
        MCMCdraws    = 1e2;
        Nfedraws     = 10;
    else
        MCMCdraws    = 3e3;
        Nfedraws     = 100;
    end

    NfcstDraws = Nfedraws * MCMCdraws;


    %% load data
    matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
    load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', 'YBARlabel', 'Nbar', 'YBARfuture', 'YCURRENTfuture', ...
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

    %% prepare wrapper
    modellabel = strcat(datalabel, sprintf('VARprior%dSTATEconst', VARpriorshrinkage));


    wrap = [];
    titlename = modellabel;
    if ~quicky
        initwrap
    end

    fprintf('Processing %s ... \n', modellabel)

    %% record location of missing values (NaN)

    samStart   = find(any(~Znanny,2),1); % for CPI/TBILL: leaves early sample with only Yrealized ...
    % samStart   = find(sum(~Znanny,2) > 1,1); % first availability of SPF

    dates      = dates(samStart:end);
    T          = length(dates);
    Nhorizons  = Ny - 1;

    Zdata      = Zdata(samStart:end,:);
    Znanny     = Znanny(samStart:end,:);
    Cz         = Cz(:,:,samStart:end);

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

    %% prepare prior for G  
    G0var  = NaN(Ny*Ny,1); % diagonal elements
    Gtheta    = [0 0];
    Gtheta(1) = VARpriorshrinkage / 1000;
    Gtheta(2) = 0.5;
    ndx = 0;
    for i=1:Ny
        for j=1:Ny
            ndx  = ndx + 1;
            if (i==j)
                G0var(ndx)=Gtheta(1);
            else
                G0var(ndx)=Gtheta(1)*Gtheta(2);
            end
        end
    end
    G0var   = diag(G0var);

    %% allocate memory for QRT runs
    Ydraws          = NaN(Ny,T,MCMCdraws);
    ETAdraws        = NaN(Ny,T,MCMCdraws);
    SVdraws         = NaN(T,MCMCdraws);

    Gdraws          = NaN(Ny,Ny,MCMCdraws,T);

    % in keeping with the CMM template, allocate memory for T jumpoffs,
    % even though we have only T-Tstart+1 jumpoffs
    YdensityDraws     = NaN(T,Nhorizons,NfcstDraws);
    fcstYhat          = NaN(T,Nhorizons);  % predictive mean
    fcstYhatRB        = NaN(T,Nhorizons); % predictive mean (linear RB)
    fcstYmedian       = NaN(T,Nhorizons); % predictive median
    fcstYhaterror     = NaN(T,Nhorizons);
    fcstYhatRBerror   = NaN(T,Nhorizons);
    fcstYmederror     = NaN(T,Nhorizons);
    fcstYlogscore     = NaN(T,Nhorizons);
    fcstYlogscoreRB   = NaN(T,Nhorizons);
    fcstYcrps         = NaN(T,Nhorizons);
    fcstYdrps         = NaN(T,Nhorizons);
    fcstYquantiles    = NaN(T,Nhorizons,Nquantiles);
    fcstYpercentiles  = NaN(T,Nhorizons,Npercentiles);

    fcstYvol        = NaN(T,Nhorizons);
    fcstYskew       = NaN(T,Nhorizons);

    YBARdensityDraws     = NaN(T,Nbar,NfcstDraws);
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

    YCURRENTdensityDraws     = NaN(T,1,NfcstDraws);
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


    fcstYmvlogscoreDraws = NaN(T,MCMCdraws); % one-step ahead only NOTE: MCMCdraws sufficient for CONST model
    fcstYmvlogscore      = NaN(T,1); % one-step ahead only


    [YFINALdraws, ETAFINALdraws] = deal(cell(T,1)); % trick to collect only thisT==T output inside parfor

    logtwopi = log(2 * pi);

    %% parfor loop over QRT jumpoffs
    parfor thisT = Tstart : T % @parfor

        TID = parid;

        %% call mcmc sampler, and catch any numerical errors
        OK = false;
        while ~OK
            try
                [mcmc_Ydensitydraws, mcmc_Y, mcmc_YhatRB, mcmc_ETA, ...
                    mcmc_G, ... % Gdraws
                    mcmc_sqrtSIGMA, ...
                    YstateTp1mean, YstateTp1var, yTpHmean, yTpHvar ] = ...
                    mcmcsamplerVARSTATEconst(Zdata, Znanny, Cz, ...
                    G0var, ...
                    thisT, MCMCdraws, Nfedraws, rndStreams{TID}, false);
                OK = true;
            catch mcmcME
                warning('mcmc aborted at thisT=%d, message: %s', thisT, mcmcME.getReport)
            end
        end

        %% collect mcmc objects
        ETAdraws(:,thisT,:)              = mcmc_ETA(:,end,:);
        Ydraws(:,thisT,:)                = mcmc_Y(:,end,:);
        Gdraws(:,:,:,thisT)              = mcmc_G;

        % predictive densities for quarterly horizons
        YdensityDraws(thisT,:,:)  = mcmc_Ydensitydraws;
        fcstYquantiles(thisT,:,:)   = prctile(mcmc_Ydensitydraws, quantileP, 2);
        fcstYpercentiles(thisT,:,:) = prctile(mcmc_Ydensitydraws, percentiles, 2);
        fcstYhat(thisT,:)         = mean(mcmc_Ydensitydraws, 2);
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

        %% RB logscore
        for hh = 1 : Nhorizons
            thesevariances            = yTpHvar(:,hh);
            thesemeans                = yTpHmean(:,hh);
            thesescores               = -.5 * (logtwopi + log(thesevariances) + ...
                (Yfuture(thisT,hh) - thesemeans).^2 ./ thesevariances);
            maxscore                  = max(thesescores);
            fcstYlogscoreRB(thisT,hh) = log(mean(exp(thesescores - maxscore))) + maxscore;
        end

        % density score
        if thisT < T

            thatT  = thisT + 1;
            thatC  = Cz(~Znanny(thatT,:),:,thatT);
            thatZ  = Zdata(thatT,~Znanny(thatT,:))';
            thatNz = size(thatZ,1);

            logscoredraws = NaN(MCMCdraws,1); % note: for const, RB works over MCMCdraws, not Nfcstdraws
            for ii = 1 : MCMCdraws
                Zvar       = thatC * YstateTp1var(:,:,ii) * thatC';
                sqrtZvar   = chol(Zvar, 'lower');
                logdetvar  = 2 * sum(log(diag(sqrtZvar)));
                thatZtilde = sqrtZvar \ (thatZ - thatC * YstateTp1mean(:,ii));
                logscoredraws(ii) = - .5 * (thatNz * logtwopi + logdetvar + sum(thatZtilde.^2));
            end


            fcstYmvlogscoreDraws(thisT, :) = logscoredraws;
            maxlogscoredraw                = max(logscoredraws);
            fcstYmvlogscore(thisT)         = log(mean(exp(logscoredraws - maxlogscoredraw))) ...
                + maxlogscoredraw;
        end


        %% predictive densities for calendar years

        % prepare RTdata
        thisDate      = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI','TBILL'}
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
            mcmc_Ydensitydraws, RT_vec, NfcstDraws, Nbar);

        YBARdensityDraws(thisT,:,:)    = theseYBARdraws;
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
            YCURRENTdensityDraws(thisT,:,:)    = theseYCURRENTdraws;
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
        else
            YFINALdraws{thisT}      = [];
            ETAFINALdraws{thisT}    = [];
        end

        fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)


    end % parfor

    YFINALdraws     = YFINALdraws{T};
    ETAFINALdraws   = ETAFINALdraws{T};


    %% YRECESS (for quarterly 0:4)
    if strcmpi(datalabel, 'RGDP') 
        theseDraws  = (exp((YdensityDraws(:,1:5,:) ./ 100)) - 1) .* 100;
        % note: undoing the logs should be an abundance of caution; at the
        % zero-threshold both log- and simple changes should be identical
        fcstYrecess = mean(theseDraws < 0, 3) * 100;

        checkdiff(fcstYrecess, mean(YdensityDraws(:,1:5,:) < 0, 3) * 100);
    end

    %% patch YAVG: draws and future
    % note: the construction of YAVG could also be moved inside the parfor
    % loop (kept here for legacy reasons)

    % collect sequence of historical realizatins to complete the current year
    ylags = NaN(T,1,Kavg-1);
    parfor thisT = Tstart : T
      thisDate      = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI','TBILL'}
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
        
        % collect three lags, note that RT_obsdate_idx is one lag realtive to thisT
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

    % construct YAVGyfuture
    YAVGfuture = NaN(size(Yfuture));
    ylags      = permute(ylags, [1 3 2]);
    % h < Kavg
    for h = 1 : (Kavg - 1)
        YAVGfuture(:,h) = (sum(Yfuture(:,1:h), 2) + sum(ylags(:,end-(Kavg-1)+h:end), 2)) / Kavg;
    end
    % h >= Kavg
    for h = Kavg : size(YAVGfuture, 2)
        YAVGfuture(:,h) = sum(Yfuture(:,h+(-(Kavg-1):0)), 2) / Kavg;
    end

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

    %     fcstYAVGlogscoreRB = NaN(T,Nhorizons);

    fcstYAVGlogscore   = NaN(T,Nhorizons);
    fcstYAVGcrps        = NaN(T,Nhorizons);
    for thisT = Tstart : T
        thesedraws = permute(YAVGdensityDraws(thisT,:,:), [3 2 1]);
        for hh = Kavg : Nhorizons
            fcstYAVGcrps(thisT,hh)   = crpsDraws(YAVGfuture(thisT,hh), thesedraws(:,hh));
            fcstYlogscore(thisT,hh)  = log(ksdensity(thesedraws(:,hh), YAVGfuture(thisT,hh)));
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
    Ng = Ny * Ny;

    % prior
    Gprior = chol(G0var, 'lower') * randn(Ng, MCMCdraws);
    Gprior = reshape(Gprior, Ny, Ny, MCMCdraws);
    Gprior = permute(Gprior, [2 1 3]);

    Gpriorlambda = NaN(MCMCdraws, 1);
    parfor mm = 1 : MCMCdraws
        Gpriorlambda(mm) = max(abs(eig(Gprior(:,:,mm))));
    end
    Gpriorlambda = repmat(Gpriorlambda, 1, T);

    % posterior
    Glambda = NaN(MCMCdraws, T);
    parfor tt = Tstart : T
        for mm = 1 : MCMCdraws
            Glambda(mm,tt) = max(abs(eig(Gdraws(:,:,mm,tt))));
        end
    end

    thisfig = figure;
    hold on
    plot(dates, median(Gpriorlambda,1), 'k--', 'linewidth', 2)
    plot(dates, prctile(Gpriorlambda, [25 75], 1), 'k--')

    plot(dates, median(Glambda,1), 'b-', 'linewidth', 2)
    plot(dates, prctile(Glambda, [25 75], 1), 'b-')


    ylim([min([ylim, 0]) max(ylim)])
    xtickdates(dates([Tstart T]))
    wrapthisfigure(thisfig, sprintf('Geigenvalues-%s', modellabel), wrap)


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

    %% YBAR DRPS
    if ~isempty(PROBdata)
        for nn = 1 : Nbar
            thisfig = figure;
            hanni = plot(dates, [fcstYBARdrps(:,nn) fcstYBARcrps(:,nn) PROBdata.fcstSPFdrps(:,nn)]);
            xtickdates(dates(Tstart:T))
            legend(hanni, 'DRPS', 'CRPS', 'SPF-DRPS')
            hanni(1).LineWidth = 2;
            hanni(2).LineStyle = ':';
            hanni(2).LineWidth = 2;
            hanni(3).LineStyle = '-.';
            hanni(3).Color     = [0 0 0];
            hanni(3).LineWidth = 2;
            title(sprintf('h=%d', nn))
            wrapthisfigure(thisfig, sprintf('YBAR-DRPS-CRPS-h%d-%s', nn, modellabel), wrap)
        end
    end

    %% store results
    diary off
    if  doStore && ~quicky
        if doStoreXXL
            save(fullfile(localstore, modellabel), ...
                '-v7.3');
        end
        clear foo theseDraws
        close all
        matname = sprintf('CGMmcmc-%s-Ndraws%d', modellabel, MCMCdraws);
        save(fullfile(localstore, matname), '-v7.3');
        clear fcstYmvlogscoreDraws Y*density* % YBARdensityDraws YCURRENTdensityDraws
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
