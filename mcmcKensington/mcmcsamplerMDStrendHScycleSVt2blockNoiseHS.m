function [YdensityDraws, Ydraws, ETAdraws, ETASVdraws, ...
    SVdraws, SVtdraws, hVCV, RHO, tDOF, tSCALElog2, sqrtSIGMA, KslopesDraws, mustarGlobalSTATE, ...
    NOISEdraws, NOISEvarstate, ...
    Ztp1mean, Ztp1meanerror, Zmvlogscore  ...
    ] = mcmcsamplerMDStrendHScycleSVt2blockNoiseHS(Z, Znanny, Cz, Ngap, ...
    ndxZbar, dates, ...
    T, MCMCdraws, Nfedraws, rndStream, showProgress, thinMax)

%#ok<*UNRCH>
%#ok<*DATNM>
%#ok<*DATST>

if nargin < 11 || isempty(showProgress)
    showProgress = false;
end
if nargin < 12 || isempty(thinMax)
    thinMax = 1;
end

%% get dimensions and process Ngap setting
[Tall, Nz] = size(Z); % note: data matrix to be transposed later
Ny         = size(Cz,2);
sumCz      = sum(Cz,2);
if ~all(sumCz(~transpose(Znanny)) == 1)
    error('Cz gap loadings should sum to one')
end
if isempty(Ngap)
    Ngap = Ny;
end
Nstates = Ngap + 1; % Ytilde and trend
if Ngap < Ny
    Cz = Cz(:,1:Ngap,:);
end
Cz        = cat(2,Cz,ones(Nz,1,Tall)); % add trend loadings.
Nhorizons = max(Ny, Nstates - 1); % Ny includes lagged data, Nhorizons counts forecast steps

Nhorizons = 17;

maxShake  = 1e3;

%% get one-step ahead objects prior to cutting the sample
if T < Tall
    Znannytp1 = Znanny(T+1,:)';
    Nztp1     = sum(~Znannytp1);
    Ztp1      = Z(T+1,~Znannytp1)';
    Ctp1      = Cz(~Znannytp1,:,T+1);
    % Note: with BOP, we can have Nztp1 > Nstates.
    % In that case, logscore will only be computed for first Nstates entries of Ztp1
    
else
    Ztp1      = [];
    Znannytp1 = [];
    Ctp1      = [];
end



%% cut things down to sample of length T
Z         = Z(1:T,:)'; % note transpose
Znanny    = Znanny(1:T,:)';
Cz        = Cz(:,:,1:T);
Z(Znanny) = 0;
for t = 1 : T
    Cz(Znanny(:,t),:,t) = 0;
end

dates     = dates(1:T);


%% prepare parameters for noise processes
if ~islogical(ndxZbar)
    error('ndxZbar should be logical')
end
Nzbar     = sum(ndxZbar);
ndxZBAR   = repmat(ndxZbar, T, 1);

Nquarters       = 4;
years           = unique(year(dates));
Nyears          = length(years);
NyearsNquarters = Nyears * Nquarters;
datesYQ         = datenum(repelem(years,Nquarters), repmat(1:3:12,1,Nyears)', ones(NyearsNquarters,1));
ndxYQ           = ismember(datesYQ, dates); % maps datesYQ to dates


%% SV
Nsv         = 2;
ndxSV11     = 1:5;
ndxSV22     = ndxSV11(end)+1:Ngap; % OK, let's hard code this. Could use cell arrays for easier extension to N blocks, but leave that for later
svN11       = length(ndxSV11);
svN22       = length(ndxSV22);
I11         = eye(svN11);
I22         = eye(svN22);

%% get KSC values
[KSC, KSCt, logy2offset] = getKSC10values(T, Ngap); % using 10-point grid as in Omori et al (2007, JoE)
% note: block SV, but Ngap signals

%% prepare t
tdofGrid           = 3 : 40;
tdof.Ndof          = length(tdofGrid);
tdof.values        = transpose(tdofGrid(:)); % enforce row vector
tdof.logprior      = repmat(-log(tdof.Ndof), 1, tdof.Ndof); % uniform
% define constant part of loglik:
svNgap             = [svN11;svN22];
tdof.loglike0      = T * (gammaln(.5 * (tdofGrid + svNgap)) - gammaln(.5 * tdofGrid) -.5 .* svNgap .* log(tdofGrid .* pi));

%% priors
% scale SV
hvcvDof = Nsv;
hvcvT   = eye(Nsv);

% priors for each Sigma Block
Sigma11Dof    = svN11;
Sigma22Dof    = svN22;
Sigma11T0     = .01 .* eye(svN11);
Sigma22T0     = .01 .* eye(svN22);

Nslopes12       = svN22 * svN11;
E0slopes12      = zeros(Nslopes12,1);
V0slopes12      = eye(Nslopes12);
invV0slopes12   = V0slopes12 \ eye(Nslopes12);
invV0E0slopes12 = invV0slopes12 * E0slopes12;


rho0   = 0.8 .* ones(Nsv,1);
rhoV0i = (1 / 0.2^2) * eye(Nsv);

%% Noise variance: priors

maxNoiseState     = 1e4;
minNoiseState     = 1e-4;

%% prepare State Space
Fgap       = diag(ones(Ngap-1,1),1);
F          = blkdiag(Fgap, 1);
FF         = repmat(F, [1 1 T]); % note: no 3D sparse matrices


%% prepare precision-based sampler
% parameters
NzT   = Nz * T;
NsT   = Nstates * T;
NsT0  = NsT + Nstates;

% initial levels
Y00                      = zeros(Nstates, 1);

invcholSigmaY00          = (1 ./ 5) .* eye(Nstates);
invcholSigmaY00(end,end) = 1 ./ 100; % diffuse trend level Ytrend(0) (note Ytrend(-1) does not matter)
invcholSigmaY00vec       = invcholSigmaY00(:);

YY00  = sparse(1:Nstates, 1, Y00, NsT0, 1);


% AA
ndx1NsTp1 = (1 : NsT0)';
arows     = repmat((1 : Nstates)', 1 , Nstates, T);
arows     = arows + permute(Nstates * (1 : T), [1 3 2]);
rowndx    = [ndx1NsTp1; arows(:)];
colndx    = [ndx1NsTp1; vec(repmat(1 : NsT, Nstates, 1))];
values    = [ones(NsT0,1); -FF(:)];
AA        = sparse(rowndx, colndx, values);

% BB
brows1  = repmat((1 : Nstates)', 1 , Nstates);
brows2  = repmat((1 : Nstates)', 1 , Nstates, T);
brows2  = brows2 + permute(Nstates * (0 : T-1), [1 3 2]);
brows2  = Nstates + brows2;
brows   = [brows1(:); brows2(:)];

bcols1  = repmat(1 : Nstates, Nstates , 1);
bcols2  = Nstates + repmat(1 : NsT, Nstates, 1);
bcols   = [bcols1(:); bcols2(:)];


% CC
crows     = repmat((1 : Nz)', 1 , Nstates, T);
crows     = crows + permute(Nz * (0 : T-1), [1 3 2]);
ccols     = repmat(Nstates + (1 : NsT), Nz, 1);
CC        = sparse(crows(:), ccols(:), Cz(:), NzT, NsT0);


CCnonoise = CC(~ndxZBAR,:);
CCnoise   = CC(ndxZBAR,:);

% ZZ
ZZnonoise = Z(~ndxZBAR);
ZZnoise   = Z(ndxZBAR);

ZZnonoiseNanny  = Znanny(~ndxZBAR);
ZZnoiseNanny    = Znanny(ndxZBAR);

ZZnonoise     = ZZnonoise(~ZZnonoiseNanny);
CCnonoise     = CCnonoise(~ZZnonoiseNanny,:);
ZZnoise       = ZZnoise(~ZZnoiseNanny);
CCnoise       = CCnoise(~ZZnoiseNanny,:);


drows  = 1:length(ZZnoise);
dcols  = 1:length(ZZnoise);

% one-off computations for nonoise case
[N1, N2]   = size(CCnonoise);
N2         = N2 - N1;

[QQ,RR,PP]  = qr(CCnonoise'); % permuted QR for (slight) improvement in speed
QQ1         = QQ(:,1:N1)';
QQ2         = QQ(:,N1+1:end)';
RR1         = RR(1:N1,1:N1)';
if spamaxabs(CCnonoise -  PP * RR1 * QQ1) > 1e-8
    warning('qr off')
end

EY        = AA \ YY00;
EZnonoise = CCnonoise * EY;
Y1tilde   = RR1 \ (transpose(PP) * (ZZnonoise - EZnonoise));
QQY1tilde = QQ1' * Y1tilde;

if any(isnan(Y1tilde))
    error('NaN in Y1tilde')
end

%% allocate memory for MCMC

Ydraws      = zeros(Nstates, T, MCMCdraws);
ETAdraws    = NaN(Nstates, T, MCMCdraws);
ETASVdraws  = NaN(Nstates, T, MCMCdraws);

NOISEdraws    = NaN(Nzbar, T, MCMCdraws);
NOISEvarstate = NaN(Nzbar, T, MCMCdraws);
noiseGlobalSTATE = NaN(Nzbar, Nquarters, MCMCdraws);
noiseGlobalSCALE = NaN(Nzbar, Nquarters, MCMCdraws);
noiseLocalSTATE  = NaN(Nzbar, Nquarters, Nyears, MCMCdraws);
noiseLocalSCALE  = NaN(Nzbar, Nquarters, Nyears, MCMCdraws);

hh           = NaN(Nsv, T, MCMCdraws);
hVCV         = NaN(Nsv, Nsv, MCMCdraws);
RHO          = NaN(Nsv, MCMCdraws);
tDOF         = NaN(Nsv, MCMCdraws);
tSCALElog2   = NaN(Nsv, T, MCMCdraws);
KslopesDraws = NaN(svN11,svN22,MCMCdraws);
sqrtSIGMA    = zeros(Ngap,Ngap,MCMCdraws); % zeros, not NaN, to prepare upper-triangular structure

mustarGlobalSTATE  = NaN(1, MCMCdraws);
mustarGlobalSCALE  = NaN(1, MCMCdraws);
mustarLocalSTATE   = NaN(1, T, MCMCdraws);
mustarLocalSCALE   = NaN(1, T, MCMCdraws);
mustarVARSTATE     = NaN(1, T, MCMCdraws);

ndxNyTT     = NsT + (1:Nstates);
YstateTmean = NaN(Nstates, MCMCdraws);
YstateTvar  = NaN(Nstates, Nstates, MCMCdraws);

igamrnd  = @(alpha,beta) 1 ./ gamrnd(alpha, 1 ./ beta);



%% loop over MCMC steps
if showProgress
    progressbar(0)
end

MCMCburnin              = 2 * MCMCdraws;
mm                      = 0;
thisMCMCdraw            = 0;
thinStep                = 0;
maxRepeatedRejections   = 5;
countRepeatedRejections = 0; % count number of rejections *that immediately follow each other*

while thisMCMCdraw < MCMCdraws
    
    
    %% init loop
    if mm == 0
        % draw initial values
        hvcvPREV      = iwishdraw(hvcvT, hvcvDof, 1, [], rndStream);
        
        mustarGlobalScalePREV = igamrnd(.5,1);
        mustarGlobalStatePREV = igamrnd(.5,1 ./ mustarGlobalScalePREV);
        mustarLocalScalePREV  = igamrnd(.5,ones(1,T));
        mustarLocalStatePREV  = igamrnd(.5,1 ./ mustarLocalScalePREV);
        
        % simulate initial values from RW
        h0PREV        = zeros(Nsv,1);
        hPREV         = chol(hvcvPREV, "lower") * randn(rndStream, Nsv, T);
        hPREV         = h0PREV + cumsum(hPREV,2);
        % impose upper bound (for initial values)
        hPREV(hPREV > 5) = 5;
        
        % hPREV           = zeros(Nsv, T);
        % tscalelog2PREV  = zeros(Nsv,T);
        tdofPREV        = randi(rndStream, [3 40], Nsv,1); % only used for init
        tdofPREV        = repmat(tdofPREV, 1, T) .* .5;
        tscalelog2PREV  = igamdraw(tdofPREV, tdofPREV);
        
        rhoPREV       = rand(rndStream,Nsv,1) * 2 - 1;
        
        sqrtsigma11PREV = iwishdraw(I11, svN11, 1, [], rndStream);
        sqrtsigma22PREV = iwishdraw(I22, svN22, 1, [], rndStream);
        slopes12PREV    = E0slopes12; %  + chol(V0slopes12, 'lower') * randn(rndStream, Nslopes12, 1);
        slopes12PREV    = reshape(slopes12PREV, svN11, svN22);
        
        % noise
        noiseGlobalScalePREV = igamrnd(.5,ones(Nzbar,Nquarters));
        noiseGlobalStatePREV = igamrnd(.5,1 ./ noiseGlobalScalePREV);
        noiseLocalScalePREV  = igamrnd(.5,ones(Nzbar,Nquarters,Nyears));  % this includes irrelevant draws when sample doe not start in Q1 or does not end in Q4
        noiseLocalStatePREV  = igamrnd(.5,1 ./ noiseLocalScalePREV);
        
    end
    
    
    %% store last values prior to rejection sampling steps
    sqrtsigma11LAST = sqrtsigma11PREV;
    sqrtsigma22LAST = sqrtsigma22PREV;
    slopes12LAST    = slopes12PREV;
    hLAST           = hPREV;
    h0LAST          = h0PREV;
    tscalelog2LAST  = tscalelog2PREV;
    rhoLAST         = rhoPREV;
    
    %% init rejection sampling
    invSVt   = transpose(exp(-.5 .* (hPREV + tscalelog2PREV)));
    OK       = false; % flag for rejection sampling
    notOKtxt = '';
    
    %% precision-based state space sampler (with separate noise)
    % exploit scale SV structure
    invsqrtsigma11 = I11 / sqrtsigma11PREV;
    invbgapsv11    = invsqrtsigma11 .* permute(invSVt(:,1), [3 2 1]);
    
    invsqrtsigma22 = I22 / sqrtsigma22PREV;
    invbgapsv22    = invsqrtsigma22 .* permute(invSVt(:,2), [3 2 1]);
    
    invsigma12     = - invsqrtsigma11 * slopes12PREV;
    invbgapsv12    = invsigma12 .* permute(invSVt(:,1), [3 2 1]);
    
    invbgapsv                       = zeros(Ngap, Ngap, T);
    invbgapsv(ndxSV11, ndxSV11, :)  = invbgapsv11;
    invbgapsv(ndxSV11, ndxSV22, :)  = invbgapsv12;
    invbgapsv(ndxSV22, ndxSV22, :)  = invbgapsv22;
    
    invbbsv                         = zeros(Nstates, Nstates, T);
    invbbsv(1:Ngap,1:Ngap,:)        = invbgapsv;
    
    thisVARSTATE          = permute(mustarGlobalStatePREV .* mustarLocalStatePREV, [3 1 2]);
    invbbsv(end, end, :)  = 1 ./ sqrt(thisVARSTATE);
    
    bvalues          = [invcholSigmaY00vec; invbbsv(:)];
    invBB            = sparse(brows, bcols, bvalues);
    
    
    AAtilde       = invBB * AA;
    AAtildeQQY1   = AAtilde * QQY1tilde;
    AAtildeQQ2    = AAtilde * QQ2';
    P22          = transpose(AAtildeQQ2) * AAtildeQQ2;
    Pmu22        = -AAtildeQQ2' * AAtildeQQY1;
    
    
    % conventional PBS with noise for Znoise
    % noise terms
    noiseState       = noiseGlobalStatePREV .* noiseLocalStatePREV;
    noiseState       = reshape(noiseState, Nzbar, NyearsNquarters);
    noiseState       = noiseState(:,ndxYQ);
    noiseState(noiseState > maxNoiseState) = maxNoiseState;
    noiseState(noiseState < minNoiseState) = minNoiseState;
    dvalues         = 1 ./ sqrt(noiseState);
    invDD           = sparse(drows, dcols, dvalues(~ZZnoiseNanny));
    
    CCnoisetilde     = invDD * CCnoise * QQ2';
    ZZnoisetilde     = invDD * (ZZnoise - CCnoise * QQY1tilde);
    P22bar           = P22   + CCnoisetilde' * CCnoisetilde;
    % choleski
    [cholP22bar, flag] = chol(P22bar, 'lower');
    if flag == 0
        OK = true;
    else
        notOKtxt = 'PSchol';
    end
    
    if OK
        
        Pmu22bar         = Pmu22 + CCnoisetilde' * ZZnoisetilde;
        y2hat            = cholP22bar \ Pmu22bar;
        y2hat            = transpose(cholP22bar) \ y2hat;
        YHAT             = EY + QQ1' * Y1tilde + QQ2' * y2hat;
        thisYhat         = YHAT(ndxNyTT);
        yvol             = cholP22bar \ QQ2(:,ndxNyTT); % note: yvol is not square!
        thisYvar         = yvol' * yvol;
        
        y2draw           = transpose(cholP22bar) \ randn(rndStream, N2, 1);
        YDRAW            = YHAT + QQ2' * y2draw;
        thisY            = reshape(YDRAW(Nstates+1:NsT0), Nstates, T);
        
        ETADRAW       = AA * YDRAW - YY00;
        eta           = reshape(ETADRAW(Nstates+1:NsT0), Nstates, T);
        
        residGAP      = transpose(eta(1:Ngap,:));
        mustar        = eta(end,:);
        
        noisevec                  = ZZnoise - CCnoise * YDRAW;
        thisNoise                 = NaN(Nzbar, T);
        thisNoise(~ZZnoiseNanny)  = noisevec;
        noise                     = NaN(Nzbar,NyearsNquarters);
        noise(:,ndxYQ)            = thisNoise;
        noise                     = reshape(noise, Nzbar, Nquarters, Nyears);
        
    end % OK precsam
    
    
    
    if OK
        %% BLOCK SV PART
        
        % draw slopes12
        lhs                            = residGAP(:,ndxSV11) .* invSVt(:,1);
        rhs                            = residGAP(:,ndxSV22) .* invSVt(:,1); % scale by SV1 for purpose of regression
        invsigma11                     = transpose(invsqrtsigma11) * invsqrtsigma11;
        [slopes12PREV, resid11scaled]  = bayesVectorRegressionGibbsDraw1(lhs, rhs, invsigma11, invV0E0slopes12, invV0slopes12, rndStream);
        slopes12PREV                   = transpose(slopes12PREV);
        
        % draw constant VCV11
        sqrtsigma11PREV = bayesCHOLVCVgibbsDraw1(Sigma11T0, Sigma11Dof, resid11scaled, rndStream);
        
        % draw constant VCV22
        resid22scaled   = residGAP(:,ndxSV22) .* invSVt(:,2);
        sqrtsigma22PREV = bayesCHOLVCVgibbsDraw1(Sigma22T0, Sigma22Dof, resid22scaled, rndStream);
        
        % common SV draw
        z11        = sqrtsigma11PREV \ transpose(resid11scaled ./ invSVt(:,1));
        z22        = sqrtsigma22PREV \ transpose(residGAP(:,ndxSV22));
        y2         = [z11; z22].^2;
        logy2      = log(y2 + logy2offset);
        
        [hPREV, ~, h0PREV, ~, tscalelog2PREV, tdofPREV] = Block2StochVoltAR1(y2, logy2, svN11, ndxSV11, svN22, ndxSV22, hPREV, rhoPREV, chol(hvcvPREV, "lower"), KSC, KSCt, tdof, Ngap, T, rndStream);
        
        %% SV-AR1 rho draw
        % TODO: consider storing chol(hvcvPREV, "lower") ?
        hLAG     = [h0PREV, hPREV(:,1:end-1)];
        rhodraws = bayesAR1SURdraw(hPREV', hLAG', hvcvPREV, rho0, rhoV0i, maxShake, rndStream);
        shake    = 0;
        thisOK   = false;
        while ~thisOK && shake < maxShake
            shake     = shake + 1;
            thisOK    = all(abs(rhodraws(:,shake)) < 1);
        end
        if thisOK
            rhoPREV = rhodraws(:,shake);
        else
            notOKtxt = 'unstableSV';
        end
    end
    
    OK = OK && thisOK;
    
    if OK
        
        %% draw hvcv
        hshock   = transpose(hPREV - rhoPREV .* hLAG);
        hvcvPREV = bayesVCVgibbsDraw1(hvcvT, hvcvDof, hshock,rndStream);
        
        %% draw horseshoe components of mustarvar
        mustar2    = mustar.^2;
        [mustarGlobalStatePREV, mustarGlobalScalePREV, mustarLocalStatePREV, mustarLocalScalePREV] = ...
            horseshoePosteriorDraw(mustar2, mustarGlobalStatePREV, mustarGlobalScalePREV, ...
            [], mustarLocalScalePREV, 2);
        
        
        
        %% draw horseshoe components of noisevol
        noise2    = noise.^2;
        [noiseGlobalStatePREV, noiseGlobalScalePREV, noiseLocalStatePREV, noiseLocalScalePREV] = ...
            horseshoePosteriorDraw(noise2, noiseGlobalStatePREV, noiseGlobalScalePREV, ...
            [], noiseLocalScalePREV, 3);
        
        %% move on and store draws post burnin
        mm                      = mm + 1;
        countRepeatedRejections = 0; % count number of rejections *that immediately follow each other*
        
        if mm > MCMCburnin
            
            thinStep = thinStep + 1;
            if thinStep == thinMax
                thinStep                        = 0;
                thisMCMCdraw                    = thisMCMCdraw + 1;
                
                Ydraws(:,:,thisMCMCdraw)        = thisY;
                
                ETAdraws(:,:,thisMCMCdraw)      = eta;
                hh(:,:,thisMCMCdraw)            = hPREV;
                hVCV(:,:,thisMCMCdraw)          = hvcvPREV;
                RHO(:,thisMCMCdraw)             = rhoPREV;
                tDOF(:,thisMCMCdraw)            = tdofPREV;
                tSCALElog2(:,:,thisMCMCdraw)    = tscalelog2PREV;
                
                mustarGlobalSTATE(:,thisMCMCdraw)   = mustarGlobalStatePREV;
                mustarGlobalSCALE(:,thisMCMCdraw)   = mustarGlobalScalePREV;
                mustarLocalSTATE(:,:,thisMCMCdraw)  = mustarLocalStatePREV;
                mustarLocalSCALE(:,:,thisMCMCdraw)  = mustarLocalScalePREV;
                mustarVARSTATE(:,:,thisMCMCdraw)    = mustarGlobalStatePREV .* mustarLocalStatePREV;
                
                KslopesDraws(:,:,thisMCMCdraw)          = slopes12PREV;
                sqrtSIGMA(ndxSV11,ndxSV11,thisMCMCdraw) = sqrtsigma11PREV;
                sqrtSIGMA(ndxSV11,ndxSV22,thisMCMCdraw) = slopes12PREV * sqrtsigma22PREV;
                sqrtSIGMA(ndxSV22,ndxSV22,thisMCMCdraw) = sqrtsigma22PREV;
                
                
                noiseGlobalSTATE(:,:,thisMCMCdraw)     = noiseGlobalStatePREV;
                noiseGlobalSCALE(:,:,thisMCMCdraw)     = noiseGlobalScalePREV;
                noiseLocalSTATE(:,:,:,thisMCMCdraw)    = noiseLocalStatePREV;
                noiseLocalSCALE(:,:,:,thisMCMCdraw)    = noiseLocalScalePREV;
                noiseState                             = reshape(noiseGlobalStatePREV .* noiseLocalStatePREV, Nzbar, NyearsNquarters);
                % noiseState(noiseState > maxNoiseState) = maxNoiseState;
                NOISEvarstate(:,:,thisMCMCdraw)        = noiseState(:,ndxYQ);
                NOISEdraws(:,:,thisMCMCdraw)           = thisNoise;

                
                YstateTmean(:,thisMCMCdraw)  = thisYhat;
                YstateTvar(:,:,thisMCMCdraw) = thisYvar;
            end
        end
        
        if showProgress
            progressbar(mm / (MCMCburnin + MCMCdraws * thinMax))
        end
        
    else % not OK
        
        % prepare to redo this MCMC step
        
        warning('Rejected MCMC step %4d, (T=%d, count=%d): %s', mm, T, countRepeatedRejections, notOKtxt)
        % note: if sampler gets stuck here, need to store stack of multiple past steps to walk back
        
        rhoPREV         = rhoLAST;
        hPREV           = hLAST;
        h0PREV          = h0LAST;
        tscalelog2PREV  = tscalelog2LAST;
        sqrtsigma11PREV = sqrtsigma11LAST;
        sqrtsigma22PREV = sqrtsigma22LAST;
        slopes12PREV    = slopes12LAST;
        countRepeatedRejections = countRepeatedRejections + 1;
        
        if countRepeatedRejections > maxRepeatedRejections
            warning('Too many repeated rejections (MCMC step %d, T=%d), restarting MCMC chain', mm, T) % if this occurs too often, need to track a stack of past steps to walk back
            mm           = 0;
            thisMCMCdraw = 0;
        end
    end
end

%% collect SV
SVdraws  = squeeze(exp(hh * 0.5));
SVtdraws = squeeze(exp((hh + tSCALElog2) .* 0.5));
for mm = 1 : MCMCdraws
    for jj = 1 : T
        
        thisBgap               = sqrtSIGMA(:,:,mm);
        thisBgap(:,ndxSV11)    = thisBgap(:,ndxSV11) .* SVtdraws(1,jj,mm);
        thisBgap(:,ndxSV22)    = thisBgap(:,ndxSV22) .* SVtdraws(2,jj,mm);
        
        thisB                  = zeros(Nstates,Nstates);
        thisB(Nstates,Nstates) = sqrt(mustarVARSTATE(1,jj,mm));
        thisB(1:Ngap,1:Ngap)   = thisBgap;
        ETASVdraws(:,jj,mm)    = sqrt(diag(thisB * thisB'));
    end
end

%% simulations: allocate memory (note: will all be reshaped and permuted as needed below)

% Ytp1 stores first and second moment of term structure vector Y (not the trend-cycle state vector)
Ytp1mean  = NaN(Nstates, MCMCdraws);
if ~isempty(Ztp1)
    Ytp1var      = NaN(Nstates, Nstates, Nfedraws, MCMCdraws);
    Ztp1noisevar = zeros(Nztp1, Nfedraws, MCMCdraws); % stores only the diagonals
end
YdensityDraws = NaN(2 * Nfedraws, Nhorizons, MCMCdraws); % antithetic simulation for linear shocks

% specify bounds for horseshoe state draws
maxTrendVarState = 20^2;

%% run simulations for each MCMCdraw
for mm = 1 : MCMCdraws
    
    %% collect MCMC parameter draws
    rho           = RHO(:,mm);
    hT            = hh(:,end,mm);
    hsqrtvcv      = chol(hVCV(:,:,mm), "lower");
    tdof          = tDOF(:,mm);
    Bgap          = sqrtSIGMA(:,:,mm);
    thisY         = Ydraws(:,end,mm);
    
    %% draw future mustar variance states from horseshoe
    mustarVariances    = mustarGlobalSTATE(:,mm) .* halfcauchydraws(1, Nfedraws, Nhorizons);
    if any(mustarVariances > maxTrendVarState, 'all')
        % warning('out-of-sample outlier draws exceed bounds (trim applied,  T=%d)', T)
    end
    % rule out highly extreme draws
    mustarVariances(mustarVariances > maxTrendVarState) = maxTrendVarState;
    
    %% simulate SVt
    SVt = simulateSVt(rho, hT, hsqrtvcv, tdof, Nfedraws, Nhorizons, rndStream);
    
    %% one-step ahead moments of termstructure
    % mean
    if ~isempty(Ztp1)
        Ytp1mean(:,mm)    = F * YstateTmean(:,mm);
        Ytp1var(:,:,:,mm) = simulateYtp1varSV2block(ndxSV11, ndxSV22, Nstates, Nfedraws, YstateTvar(:,:,mm), Bgap, SVt, mustarVariances(:,:,1), F);
        
        qNext = quarter(dates(T)) + 1; % note: dates cut to length T already
        if qNext > 4
            qNext = 1;
        end
        noiseVariances = noiseGlobalSTATE(:,qNext,mm) .* halfcauchydraws(Nzbar, Nfedraws);
        noiseVariances(noiseVariances > maxNoiseState) = maxNoiseState;
        
        for ii = 1 : Nfedraws
            tmp                   = zeros(Nz,1);
            tmp(ndxZbar)          = noiseVariances(:,ii);
            Ztp1noisevar(:,ii,mm) = tmp(~Znannytp1);
        end
    end
    
    %% simulate density of Y
    YdensityDraws(:,:,mm) = simulateMDS2block(ndxSV11, ndxSV22, Nstates, Nfedraws, Nhorizons, Bgap, SVt, mustarVariances, F, thisY, rndStream);
    
end

%% logscore and mean for one-step ahead data vector
if ~isempty(Ztp1)
    % mean (per MCMC node)
    Ztp1mean      = Ctp1 * Ytp1mean;
    Ztp1meanerror = Ztp1 - Ztp1mean;
    % logscore
    Zmvlogscore = calculateLogscore(Ytp1var, Ctp1, Ztp1meanerror, Ztp1noisevar);
    % mean (integrated across MCMC)
    Ztp1mean      = mean(Ztp1mean,2);
    Ztp1meanerror = mean(Ztp1meanerror,2);
    
else
    % defaults: NaN instead of empty to facilitate assignment into vectors
    Zmvlogscore   = NaN;
    Ztp1mean      = NaN;
    Ztp1meanerror = NaN;
    
end


%% collect output and transform Ydraws from trend-cycle form into term structure
% ETA and ETASVdraws remain defined for trend-cycle states
% YdensityDraws are already in term structure form
Ndraws        = 2 * Nfedraws * MCMCdraws;
YdensityDraws = permute(YdensityDraws, [2 1 3]);
YdensityDraws = reshape(YdensityDraws, Nhorizons, Ndraws);
% add trend to all elements of Ydraws
Ydraws(1:Nstates-1,:,:) = Ydraws(1:Nstates-1,:,:) + Ydraws(Nstates,:,:);

end


