function [YdensityDraws, Ydraws, YhatRB, ETAdraws, RESIDSVdraws,ZRESIDdraws, ...
    Gdraws, Glambda, G0var, YgapCONSTdraws, ...
    SVdraws, SVtdraws, hVCV, RHO, tDOF, tSCALElog2, sqrtSIGMA, KslopesDraws, mustarGlobalSTATE, ...
    Ztp1mean, Ztp1meanerror, Zmvlogscore  ...
    ] = mcmcsamplerVAR0trendHScycleSVt2block(Z, Znanny, Cz, Ngap, ...
    G0var, ...
    T, MCMCdraws, Nfedraws, rndStream, showProgress)

%#ok<*UNRCH>
if nargin < 10
    showProgress = false;
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
thinMax   = 1;

%% get one-step ahead objects prior to cutting the sample
if T < Tall
    Znannytp1 = Znanny(T+1,:)';
    % Nztp1     = sum(~Znannytp1);
    Ztp1      = Z(T+1,~Znannytp1)';
    Ctp1      = Cz(~Znannytp1,:,T+1);
    % Note: with BOP, we can have Nztp1 > Nstates.
    % In that case, logscore will only be computed for first Nstates entries of Ztp1
    
else
    Ztp1      = [];
    Znannytp1 = []; %#ok<NASGU>
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
Igap        = eye(Ngap);


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



%% prepare State Space: Y transition
Fgap       = diag(ones(Ngap-1,1),1);
F          = blkdiag(Fgap, 1);

IgapMinusFgap = Igap - Fgap;

%% prepare State Space: eta and eta-tilde transition
G0         = zeros(Ngap); % dummy values
Ng         = numel(G0);

% minnesota-style prior (assuming equal residuals variances)
if isempty(G0var)
    G0var  = NaN(Ng,1); % diagonal elements
    Gtheta    = [.2^2, 0.5^2];
    ndx = 0;
    for i = 1 : Ngap
        for j = 1 : Ngap
            ndx  = ndx + 1;
            if (i==j)
                G0var(ndx)=Gtheta(1);
            else
                G0var(ndx)=Gtheta(1)*Gtheta(2);
            end
        end
    end
    invG0var   = diag(1 ./ G0var);
else
    invG0var = eye(Ng) / G0var;
end

invG0varG0 = invG0var * G0(:);




%% prepare precision-based sampler
% State vector is S(t) = [Ygap(t); Ytrend(t)];
% initial state vector loads on two lags of Ygap, i.e. Ygap(0) and Ygap(-1) (lagged Ytrend is superfluous; fixed to zero in model)

% parameters
NzT      = Nz * T;
NsT      = Nstates * T;
Nstates2 = Nstates * 2; % elements of VAR state vector (note: contains one superfluous lag of Ytrend)
N0       = Nstates2; % time 0 states
Nconst   = Ngap;
NSTATES  = NsT + N0 + Nconst;

% initial levels
invcholSigmaY00          = (1 ./ 5) .* eye(N0);
invcholSigmaY00(1,1)     = 1 ./ 100; % irrelevant
invcholSigmaY00(Nstates+1,Nstates+1) = 1 ./ 100;
invcholSigmaY00vec       = invcholSigmaY00(:);

YY00  = sparse(NSTATES, 1); % all zeros

cholSigmaYconst          = eye(Nconst); % placeholder
for nn=1:Nconst
    cholSigmaYconst(nn,nn) = 5/sqrt(nn); % linear decay
end
invcholSigmaYconst       = Igap / cholSigmaYconst;
invcholSigmaYconstvec    = invcholSigmaYconst(:);

ndx1NSTATES = (1 : NSTATES)';
onesNSTATES = ones(NSTATES,1);

ndxYgapCONST = NSTATES - Nconst + (1:Nconst);
ndxY         = N0+1:NSTATES-Nconst;

%% AAy as in AAy * Y = Y0 + ETA (note: tracking also the initial lag)
rowsF    = repmat((1 : Nstates)', 1 , Nstates, T+1);
rowsF    = rowsF + permute(Nstates * (0 : T), [1 3 2]);
rowsF    = rowsF + Nstates; % to account for Y(-1)

rowsCONST    = repmat((1 : Nstates)', 1 , Ngap, T+1); % note: only Ngap columns
rowsCONST    = rowsCONST + permute(Nstates * (0 : T), [1 3 2]);
rowsCONST    = rowsCONST + Nstates; % to account for Y00(-1)

rowsAAy  = [ndx1NSTATES; rowsF(:); rowsCONST(:)];

colsF     = repmat(1 : NsT + Nstates, Nstates, 1);
colsCONST = (NSTATES - Ngap) + repmat(1:Ngap,Nstates,1,T+1); % 3D to generate correct ordering
colsAAy   = [ndx1NSTATES; colsF(:); colsCONST(:)];

FF            = repmat(F, [1 1 T+1]); % T+1 to generate also ETAlag

GGCONST           = zeros(Nstates,Ngap);  % to capture loadings on constants
GGCONST(1:Ngap,:) = IgapMinusFgap;
aayCONST          = repmat(GGCONST, [1 1 T+1]);
avalues           = [onesNSTATES; -FF(:); -aayCONST(:)];
AAy               = sparse(rowsAAy, colsAAy, avalues, NSTATES, NSTATES);


%% AA as in AA * Y = Y0 + RESID
rowsAAlags     = repmat((1 : Nstates)', 1 , Nstates, T);
rowsAAlags     = rowsAAlags + permute(Nstates * (0 : T-1), [1 3 2]);
rowsAAlags     = rowsAAlags + N0; % to account for Y00(-1)
rowsAA         = [ndx1NSTATES; rowsAAlags(:); rowsAAlags(:); rowsCONST(:)];

colsAAlag2  = repmat(1 : NsT, Nstates, 1);
colsAAlag1  = colsAAlag2 +  Nstates; % to account for Y00(-1)
colsAA      = [ndx1NSTATES; colsAAlag1(:); colsAAlag2(:); colsCONST(:)];

% sortndx for AA
ndx           = sub2ind([NSTATES, NSTATES], rowsAA, colsAA);
[~,sortndxAA] = sort(ndx);
rowsAA        = rowsAA(sortndxAA);
colsAA        = colsAA(sortndxAA);

%% BB
brows1  = repmat((1 : N0)', 1 , N0);
brows2  = repmat((1 : Nstates)', 1 , Nstates, T);
brows2  = brows2 + permute(Nstates * (0 : T-1), [1 3 2]);
brows2  = brows2 + N0;
brows3  = repmat((1 : Nconst)', 1 , Nconst);
brows3  = brows3 + NSTATES - Nconst;
brows   = [brows1(:); brows2(:); brows3(:)];

bcols1  = repmat(1 : N0, N0 , 1);
bcols2  = repmat(1 : NsT, Nstates, 1) + N0;
bcols3  = repmat(1 : Nconst, Nconst , 1) + NSTATES - Nconst;
bcols   = [bcols1(:); bcols2(:); bcols3(:)];

% CC
crows     = repmat((1 : Nz)', 1 , Nstates, T);
crows     = crows + permute(Nz * (0 : T-1), [1 3 2]);
ccols     = repmat(N0 + (1 : NsT), Nz, 1);
CC        = sparse(crows(:), ccols(:), Cz(:), NzT, NSTATES);

% ZZ
ZZ        = Z(:);
ZZnanny   = Znanny(:);

ZZ        = ZZ(~ZZnanny);
CC        = CC(~ZZnanny,:);


% one-off computations
[N1, N2]  = size(CC);
N2        = N2 - N1;

[QQ,RR,PP]  = qr(CC'); % permuted QR for (slight) improvement in speed
QQ1         = QQ(:,1:N1)';
QQ2         = QQ(:,N1+1:end)';
RR1         = RR(1:N1,1:N1)';
if spamaxabs(CC -  PP * RR1 * QQ1) > 1e-8
    warning('qr off')
end

%% allocate memory for MCMC

Ydraws         = zeros(Nstates, T, MCMCdraws);
ETAdraws       = NaN(Nstates, T, MCMCdraws);
RESIDSVdraws   = NaN(Nstates, T, MCMCdraws);
ZRESIDdraws    = NaN(Nstates, T, MCMCdraws);

Gdraws         = NaN(Ngap, Ngap, MCMCdraws);
Glambda        = NaN(1, MCMCdraws);
YgapCONSTdraws = NaN(Ngap, MCMCdraws);

hh           = NaN(Nsv, T, MCMCdraws);
hVCV         = NaN(Nsv, Nsv, MCMCdraws);
RHO          = NaN(Nsv, MCMCdraws);
tDOF         = NaN(Nsv, MCMCdraws);
tSCALElog2   = NaN(Nsv, T, MCMCdraws);
KslopesDraws = NaN(svN11,svN22,MCMCdraws);
sqrtSIGMA    = zeros(Ngap,Ngap,MCMCdraws); % zeros, not NaN, to prepare upper-triangular structure

% mustarSIG          = NaN(1, MCMCdraws);
mustarGlobalSTATE  = NaN(1, MCMCdraws);
mustarGlobalSCALE  = NaN(1, MCMCdraws);
mustarLocalSTATE   = NaN(1, T, MCMCdraws);
mustarLocalSCALE   = NaN(1, T, MCMCdraws);
mustarVARSTATE     = NaN(1, T, MCMCdraws);

ndxStateT  = NSTATES - Nconst - Nstates2 + (1:Nstates2); % note: includes one superfluous lag of Ytrend
StateTmean = NaN(Nstates2, MCMCdraws);
StateTvar  = NaN(Nstates2, Nstates2, MCMCdraws);

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
        
        GGPREV    = zeros(Nstates); % We mock up a G matrix for the entire matrix (inclusive of trend), but only use the Ngap x Ngap block
        GPREV     = G0;
        
        GGCONST   = zeros(Nstates,Ngap);  % to capture loadings on constants
        
    end
    
    %% store last values prior to rejection sampling steps
    sqrtsigma11LAST = sqrtsigma11PREV;
    sqrtsigma22LAST = sqrtsigma22PREV;
    slopes12LAST    = slopes12PREV;
    hLAST           = hPREV;
    h0LAST          = h0PREV;
    tscalelog2LAST  = tscalelog2PREV;
    rhoLAST         = rhoPREV;
    GLAST           = GPREV;
    
    %% init rejection sampling
    invSVt     = transpose(exp(-.5 .* (hPREV + tscalelog2PREV)));
    OK       = false; % flag for rejection sampling
    notOKtxt = '';
    
    %% precision-based state space sampler
    
    
    % AA
    GGPREV(1:Ngap,1:Ngap) = GPREV;
    GGCONST(1:Ngap,:)     = (Igap - GPREV) * IgapMinusFgap;
    aaLag1                = repmat(GGPREV + F, [1 1 T]);
    aaLag2                = -repmat(GGPREV * F, [1 1 T]);
    aaCONST               = repmat(GGCONST, [1 1 T+1]);
    avalues               = [onesNSTATES; -aaLag1(:); -aaLag2(:); -aaCONST(:)];
    avalues               = avalues(sortndxAA);
    AA                    = sparse(rowsAA, colsAA, avalues, NSTATES, NSTATES);
    
    EY        = AA \ YY00;
    EZ        = CC * EY;
    Y1tilde   = RR1 \ (transpose(PP) * (ZZ - EZ));
    QQY1tilde = QQ1' * Y1tilde;
    
    if any(isnan(Y1tilde))
        error('NaN in Y1tilde')
    end
    
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
    
    bvalues          = [invcholSigmaY00vec; invbbsv(:); invcholSigmaYconstvec];
    invBB            = sparse(brows, bcols, bvalues);
    
    AAtilde       = invBB * AA;
    AAtildeQQY1   = AAtilde * QQY1tilde;
    AAtildeQQ2    = AAtilde * QQ2';
    invQSIG22     = transpose(AAtildeQQ2) * AAtildeQQ2;
    
    % choleski
    [cholinvQSIG22, flag] = chol(invQSIG22, 'lower');
    if flag == 0
        OK = true;
    else
        notOKtxt = 'PSchol';
    end
    
    if OK
        % compute y2hat and y2draw separately to store mean
        y2hat         = - cholinvQSIG22 \ (AAtildeQQ2' * AAtildeQQY1);
        y2hat         = transpose(cholinvQSIG22) \ y2hat;
        
        N2draws       = randn(rndStream, N2, 1);
        y2draw        = transpose(cholinvQSIG22) \ N2draws;
        YHAT          = EY + QQ1' * Y1tilde + QQ2' * y2hat;
        YDRAW         = YHAT + QQ2' * y2draw;
        
        
        thisShat      = YHAT(ndxStateT);
        statevol      = cholinvQSIG22 \ QQ2(:,ndxStateT); % note: statevol is not square!
        thisSvar      = statevol' * statevol;
        
        thisYconst      = YDRAW(ndxYgapCONST);
        thisY           = reshape(YDRAW(ndxY), Nstates, T); % Note: includes constants
        
        ETADRAW       = AAy * YDRAW;
        eta           = reshape(ETADRAW(ndxY), Nstates, T);
        etagap0       = ETADRAW(Nstates+(1:Ngap));
        etaGAP        = transpose(eta(1:Ngap,:));
        etaGAPlag     = transpose([etagap0, eta(1:Ngap,1:end-1)]);
        mustar        = eta(end,:);
        
    end % OK precsam
    
    if OK
        
        %% draw VAR for eta
        % for now: standard call to system estimation with SV
        
        % construct inverse VCV
        for tt = 1 : T
            invbgapsv(:,:,tt) = transpose(invbgapsv(:,:,tt)) * invbgapsv(:,:,tt);
        end
        
        [GPREV, residGAP]  = bayesVectorSVregressionGibbsDraw1(etaGAP, etaGAPlag, invbgapsv, ...
            invG0varG0, invG0var, rndStream);
        
        GPREV   = transpose(GPREV);
        maxabsG = max(abs(eig(GPREV)));
        if  maxabsG < 1
            thisOK = true;
        else
            thisOK   = false;
            notOKtxt = 'unstableVAR';
        end
    end % OK VAR
    
    OK = OK && thisOK;
    
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
            
        
        %% move on and store draws post burnin
        mm                      = mm + 1;
        countRepeatedRejections = 0; % count number of rejections *that immediately follow each other*
        
        if mm > MCMCburnin
            
            thinStep = thinStep + 1;
            if thinStep == thinMax
                thinStep                        = 0;
                thisMCMCdraw                    = thisMCMCdraw + 1;
                
                Ydraws(:,:,thisMCMCdraw)        = thisY;
                YgapCONSTdraws(:,thisMCMCdraw)  = thisYconst;
                
                ETAdraws(:,:,thisMCMCdraw)      = eta;
                ZRESIDdraws(:,:,thisMCMCdraw)   = [transpose(residGAP); mustar]; % will be scaled later scaled
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
                
                Gdraws(:,:,thisMCMCdraw)    = GPREV;
                Glambda(:,thisMCMCdraw)     = maxabsG;
                
                StateTmean(:,thisMCMCdraw)  = thisShat;
                StateTvar(:,:,thisMCMCdraw) = thisSvar;
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
        GPREV           = GLAST;
        
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
        RESIDSVdraws(:,jj,mm)  = sqrt(diag(thisB * thisB'));
    end
end
ZRESIDdraws = ZRESIDdraws ./ RESIDSVdraws;

%% simulations: allocate memory (note: will all be reshaped and permuted as needed below)

% Ytp1 stores first and second moment of term structure vector Y (not the trend-cycle state vector)
Ytp1mean = NaN(Nstates, MCMCdraws);
Ytp1var  = NaN(Nstates, Nstates, Nfedraws, MCMCdraws);

YdensityDraws = NaN(2 * Nfedraws, Nhorizons, MCMCdraws); % antithetic simulation for linear shocks
YhatRB        = NaN(Nhorizons, MCMCdraws);

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
    G             = Gdraws(:,:,mm);
    Y             = Ydraws(:,end,mm);
    Yconst        = YgapCONSTdraws(:,mm);
    ETA           = ETAdraws(:,end,mm);
    
    GG                = zeros(Nstates); % We mock up a G matrix for the entire matrix (inclusive of trend), but only use the Ngap x Ngap block
    GG(1:Ngap,1:Ngap) = G;
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
        Ycompanion           = [- F * GG, F + GG];
        thisYmean            = Ycompanion * StateTmean(:,mm);
        thisYmean(1:Ngap)    = thisYmean(1:Ngap) + (Igap - G) * (IgapMinusFgap) * Yconst;
        
        Ytp1mean(:,mm)       = thisYmean;
        Ytp1var(:,:,:,mm)    = simulateYtp1varSV2block(ndxSV11, ndxSV22, Nstates, Nfedraws, StateTvar(:,:,mm), Bgap, SVt, mustarVariances(:,:,1), Ycompanion);
    end
    
    
    %% simulate density of Y
    [YdensityDraws(:,:,mm), YhatRB(:,mm)] = simulateVAR02block(ndxSV11, ndxSV22, Nstates, Nfedraws, Nhorizons, Bgap, SVt, mustarVariances, F, Y, GG, Yconst, ETA, rndStream);
    
end


%% logscore and mean for one-step ahead data vector
if ~isempty(Ztp1)
    % mean (per MCMC node)
    Ztp1mean      = Ctp1 * Ytp1mean;
    Ztp1meanerror = Ztp1 - Ztp1mean;
    % logscore
    Zmvlogscore = calculateLogscore(Ytp1var, Ctp1, Ztp1meanerror);
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
% ETA and RESIDSVdraws remain defined for trend-cycle states
% YdensityDraws are already in term structure form
Ndraws        = 2 * Nfedraws * MCMCdraws;
YdensityDraws = permute(YdensityDraws, [2 1 3]);
YdensityDraws = reshape(YdensityDraws, Nhorizons, Ndraws);
% add trend to all elements of Ydraws
Ydraws(1:Nstates-1,:,:) = Ydraws(1:Nstates-1,:,:) + Ydraws(Nstates,:,:);
Ydraws(Nstates,:,:)     = Ydraws(Nstates,:,:) + permute(YgapCONSTdraws(1,:), [1 3 2]); % align trend with constant of outcome process

end
