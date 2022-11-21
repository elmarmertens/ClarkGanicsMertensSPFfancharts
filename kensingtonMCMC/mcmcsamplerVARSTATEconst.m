function [YdensityDraws, Ydraws, YhatRB, ETAdraws, ...
    Gdraws, ...
    sqrtSIGMA, ...
    YstateTp1mean, YstateTp1var, yTpHmean, yTpHvar  ...
    ] = mcmcsamplerVARSTATEconst(Z, Znanny, Cz, ...
    G0var, ...
    T, MCMCdraws, Nfedraws, rndStream, showProgress)

%#ok<*UNRCH>
if nargin < 8
    showProgress = false;
end


%% cut sample, set parameters
Z         = Z(1:T,:)'; % note transpose
Znanny    = Znanny(1:T,:)';
Cz        = Cz(:,:,1:T);

Nz        = size(Z,1);
Ny        = size(Cz,2);
Nhorizons = Ny - 1;


%% set missing obs to zero
Z(Znanny) = 0;

%% priors
Iy          = eye(Ny);
% SigmaY      = Iy;
decayY      = .95^2;
varY        = nanvar(Z(1,:)) * (1 - decayY); %#ok<NANVAR>
SigmaY      = varY * diag([.5, decayY.^(1:Ny-1)]); % set prior for nowcast error to half the nowcast update
SigmaDof    = Ny + 2;
SigmaT0     = SigmaY * (SigmaDof - Ny - 1); % if doDiffusePrior is true, this prior only needed to init missing values
cholSigmaT0 = chol(SigmaT0, 'lower'); % = sqrt(SigmaT) at present


%% prepare State Space: Y transition
F          = diag(ones(Ny-1,1),1);
F(end,end) = 1;
fff        = repmat(F, [1 1 T]);

%% prepare State Space: eta and eta-tilde transition

% H          = [Igap ones(Ngap, 1); zeros(1,Ngap) 1];
% invH       = eye(Ny) / H;

G0         = zeros(Ny); % dummy values
GPREV      = G0; % init
Ng         = numel(G0);

% minnesota-style prior (assuming equal residuals variances)
if isempty(G0var)
    G0var  = NaN(Ng,1); % diagonal elements
    %     Gtheta = [0.1 0.5];
    Gtheta = [0.05 0.5];
    ndx = 0;
    for i = 1 : Ny
        for j = 1 : Ny
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



%% prepare State Space: initial levels

Y0          = zeros(Ny, 1);
cholSigmaY0 = [20 * ones(Ny, 1), 5 * Iy];  % correlated initial levels
cholSigmaY0 = chol(cholSigmaY0 * cholSigmaY0', 'lower');

invcholSigmaY0    = Iy / cholSigmaY0;
invcholSigmaY0vec = invcholSigmaY0(:);

eta00          = zeros(Ny, 1);
cholSigmaEta00 = 5 * Iy;

invcholSigmaEta00    = Iy / cholSigmaEta00;
invcholSigmaEta00vec = invcholSigmaEta00(:);

%% prepare precision-based sampler
% parameters
NzT   = Nz * T;
NyT   = Ny * T;
Ny2   = Ny * 2;
% NyTp1 = NyT + Ny;
NyTp2 = NyT + Ny + Ny;
YY00  = sparse(1:Ny2, 1, [eta00; Y0], NyTp2, 1);


ndx1NyTp2 = (1 : NyTp2)';
onesNyTp2 = ones(NyTp2,1);


% AAy for AAy * Y = Y0 + ETA
arowsF    = repmat((1 : Ny)', 1 , Ny, T);
arowsF    = arowsF + permute(Ny * (1 : T), [1 3 2]);
arowsF    = arowsF + Ny; % to account for eta00
arows     = [ndx1NyTp2; arowsF(:)];

acolsF    = repmat(1 : NyT, Ny, 1);
acolsF    = acolsF + Ny; % to account for eta00
acols     = [ndx1NyTp2; acolsF(:)];

avalues   = [onesNyTp2; -fff(:)];

AAy       = sparse(arows, acols, avalues);


% prepare indices for AA
% AA for AA * Y = Y0 + RESID
arowsLag1     = repmat((1 : Ny)', 1 , Ny, T);
arowsLag1     = arowsLag1 + permute(Ny * (0 : T-1), [1 3 2]);
arowsLag1     = arowsLag1 + Ny2; % to account for eta00 and Y0

% arowsLag2 identical to arowsLag1
arows      = [ndx1NyTp2; arowsLag1(:); arowsLag1(:)];

acolsLag2  = repmat(1 : NyT, Ny, 1);
acolsLag1  = Ny + acolsLag2;
acols      = [ndx1NyTp2; acolsLag1(:); acolsLag2(:)];


%% BB
brows1  = repmat((1 : Ny)', 1 , Ny); % eta00
brows2  = Ny + brows1; % Y0
brows3  = repmat((1 : Ny)', 1 , Ny, T);
brows3  = brows3 + permute(Ny * (0 : T-1), [1 3 2]);
brows3  = Ny2 + brows3;
brows   = [brows1(:); brows2(:); brows3(:)];

bcols1  = repmat(1 : Ny, Ny , 1);
bcols2  = Ny  + bcols1;
bcols3  = Ny2 + repmat(1 : NyT, Ny, 1);
bcols   = [bcols1(:); bcols2(:); bcols3(:)];

% CC
crows     = repmat((1 : Nz)', 1 , Ny, T);
crows     = crows + permute(Nz * (0 : T-1), [1 3 2]);
ccols     = repmat(Ny2 + (1 : NyT), Nz, 1);
CC        = sparse(crows(:), ccols(:), Cz(:), NzT, NyTp2);

% ZZ
ZZ        = Z(:);
ZZnanny   = Znanny(:);

ZZ        = ZZ(~ZZnanny);
CC        = CC(~ZZnanny,:);


% one-off computations
[N1, N2]  = size(CC);
N2        = N2 - N1;

[QQ,RR]   = qr(CC');
QQ1       = QQ(:,1:N1)';
QQ2       = QQ(:,N1+1:end)';
RR1       = RR(1:N1,1:N1)';
% if spamaxabs(CC -  RR1 * QQ1) > 1e-8
%     warning('qr off')
% end

%% allocate memory for MCMC

Ydraws    = NaN(Ny, T, MCMCdraws);
YhatRB    = NaN(Nhorizons, MCMCdraws);
ETAdraws  = NaN(Ny, T, MCMCdraws);

Gdraws    = NaN(Ny, Ny, MCMCdraws);

sqrtSIGMA = NaN(Ny,Ny,MCMCdraws);

N2draws   = randn(rndStream, N2, 2 * MCMCdraws + 1);

ndxNsTT    = NyT - Ny + (1:Ny2);
StateTmean = NaN(Ny2, MCMCdraws);
StateTvar  = NaN(Ny2, Ny2, MCMCdraws);


%% draw initial values
sqrtsigmaPREV = iwishcholdraw(cholSigmaT0, SigmaDof, 1, rndStream);


%% draw blocks of random variables
z4vcvdraws = randn(rndStream, Ny, SigmaDof + T, 2 * MCMCdraws + 1);
z4gdraws   = randn(rndStream, Ng, 2 * MCMCdraws + 1);


%% loop over MCMC steps
if showProgress
    progressbar(0)
end

for mm = -MCMCdraws : MCMCdraws

    iterMCMC = mm + MCMCdraws + 1;

    %% precision-based state space sampler

    % update transition matrix
    aaLag1        = repmat(GPREV + F, [1 1 T]);
    aaLag2        = -repmat(GPREV * F, [1 1 T]);
    %  t = 1
    aaLag1(:,:,1) = F;
    aaLag2(:,:,1) = GPREV;
    avalues    = [onesNyTp2; -aaLag1(:); -aaLag2(:)];
    AA         = sparse(arows, acols, avalues, NyTp2, NyTp2);

    EY        = AA \ YY00;
    EZ        = CC * EY;
    Y1tilde   = RR1 \ (ZZ - EZ);

    if any(isnan(Y1tilde))
        error('NaN in Y1tilde')
    end

    % exploit scale SV structure
    invsqrtsigmaPREV = Iy / sqrtsigmaPREV;
    invbbb           = repmat(invsqrtsigmaPREV, [1 1 T]);
    invBB            = sparse(brows, bcols, [invcholSigmaEta00vec; invcholSigmaY0vec; invbbb(:)]);

    AAtildeQQ   = invBB * AA * QQ;
    AAtildeQQ2  = AAtildeQQ(:,N1+1:end);
    QQ2invSIG   = AAtildeQQ2' * AAtildeQQ;
    invQSIG21   = QQ2invSIG(:,1:N1);
    invQSIG22   = QQ2invSIG(:,1+N1:end);

    cholinvQSIG22 = chol(invQSIG22, 'lower');

    y2hat         = - cholinvQSIG22 \ (invQSIG21 * Y1tilde);
   
    % compute y2hat and y2draw separately to store mean
    y2hat         = cholinvQSIG22' \ y2hat;
    y2draw        = cholinvQSIG22' \ N2draws(:,mm + MCMCdraws + 1);
    
    YHAT          = EY + QQ1' * Y1tilde + QQ2' * y2hat;
    YDRAW         = YHAT + QQ2' * y2draw;

    thisShat      = YHAT(ndxNsTT);

    statevol      = cholinvQSIG22 \ QQ2(:,ndxNsTT); % note: statevol is not square!
    thisSvar      = statevol' * statevol;

    thisY         = reshape(YDRAW(Ny2+1:end), Ny, T);

    ETADRAW       = AAy * YDRAW - YY00;
    eta           = reshape(ETADRAW(Ny2+1:end), Ny, T);
    etaLag        = cat(2, ETADRAW(1:Ny), eta(:,1:end-1));



    %% draw VAR for eta
    iSigmaResid = invsqrtsigmaPREV' * invsqrtsigmaPREV;
    checkdiff(iSigmaResid, inv(sqrtsigmaPREV * sqrtsigmaPREV'));
    [GPREV, resid]  = bayesVectorRegressionGibbsDraw1(eta', etaLag', ...
        iSigmaResid, invG0varG0, invG0var,  z4gdraws(:, iterMCMC));

    GPREV           = GPREV';
    maxabsG = max(abs(eig(GPREV)));
    if  maxabsG > 1
        warning('G unstable: max root is %6.2f (T=%d)', maxabsG, T)
    end

    %% draw constant VCV
    sqrtsigmaPREV = bayesSQRTVCVgibbsDraw1(SigmaT0, SigmaDof, resid, z4vcvdraws(:,:,iterMCMC));

    %% store draws post burnin
    if mm > 0
        Ydraws(:,:,mm)    = thisY;
        ETAdraws(:,:,mm)  = eta;
        sqrtSIGMA(:,:,mm) = sqrtsigmaPREV;

        Gdraws(:,:,mm)    = GPREV;
    
        StateTmean(:,mm)  = thisShat;
        StateTvar(:,:,mm) = thisSvar;
    end

    if showProgress
        progressbar((mm + MCMCdraws + 1) / (2 * MCMCdraws + 1))
    end
end

%% simulate draws from predictive density

% allocate memory (note: will all be reshaped and permuted as needed below)
YdensityDraws = NaN(Nfedraws, Nhorizons, MCMCdraws);
YstateTp1mean = NaN(Ny, MCMCdraws);
YstateTp1var  = NaN(Ny, Ny, MCMCdraws);
yTpHmean      = NaN(1, Nhorizons, MCMCdraws);
yTpHvar       = NaN(1, Nhorizons, MCMCdraws);


for mm = 1 : MCMCdraws

    % collect MCMC parameter draws
    B             = sqrtSIGMA(:,:,mm);
    BB            = B * B';
    thisY         = Ydraws(:,end,mm);
    thisETA       = ETAdraws(:,end,mm);
    
    % collect one-step ahead moments for state vector
    % build VAR2 companion matrix (and recall that order is Y(t-1), Y(t)
    Scompanion                 = zeros(Ny2);
    Scompanion(1:Ny,Ny+(1:Ny)) = Iy;
    Scompanion(Ny+(1:Ny),:)    = [- F * GPREV, F + GPREV];
    Ycompanion                 = Scompanion(Ny+(1:Ny),:);
    % pick only loadings on Y(t)
    YstateTp1mean(:,mm)    = Ycompanion * StateTmean(:,mm);
    YvarTT                 = Ycompanion * StateTvar(:,:,mm) * Ycompanion';
    YstateTp1var(:,:,mm)   = YvarTT + BB;

    % collect TplusH moments of outcome variable
    Smean  = StateTmean(:,mm);
    Svar   = StateTvar(:,:,mm);
    thisBBB                    = zeros(Ny2);
    thisBBB(Ny+1:Ny2,Ny+1:Ny2) = BB;
    for jj = 1 : Nhorizons
        Smean             = Scompanion * Smean;
        Svar              = Scompanion * Svar * Scompanion' + thisBBB;
        yTpHmean(:,jj,mm) = Smean(Ny+1);
        yTpHvar(1,jj,mm)  = Svar(Ny+1,Ny+1);
    end

    % construct shocks to eta
    z      = randn(rndStream,Ny,Nhorizons * Nfedraws);
    shocks = B * z;
    shocks = reshape(shocks, Ny, Nhorizons, Nfedraws);
    shocks = permute(shocks, [1 3 2]);

    % simulate density draws
    for jj = 1 : Nhorizons
        thisETA                = GPREV * thisETA + shocks(:,:,jj);
        thisY                  = F * thisY + thisETA;
        YdensityDraws(:,jj,mm) = thisY(1,:);
    end

    % simulate RB means
    thisY         = Ydraws(:,end,mm);
    thisETA       = ETAdraws(:,end,mm);
    for jj = 1 : Nhorizons
        thisETA                = GPREV * thisETA;
        thisY                  = F * thisY + thisETA;
        YhatRB(jj,mm)          = thisY(1);
    end

end


Ndraws = Nfedraws * MCMCdraws;
YdensityDraws = permute(YdensityDraws, [2 1 3]);
YdensityDraws = reshape(YdensityDraws, Nhorizons, Ndraws);

yTpHmean  = permute(yTpHmean, [3 2 1]);
yTpHvar   = permute(yTpHvar, [3 2 1]);

