function [YdensityDraws, Ydraws, YhatRB, ETAdraws, ...
    Gdraws, ...
    SVSCALEdraws, hVAR, RHO, sqrtSIGMA, mustarVAR, ...
    YstateTp1mean, YstateTp1var, yTpHmean, yTpHvar  ...
    ] = mcmcsamplerVARSTATEtrendgapSV(Z, Znanny, Cz, ...
    G0var, ...
    T, MCMCdraws, Nfedraws, rndStream, showProgress)

%#ok<*UNRCH>
if nargin < 8
    showProgress = false;
end


%% cut sample, set parameters
maxShake  = 1e3;

Z         = Z(1:T,:)'; % note transpose
Znanny    = Znanny(1:T,:)';
Cz        = Cz(:,:,1:T);

Nz        = size(Z,1);
Ny        = size(Cz,2);
Ngap      = Ny - 1;
Nhorizons = Ny - 1;

[KSC, KSCt, logy2offset] = getKSC10values(T, Ngap); % using 10-point grid as in Omori et al (2007, JoE)
% note: single SV, but Ngap signals

%% set missing obs to zero 
Z(Znanny) = 0;

%% priors
Iy          = eye(Ny);
Igap        = eye(Ngap);
SigmaYgap   = Igap;
SigmaDof    = Ngap + 2;
SigmaT0     = SigmaYgap * (SigmaDof - Ngap - 1); 
cholSigmaT0 = chol(SigmaT0, 'lower'); % = sqrt(SigmaT) at present

% scale SV
% Eh0  = 0;
% Vh0  = 0;
hvarDof   = 3;
hvarT     = 0.2^2 * (hvarDof - 1 - 1);

rho0  = 0.8;
rhoV0 = 0.2^2;

mustarvarDof = 3;
mustarvarT   = (.1 / 4)^2 * (mustarvarDof - 1 - 1); 

%% prepare State Space: Y transition
F          = diag(ones(Ny-1,1),1);
F(end,end) = 1;
fff        = repmat(F, [1 1 T]);

%% prepare State Space: eta and eta-tilde transition

% H          = [Igap ones(Ngap, 1); zeros(1,Ngap) 1];
% invH       = eye(Ny) / H;

G0         = zeros(Ngap); % dummy values
GPREV      = G0; % init todo
Ng         = numel(G0);

% minnesota-style prior (assuming equal residuals variances)
if isempty(G0var)
    G0var  = NaN(Ng,1); % diagonal elements
    %     Gtheta = [0.1 0.5];
    Gtheta = [0.05 0.5];
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

Gdraws    = NaN(Ngap, Ngap, MCMCdraws);

hh        = NaN([1 T MCMCdraws]);
hVAR      = NaN(1, MCMCdraws);
RHO       = NaN(1, MCMCdraws);
sqrtSIGMA = NaN(Ngap,Ngap,MCMCdraws);
mustarVAR = NaN(1, MCMCdraws);

N2draws   = randn(rndStream, N2, 2 * MCMCdraws + 1);

ndxNsTT    = NyT - Ny + (1:Ny2);
StateTmean = NaN(Ny2, MCMCdraws);
StateTvar  = NaN(Ny2, Ny2, MCMCdraws);

GGPREV = zeros(Ny);

%% draw initial values
sqrtsigmaPREV = iwishcholdraw(cholSigmaT0, SigmaDof, 1, rndStream);
hvarPREV      = igamdraw(hvarT, hvarDof, 1);

if mustarvarDof == 0
    mustarvolPREV = 1;
else
    mustarvolPREV = sqrt(igamdraw(mustarvarT, mustarvarDof, 1));
end

hPREV         = sqrt(hvarPREV) .* randn(rndStream, 1, T);
hPREV         = cumsum(hPREV); % h0 = 0

rhoPREV = rand(rndStream,1) * 2 - 1;


%% draw blocks of random variables 
z4vcvdraws = randn(rndStream, Ngap, SigmaDof + T, 2 * MCMCdraws + 1);
z4gdraws   = randn(rndStream, Ng, 2 * MCMCdraws + 1);


%% loop over MCMC steps
if showProgress
    progressbar(0)
end

for mm = -MCMCdraws : MCMCdraws

    iterMCMC = mm + MCMCdraws + 1;

    SVscale = exp(hPREV * 0.5);

    %% precision-based state space sampler

    %     GGPREV        = H * blkdiag(GPREV, 0) * invH;
    %     GGPREV        = [GPREV, -sum(GPREV,2); zeros(1,Ny)];
    GGPREV(1:Ngap,:)    = [GPREV, -sum(GPREV,2)];

    % update transition matrix
    aaLag1        = repmat(GGPREV + F, [1 1 T]);
    aaLag2        = -repmat(GGPREV * F, [1 1 T]);
    %  t = 1
    aaLag1(:,:,1) = F;
    aaLag2(:,:,1) = GGPREV;
    avalues    = [onesNyTp2; -aaLag1(:); -aaLag2(:)];
    AA         = sparse(arows, acols, avalues, NyTp2, NyTp2);

    EY        = AA \ YY00;
    EZ        = CC * EY;
    Y1tilde   = RR1 \ (ZZ - EZ);

    if any(isnan(Y1tilde))
        error('NaN in Y1tilde')
    end

    % exploit scale SV structure
    invsqrtsigmaPREV = Igap / sqrtsigmaPREV;
    invbgapsv        = invsqrtsigmaPREV ./ permute(SVscale, [1 3 2]);
    
    invbbsv                        = zeros(Ny, Ny, T);
    invbbsv(1:Ngap, 1:Ngap, :)     = invbgapsv;
    invbbsv(1:Ngap, end, :)        = - sum(invbgapsv,2); % partitioned inverse computation
    invbbsv(end, end, :)           = 1 / mustarvolPREV;

    invBB            = sparse(brows, bcols, [invcholSigmaEta00vec; invcholSigmaY0vec; invbbsv(:)]);

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

    mustar        = eta(end,:);
    mustarLag     = etaLag(end,:);
    
    etaGAP        = eta(1:Ngap,:) - mustar;
    etaGAPlag     = etaLag(1:Ngap,:) - mustarLag;
    
    %     check = invH * eta;
    %     checkdiff(etaGAP, check(1:Ngap,:));
    %     checkdiff(mustar, check(end,:));
    %     checkdiff(eta, H * cat(1, etaGAP, mustar));

    %% draw VAR for eta
    GLAST        = GPREV; % store in case AR1 sampling keeps hittting unit root

    etaGAPscaled    = (etaGAP ./ SVscale)';
    SVlag           = [1, SVscale(1:end-1)]; % recall that initial SVscale is fixed at 1
    etaGAPscaledLag = (etaGAPlag ./ SVlag)';

    iSigmaResid = invsqrtsigmaPREV' * invsqrtsigmaPREV;
    % checkdiff(iSigmaResid, inv(sqrtsigmaPREV * sqrtsigmaPREV'));
    [GPREV, residGAPscaled]  = bayesVectorRegressionGibbsDraw1(etaGAPscaled, etaGAPscaledLag, ...
        iSigmaResid, invG0varG0, invG0var,  z4gdraws(:, iterMCMC));

    residGAP        = residGAPscaled' .* SVscale;
    GPREV           = GPREV';
    maxabsG = max(abs(eig(GPREV)));
    if  maxabsG > 1
        warning('G unstable: max root is %6.2f (T=%d)', maxabsG, T)
    end

    %% draw constant VCV
    sqrtsigmaLAST = sqrtsigmaPREV; % store in case AR1 sampling keeps hittting unit root


    %  residScaled obtained from previous step

    sqrtsigmaPREV = bayesSQRTVCVgibbsDraw1(SigmaT0, SigmaDof, residGAPscaled, z4vcvdraws(:,:,iterMCMC));

    %% common SV draw
    hLAST           = hPREV; % store in case AR1 sampling keeps hittting unit root

    residTilde      = sqrtsigmaPREV \ residGAP;
    logy2           = log(residTilde.^2 + logy2offset);
    hPREV           = CommonStochVolAR1(logy2, hPREV, rhoPREV, sqrt(hvarPREV), KSC, KSCt, Ngap, T, rndStream);


    %% SV-AR1 rho draw
    rhoLAST  = rhoPREV; % store in case AR1 sampling keeps hittting unit root
    hLAG     = [0, hPREV(1:end-1)];
    rhodraws = bayesAR1SURdraw(hPREV', hLAG', hvarPREV, rho0, rhoV0, maxShake, rndStream);
    % note on previous line: single equation, no SUR needed ...
    % but the SUR function conveniently draws multiple rhos from posterior ....
    shake    = 0;
    OK       = false;
    while ~OK && shake < maxShake
        shake = shake + 1;
        OK    = abs(rhodraws(shake)) < 1;
    end

    if OK
        rhoPREV = rhodraws(:,shake);
    else
        warning('maxShake exhausted T=%d, mm=%d', T, mm)
        rhoPREV       = rhoLAST;
        hPREV         = hLAST;
        sqrtsigmaPREV = sqrtsigmaLAST;
        GPREV         = GLAST;
        % keep the old values for this step; not a problem if it occurs rarely / during burnin
		% need to watch output
    end

    %% draw hvar
    hshock   = hPREV - rhoPREV * hLAG;
    %     hvarPREV = bayesVCVgibbsDraw1(hvarT, hvarDof, hshock', rndStream, true);
    hvarPREV = igamVarianceDraw(hshock, hvarT, hvarDof);

    %% draw mustarvar
    %     mustarvolPREV = sqrt(bayesVCVgibbsDraw1(mustarvarT, mustarvarDof, mustar', rndStream, true));
    mustarvolPREV = sqrt(igamVarianceDraw(mustar, mustarvarT, mustarvarDof));


    %% store draws post burnin
    if mm > 0
        Ydraws(:,:,mm)    = thisY;
        ETAdraws(:,:,mm)  = eta;
        hh(:,:,mm)        = hPREV;
        hVAR(:,mm)        = hvarPREV;
        RHO(:,mm)         = rhoPREV;
        mustarVAR(:,mm)   = mustarvolPREV; % will be converted later into variances
        sqrtSIGMA(:,:,mm) = sqrtsigmaPREV;

        Gdraws(:,:,mm)    = GPREV;


        StateTmean(:,mm)  = thisShat;
        StateTvar(:,:,mm) = thisSvar;
    end

    if showProgress
        progressbar((mm + MCMCdraws + 1) / (2 * MCMCdraws + 1))
    end
end

%% collect SV


SVSCALEdraws = squeeze(exp(hh * 0.5));

%% simulate draws from predictive density

% allocate memory (note: will all be reshaped and permuted as needed below)
YdensityDraws = NaN(Nfedraws, Nhorizons, MCMCdraws);
YstateTp1mean = NaN(Ny, Nfedraws, MCMCdraws); 
YstateTp1var  = NaN(Ny, Ny, Nfedraws, MCMCdraws); 
yTpHmean      = NaN(Nfedraws, Nhorizons, MCMCdraws); 
yTpHvar       = NaN(Nfedraws, Nhorizons, MCMCdraws); 

% B   = zeros(Ny);
BB  = zeros(Ny);

for mm = 1 : MCMCdraws

    % collect MCMC parameter draws
    rho           = RHO(mm);
    hT            = hh(1,end,mm);
    hvol          = sqrt(hVAR(mm));
    mustarvol     = mustarVAR(mm);
    Bgap          = sqrtSIGMA(:,:,mm);
    GPREV         = Gdraws(:,:,mm);
    %     GGPREV        = H * blkdiag(GPREV, 0) * invH;
    GGPREV(1:Ngap,:)    = [GPREV, -sum(GPREV,2)];
    thisETA       = ETAdraws(:,end,mm);
    thisY         = Ydraws(:,end,mm);

    % generate SV (note the second dimension are the simulated forecasts
    hshocks = randn(rndStream,1,Nhorizons,Nfedraws);
    for m = 1 : Nfedraws
        hshocks(1,:,m) = hvol * hshocks(1,:,m);
    end

    h = NaN(1, Nhorizons,Nfedraws);
    j = 1;
    h(1,j,:) = rho * hT + hshocks(1,j,:);
    for j = 2 : Nhorizons
        h(1,j,:) = rho .* h(1,j-1,:) + hshocks(1,j,:);
    end
    SVscale = exp(h * .5);

    % collect one-step ahead moments for state vector
    % build VAR2 companion matrix (and recall that order is Y(t-1), Y(t)
    Scompanion                 = zeros(Ny2);
    Scompanion(1:Ny,Ny+(1:Ny)) = Iy;
    Scompanion(Ny+(1:Ny),:)    = [- F * GGPREV, F + GGPREV]; 
    Ycompanion                 = Scompanion(Ny+(1:Ny),:);
    % pick only loadings on Y(t)
    YstateTp1mean(:,:,mm)      = repmat(Ycompanion * StateTmean(:,mm), 1, Nfedraws);
    YvarTT                 = Ycompanion * StateTvar(:,:,mm) * Ycompanion';

    BB(:)                  = mustarvol^2;
    BBgap                  = Bgap * Bgap';
    varscale               = (squeeze(SVscale(1,1,:))).^2;

    for ii = 1 : Nfedraws
        thisBB                  = BB;
        thisBB(1:Ngap,1:Ngap)   = thisBB(1:Ngap,1:Ngap) + varscale(ii) * BBgap;

        YstateTp1var(:,:,ii,mm) = YvarTT + thisBB;
    end

    % collect TplusH moments of outcome variable
    Smean  = StateTmean(:,mm);
    Svar   = repmat(StateTvar(:,:,mm), 1, 1, Nfedraws);
    BB(:)                  = mustarvol^2;
    BBgap                  = Bgap * Bgap';
    thisBBB                = zeros(Ny2);
    for jj = 1 : Nhorizons
        Smean             = Scompanion * Smean;
        yTpHmean(:,jj,mm) = Smean(Ny+1);

        varscale          = squeeze(SVscale(1,jj,:).^2);
        for ii = 1 : Nfedraws
            Svar(:,:,ii)               = Scompanion * Svar(:,:,ii) * Scompanion';
            thisBB                     = BB;
            thisBB(1:Ngap,1:Ngap)      = thisBB(1:Ngap,1:Ngap) + varscale(ii) * BBgap;
            thisBBB(Ny+1:Ny2,Ny+1:Ny2) = thisBB;
            Svar(:,:,ii)               = Svar(:,:,ii) + thisBBB;
            yTpHvar(ii,jj,mm)          = Svar(Ny+1,Ny+1,ii);
        end
    end

    % construct eta
    B        = zeros(Ny);
    B(:,Ny)  = mustarvol;
    shocks   = NaN(Ny,Nhorizons,Nfedraws);
    z        = randn(rndStream,Ny,Nhorizons,Nfedraws); % maintain same index order as in baseline model
    for ii = 1 : Nfedraws
        for jj = 1 : Nhorizons
            B(1:Ngap,1:Ngap) = SVscale(1,jj,ii) * Bgap;
            shocks(:,jj,ii)  = B * z(:,jj,ii);
        end
    end
    shocks = permute(shocks, [1 3 2]);

    % simulate density draws
    for jj = 1 : Nhorizons
        thisETA                = GGPREV * thisETA + shocks(:,:,jj);
        thisY                  = F * thisY + thisETA;
        YdensityDraws(:,jj,mm) = thisY(1,:);
    end

    % simulate RB means
    thisY         = Ydraws(:,end,mm);
    thisETA       = ETAdraws(:,end,mm);
    for jj = 1 : Nhorizons
        thisETA                = GGPREV * thisETA;
        thisY                  = F * thisY + thisETA;
        YhatRB(jj,mm)          = thisY(1);
    end

end


Ndraws = Nfedraws * MCMCdraws;
YdensityDraws = permute(YdensityDraws, [2 1 3]);
YdensityDraws = reshape(YdensityDraws, Nhorizons, Ndraws);

yTpHmean  = permute(yTpHmean, [1 3 2]);
yTpHmean  = reshape(yTpHmean, [Ndraws, Nhorizons]);

yTpHvar  = permute(yTpHvar, [1 3 2]);
yTpHvar  = reshape(yTpHvar, [Ndraws, Nhorizons]);

YstateTp1mean = reshape(YstateTp1mean, Ny, Nfedraws * MCMCdraws);
YstateTp1var  = reshape(YstateTp1var, Ny, Ny, Nfedraws * MCMCdraws);

% convert mustarVAR from vols to variances
mustarVAR = mustarVAR.^2;

