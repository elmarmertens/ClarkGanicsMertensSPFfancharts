function [YdensityDraws, Ydraws, ETAdraws, ETASVdraws, ...
    SVSCALEdraws, hVAR, RHO, sqrtSIGMA, trendnoiseVAR, mustarVAR, ...
    YstateTp1mean, YstateTp1var, yTpHmean, yTpHvar  ...
    ] = mcmcsamplerSTATEtrendgaplongSV(Ngap, Z, Znanny, Cz, ...
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
if isempty(Ngap)
    Ngap      = Ny - 1;
end
Ntrendnoise    = Ny - Ngap - 1;

Nhorizons = Ny - 1;

[KSC, KSCt, logy2offset] = getKSC10values(T, Ngap); % using 10-point grid as in Omori et al (2007, JoE)
% note: single SV, but Ngap signals

%% set missing obs to zero 
Z(Znanny) = 0;

%% priors
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

trendnoisevarDof = 3;
trendnoisevarT   = (.1 / 4)^2 * (trendnoisevarDof - 1 - 1); 

%% prepare State Space
Iy         = eye(Ny);

F          = diag(ones(Ny-1,1),1);
F(end,end) = 1;
FF         = repmat(F, [1 1 T]);

Y0          = zeros(Ny, 1);
cholSigmaY0 = [20 * ones(Ny, 1), 5 * Iy];  % correlated initial levels
cholSigmaY0 = chol(cholSigmaY0 * cholSigmaY0', 'lower');

invcholSigmaY0    = Iy / cholSigmaY0;
invcholSigmaY0vec = invcholSigmaY0(:);


%% prepare precision-based sampler
% parameters
NzT   = Nz * T;
NyT   = Ny * T;
NyTp1 = NyT + Ny;
YY00  = sparse(1:Ny, 1, Y0, NyTp1, 1);


% AA
ndx1NyTp1 = (1 : NyTp1)';
arows     = repmat((1 : Ny)', 1 , Ny, T);
arows     = arows + permute(Ny * (1 : T), [1 3 2]);
rowndx    = [ndx1NyTp1; arows(:)];
colndx    = [ndx1NyTp1; vec(repmat(1 : NyT, Ny, 1))];
values    = [ones(NyTp1,1); -FF(:)];
AA        = sparse(rowndx, colndx, values);

% BB
brows1  = repmat((1 : Ny)', 1 , Ny);
brows2  = repmat((1 : Ny)', 1 , Ny, T);
brows2  = brows2 + permute(Ny * (0 : T-1), [1 3 2]);
brows2  = Ny + brows2;
brows   = [brows1(:); brows2(:)];

bcols1  = repmat(1 : Ny, Ny , 1);
bcols2  = Ny + repmat(1 : NyT, Ny, 1);
bcols   = [bcols1(:); bcols2(:)];

% CC
crows     = repmat((1 : Nz)', 1 , Ny, T);
crows     = crows + permute(Nz * (0 : T-1), [1 3 2]);
ccols     = repmat(Ny + (1 : NyT), Nz, 1);
CC        = sparse(crows(:), ccols(:), Cz(:), NzT, NyTp1);

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
EY        = AA \ YY00;
EZ        = CC * EY;
Y1tilde   = RR1 \ (ZZ - EZ);

if any(isnan(Y1tilde))
    error('NaN in Y1tilde')
end

%% allocate memory for MCMC

Ydraws      = NaN(Ny, T, MCMCdraws);
ETAdraws    = NaN(Ny, T, MCMCdraws);
ETASVdraws  = NaN(Ny, T, MCMCdraws);

hh        = NaN([1 T MCMCdraws]);
hVAR      = NaN(1, MCMCdraws);
RHO       = NaN(1, MCMCdraws);
sqrtSIGMA = NaN(Ngap,Ngap,MCMCdraws);
mustarVAR = NaN(1, MCMCdraws);
trendnoiseVAR  = NaN(Ntrendnoise, MCMCdraws);

N2draws   = randn(rndStream, N2, 2 * MCMCdraws + 1);

ndxNyTT     = NyT + (1:Ny);
YstateTmean = NaN(Ny, MCMCdraws);
YstateTvar  = NaN(Ny, Ny, MCMCdraws);

%% draw initial values
sqrtsigmaPREV = iwishcholdraw(cholSigmaT0, SigmaDof, 1, rndStream);
hvarPREV      = igamdraw(hvarT, hvarDof, 1);

if trendnoisevarDof == 0
    trendnoisevolPREV  = ones(Ntrendnoise,1);
else
    trendnoisevolPREV  = sqrt(igamdraw(trendnoisevarT, trendnoisevarDof, Ntrendnoise, 1));
end
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

%% loop over MCMC steps
if showProgress
    progressbar(0)
end

for mm = -MCMCdraws : MCMCdraws

    iterMCMC = mm + MCMCdraws + 1;

    SVscale = exp(hPREV * 0.5);

    %% precision-based state space sampler

    % exploit scale SV structure
    invsqrtsigmaPREV = Igap / sqrtsigmaPREV;
    invbgapsv        = invsqrtsigmaPREV ./ permute(SVscale, [1 3 2]);

    
    invbbsv                                 = zeros(Ny, Ny, T);
    invbbsv(1:Ngap, 1:Ngap, :)              = invbgapsv;
    invbbsv(1:Ngap, end, :)                 = - sum(invbgapsv,2); % partitioned inverse computation
    
    if Ntrendnoise > 0
        invtrendnoisevolPREV                = 1 ./ trendnoisevolPREV;
        %         % repmat -- SLOWER
        %         invbbsv(Ngap+(1:Ntrendnoise), Ngap+(1:Ntrendnoise), :)   = ...
        %             repmat(diag(invtrendnoisevolPREV), [1 1 T]);
        %         invbbsv(Ngap+(1:Ntrendnoise), end, :) = - repmat(invtrendnoisevolPREV, [1 1 T]);
        % loop
        for jj = 1 : Ntrendnoise
            invbbsv(Ngap+jj, Ngap+jj,:) = invtrendnoisevolPREV(jj);
            invbbsv(Ngap+jj, end,:)     = -invtrendnoisevolPREV(jj);
        end
    end

    invbbsv(end, end, :)                  = 1 ./ mustarvolPREV;
    
    invBB            = sparse(brows, bcols, [invcholSigmaY0vec; invbbsv(:)]);

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

    thisYhat      = YHAT(ndxNyTT);

    yvol          = cholinvQSIG22 \ QQ2(:,ndxNyTT); % note: yvol is not square!
    thisYvar      = yvol' * yvol;

    thisY         = reshape(YDRAW(Ny+1:end), Ny, T);

    ETADRAW       = AA * YDRAW - YY00;
    eta           = reshape(ETADRAW(Ny+1:end), Ny, T);

    mustar        = eta(end,:);
    etaGAP        = eta(1:Ngap,:) - mustar;
    etaTRENDNOISE = transpose(eta(Ngap+(1:Ntrendnoise),:) - mustar);

    %% draw constant VCV
    sqrtsigmaLAST = sqrtsigmaPREV; % store in case AR1 sampling keeps hittting unit root


    etaScaled     = transpose(etaGAP ./ SVscale); % could also use SVscale with automatic matrix expansion

    sqrtsigmaPREV = bayesSQRTVCVgibbsDraw1(SigmaT0, SigmaDof, etaScaled, z4vcvdraws(:,:,iterMCMC));

    %% common SV draw
    hLAST           = hPREV; % store in case AR1 sampling keeps hittting unit root

    etaTilde        = sqrtsigmaPREV \ etaGAP;
    logy2           = log(etaTilde.^2 + logy2offset);
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
        % keep the old values for this step; not a problem if it occurs rarely / during burnin
		% need to watch output
    end

    %% draw hvar
    hshock   = hPREV - rhoPREV * hLAG;
    %     hvarPREV = bayesVCVgibbsDraw1(hvarT, hvarDof, hshock', rndStream, true);
    hvarPREV = igamVarianceDraw(hshock, hvarT, hvarDof);

    %% draw trendnoisevar
    if Ntrendnoise > 0 % ~isempty(etaTRENDNOISE)
        trendnoisevolPREV = sqrt(igamVarianceDraw(etaTRENDNOISE, trendnoisevarT, trendnoisevarDof));
    end

    %% draw mustarvar
    %     mustarvolPREV = sqrt(bayesVCVgibbsDraw1(mustarvarT, mustarvarDof, mustar', rndStream, true));
    mustarvolPREV = sqrt(igamVarianceDraw(mustar, mustarvarT, mustarvarDof));



    %% store draws post burnin
    if mm > 0
        Ydraws(:,:,mm)       = thisY;
        ETAdraws(:,:,mm)     = eta;
        hh(:,:,mm)           = hPREV;
        hVAR(:,mm)           = hvarPREV;
        RHO(:,mm)            = rhoPREV;
        trendnoiseVAR(:,mm)  = trendnoisevolPREV; % will be converted later into variances
        mustarVAR(:,mm)      = mustarvolPREV; % will be converted later into variances
        sqrtSIGMA(:,:,mm)    = sqrtsigmaPREV;


        YstateTmean(:,mm)  = thisYhat;
        YstateTvar(:,:,mm) = thisYvar;
    end

    if showProgress
        progressbar((mm + MCMCdraws + 1) / (2 * MCMCdraws + 1))
    end
end

%% collect SV


SVSCALEdraws = squeeze(exp(hh * 0.5));


for mm = 1 : MCMCdraws
    thisB                  = zeros(Ny,Ny);
    thisB(:,Ny)            = mustarVAR(mm); % recall: mustarVAR stores VOL (will be converted to VAR only later)
    thisB(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) = diag(trendnoiseVAR(:,mm)); % recall: trednoiseVAR stores VOL (will be converted to VAR only later)
    for jj = 1 : T
        thisB(1:Ngap,1:Ngap)   = SVSCALEdraws(jj,mm) * sqrtSIGMA(:,:,mm);
        ETASVdraws(:,jj,mm)    = sqrt(diag(thisB * thisB'));
    end
end

%% simulate draws from predictive density

% allocate memory (note: will all be reshaped and permuted as needed below)
YdensityDraws = NaN(Nfedraws, Nhorizons, MCMCdraws);
YstateTp1mean = NaN(Ny, Nfedraws, MCMCdraws); 
YstateTp1var  = NaN(Ny, Ny, Nfedraws, MCMCdraws); 
yTpHmean      = NaN(Nfedraws, Nhorizons, MCMCdraws); 
yTpHvar       = NaN(Nfedraws, Nhorizons, MCMCdraws); 

B   = zeros(Ny);
BB  = zeros(Ny);

for mm = 1 : MCMCdraws

    % collect MCMC parameter draws
    rho            = RHO(mm);
    hT             = hh(1,end,mm);
    hvol           = sqrt(hVAR(mm));
    mustarvol      = mustarVAR(mm);
    trendnoisevol  = trendnoiseVAR(:,mm);
    Bgap           = sqrtSIGMA(:,:,mm);
    thisY          = Ydraws(:,end,mm);

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
    YstateTp1mean(:,:,mm)  = repmat(F * YstateTmean(:,mm) , 1, Nfedraws);
    YvarTT                 = F * YstateTvar(:,:,mm) * F';

    BB(:)                  = mustarvol^2;
    BB(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) = ...
        BB(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) + diag(trendnoisevol.^2);
    BBgap                  = Bgap * Bgap';
    varscale               = (squeeze(SVscale(1,1,:))).^2;

    for ii = 1 : Nfedraws
        thisBB                  = BB;
        thisBB(1:Ngap,1:Ngap)   = thisBB(1:Ngap,1:Ngap) + varscale(ii) * BBgap;

        YstateTp1var(:,:,ii,mm) = YvarTT + thisBB;
    end

    % collect TplusH moments of outcome variable
    Ymean  = YstateTmean(:,mm);
    Yvar   = repmat(YstateTvar(:,:,mm), 1, 1, Nfedraws);
    BB(:)                  = mustarvol^2;
    BB(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) = ...
        BB(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) + diag(trendnoisevol.^2);
    BBgap                  = Bgap * Bgap';
    for jj = 1 : Nhorizons
        Ymean             = F * Ymean;
        yTpHmean(:,jj,mm) = Ymean(1);
        %         volscale          = squeeze(SVscale(1,jj,:));
        varscale          = squeeze(SVscale(1,jj,:).^2);
        for ii = 1 : Nfedraws
            Yvar(:,:,ii)         = F * Yvar(:,:,ii) * F';
            
            
            thisBB                  = BB;
            thisBB(1:Ngap,1:Ngap)   = thisBB(1:Ngap,1:Ngap) + varscale(ii) * BBgap;

            
            Yvar(:,:,ii)         = Yvar(:,:,ii) + thisBB;
            yTpHvar(ii,jj,mm)    = Yvar(1,1,ii);
        end
    end

    % construct eta
    B(:,Ny)  = mustarvol;
    B(Ngap+(1:Ntrendnoise),Ngap+(1:Ntrendnoise)) = diag(trendnoisevol);
    eta      = NaN(Ny,Nhorizons,Nfedraws);
    z        = randn(rndStream,Ny,Nhorizons,Nfedraws); % maintain same index order as in baseline model
    for ii = 1 : Nfedraws
        for jj = 1 : Nhorizons
            B(1:Ngap,1:Ngap) = SVscale(1,jj,ii) * Bgap;
            eta(:,jj,ii)     = B * z(:,jj,ii);
        end
    end
    eta = permute(eta, [1 3 2]);

    % simulate density draws
    for jj = 1 : Nhorizons
        thisY                  = F * thisY + eta(:,:,jj); % note: expansion into Ny x Nfedraws array at ii=1
        YdensityDraws(:,jj,mm) = thisY(1,:);
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
trendnoiseVAR  = trendnoiseVAR.^2;

