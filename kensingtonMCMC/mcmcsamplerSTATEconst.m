function [YdensityDraws, Ydraws, ETAdraws, ETASVdraws, ...
    sqrtSIGMA, ...
    YstateTp1mean, YstateTp1var, yTpHmean, yTpHvar  ...
    ] = mcmcsamplerSTATEconst(Z, Znanny, Cz, ...
    doDiffusePrior, ...
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

Iy = eye(Ny);
% constant VCV
zeroNy      = zeros(Ny); % only needed if doDiffusePrior equal to true

% SigmaY      = Iy;
decayY      = .95^2;
varY        = nanvar(Z(1,:)) * (1 - decayY); %#ok<NANVAR>
SigmaY      = varY * diag([.5, decayY.^(1:Ny-1)]); % set prior for nowcast error to half the nowcast update
SigmaDof    = Ny + 2;
SigmaT0     = SigmaY * (SigmaDof - Ny - 1); % if doDiffusePrior is true, this prior only needed to init missing values
cholSigmaT0 = chol(SigmaT0, 'lower'); % = sqrt(SigmaT) at present


%% prepare State Space
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

sqrtSIGMA = NaN(Ny,Ny,MCMCdraws);

N2draws   = randn(rndStream, N2, 2 * MCMCdraws + 1);

ndxNyTT     = NyT + (1:Ny);
YstateTmean = NaN(Ny, MCMCdraws);
YstateTvar  = NaN(Ny, Ny, MCMCdraws);

%% draw initial values
sqrtsigmaPREV = iwishcholdraw(cholSigmaT0, SigmaDof, 1, rndStream);


%% draw blocks of random variables 
z4vcvdraws = randn(rndStream, Ny, SigmaDof + T, 2 * MCMCdraws + 1);

%% loop over MCMC steps
if showProgress
    progressbar(0)
end

for mm = -MCMCdraws : MCMCdraws

    iterMCMC = mm + MCMCdraws + 1;


    %% precision-based state space sampler

    % exploit scale SV structure
    invsqrtsigmaPREV = Iy / sqrtsigmaPREV;
    invbbb           = repmat(invsqrtsigmaPREV, [1 1 T]);
    invBB            = sparse(brows, bcols, [invcholSigmaY0vec; invbbb(:)]);

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


    %% draw constant VCV
    if doDiffusePrior
        sqrtsigmaPREV = bayesSQRTVCVgibbsDraw1(zeroNy, 0, eta', z4vcvdraws(:,:,iterMCMC));
    else
        sqrtsigmaPREV = bayesSQRTVCVgibbsDraw1(SigmaT0, SigmaDof, eta', z4vcvdraws(:,:,iterMCMC));
    end
    %         OK    = abs(det(sqrtsigmaPREV)) > 1e-6;
    %         count = count + 1;
    %     end



    %% store draws post burnin
    if mm > 0
        Ydraws(:,:,mm)    = thisY;
        ETAdraws(:,:,mm)  = eta;
        sqrtSIGMA(:,:,mm) = sqrtsigmaPREV;


        YstateTmean(:,mm)  = thisYhat;
        YstateTvar(:,:,mm) = thisYvar;
    end

    if showProgress
        progressbar((mm + MCMCdraws + 1) / (2 * MCMCdraws + 1))
    end
end

%% collect SV

for mm = 1 : MCMCdraws
    thisB               = sqrtSIGMA(:,:,mm);
    ETASVdraws(:,:,mm)  = repmat(sqrt(diag(thisB * thisB')), 1, T);
end

%% simulate draws from predictive density

% allocate memory
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


    % collect one-step ahead moments for state vector
    YstateTp1mean(:,mm)    = F * YstateTmean(:,mm);
    YstateTp1var(:,:,mm)   = F * YstateTvar(:,:,mm) * F' + BB;


    % collect T+H-step ahead moments of outcome
    Ymean  = YstateTmean(:,mm);
    Yvar   = YstateTvar(:,:,mm);
    for jj = 1 : Nhorizons
        Ymean             = F * Ymean;
        Yvar              = F * Yvar * F' + BB;
        yTpHmean(1,jj,mm) = Ymean(1);
        yTpHvar(1,jj,mm)  = Yvar(1,1);
    end

    % construct eta
    z   = randn(rndStream,Ny,Nhorizons * Nfedraws);
    eta = B * z;
    eta = reshape(eta, Ny, Nhorizons, Nfedraws);
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

yTpHmean  = permute(yTpHmean, [3 2 1]);
yTpHvar   = permute(yTpHvar, [3 2 1]);

