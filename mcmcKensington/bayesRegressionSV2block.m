function [GPREV, residGAP] = bayesRegressionSV2block(GPREV, etaGAP, etaGAPlag, Ngibbs, invSVt, invsigma11, invsigma22, K, Ptilde11, P22, ndxSV11, ndxSV22, Igap, rndStream)

etaGAPtilde2 = etaGAP(:,ndxSV22);
etaGAPtilde1 = etaGAP(:,ndxSV11) - etaGAPtilde2 * transpose(K);

Kappa          = kron(K, Igap);
Gprime         = transpose(GPREV);
Gtilde2        = Gprime(:,ndxSV22);

for nn = 1 : Ngibbs

    % step 1
    mu12           = -Kappa * Gtilde2(:);
    % Ptilde11       = P11;
    Ptilde11mu12   = Ptilde11 * mu12;
    lhs            = etaGAPtilde1 .* invSVt(:,1);
    rhs            = etaGAPlag .* invSVt(:,1);
    Gtilde1        = bayesVectorRegressionGibbsDraw1(lhs, rhs, invsigma11, Ptilde11mu12, Ptilde11, rndStream);
    
    % step2
    Ptilde22       = Kappa' * Ptilde11 * Kappa + P22;
    Ptilde22mu21   = - Kappa' * Ptilde11 * Gtilde1(:);
    lhs            = etaGAPtilde2 .* invSVt(:,2);
    rhs            = etaGAPlag .* invSVt(:,2);
    Gtilde2        = bayesVectorRegressionGibbsDraw1(lhs, rhs, invsigma22, Ptilde22mu21, Ptilde22, rndStream);
end

% reconstruct GPREV
G2       = Gtilde2;
G1       = Gtilde1 + Gtilde2 * transpose(K);
GPREV    = [G1 G2];
residGAP = etaGAP - etaGAPlag * GPREV;

