function [densityDraws, yhatRB] = simulateVAR0(Nstates, Nfedraws, Nhorizons, ...
    Bgap, SV, mustarvariance, F, Y, G, YgapConst, ETA, rndStream)

Yconst             = cat(1, YgapConst, 0);
Ngap               = Nstates - 1;
shocksGAP          = randn(rndStream,Ngap,Nfedraws,Nhorizons);
shocksTrend        = randn(rndStream,1,Nfedraws,Nhorizons);

shocksGAP   = SV .* shocksGAP;
shocksGAP   = pagemtimes(Bgap, shocksGAP);
shocksTrend = sqrt(mustarvariance) .* shocksTrend;
shocks      = cat(1,shocksGAP, shocksTrend);


shocks       = cat(2, shocks, -1 .* shocks); % antithetic simulation
densityDraws = NaN(2 * Nfedraws, Nhorizons);
thisY        = Y - Yconst;
thisETA      = ETA;
for jj = 1 : Nhorizons
    thisETA            = G * thisETA + shocks(:,:,jj);
    thisY              = F * thisY + thisETA;
    densityDraws(:,jj) = thisY(1,:) + thisY(Nstates,:) + Yconst(1);
end

if nargout > 1
    yhatRB = NaN(Nhorizons, 1);
    thisY   = Y - Yconst;
    thisETA = ETA;
    for jj = 1 : Nhorizons
        thisETA         = G * thisETA;
        thisY           = F * thisY + thisETA;
        yhatRB(jj)      = thisY(1,:) + thisY(Nstates,:) + Yconst(1);
    end
end
