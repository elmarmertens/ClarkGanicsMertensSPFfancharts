function densityDraws = simulateMDS2block(ndxSV11, ndxSV22, Nstates, Nfedraws, Nhorizons, Bgap, SV, mustarvariance, F, Y, rndStream)
Ngap               = Nstates - 1;

shocksGAP              = randn(rndStream,Ngap,Nfedraws,Nhorizons);
shocksGAP(ndxSV11,:,:) = SV(1,:,:) .* shocksGAP(ndxSV11,:,:);
shocksGAP(ndxSV22,:,:) = SV(2,:,:) .* shocksGAP(ndxSV22,:,:);
shocksGAP              = pagemtimes(Bgap, shocksGAP);

shocksTrend  = randn(rndStream,1,Nfedraws,Nhorizons);
shocksTrend  = sqrt(mustarvariance) .* shocksTrend;

shocks       = cat(1,shocksGAP, shocksTrend);

shocks       = cat(2, shocks, -1 .* shocks); % antithetic simulation
densityDraws = NaN(2 * Nfedraws, Nhorizons);
for jj = 1 : Nhorizons
    Y                  = F * Y + shocks(:,:,jj); % note: expansion into Ny x Nfedraws array at jj=1
    densityDraws(:,jj) = Y(1,:) + Y(Nstates,:);
end
