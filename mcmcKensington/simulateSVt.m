function SVt = simulateSVt(rho, hT, hvol, tdof, Nfedraws, Nhorizons,rndStream)

Nsv = length(rho);
if Nsv == 1
    hshocks = hvol * randn(rndStream, 1, Nfedraws, Nhorizons);
else
    hshocks = hvol * randn(rndStream,Nsv,Nhorizons * Nfedraws); % hvol is assumed to be square-root factor of VCV
    hshocks = reshape(hshocks, Nsv, Nfedraws, Nhorizons);
end
h = NaN(Nsv, Nfedraws, Nhorizons);
j = 1;
h(:,:,j) = rho .* hT + hshocks(:,:,j);
for j = 2 : Nhorizons
    h(:,:,j) = rho .* h(:,:,j-1) + hshocks(:,:,j);
end
chi2draws  = chi2rnd(repmat(tdof, [1, Nfedraws, Nhorizons]));
tscalelog2 = log(tdof) - log(chi2draws);
SVt        = exp((h + tscalelog2) * .5);