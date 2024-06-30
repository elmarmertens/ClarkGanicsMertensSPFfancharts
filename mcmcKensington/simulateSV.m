function SV = simulateSV(rho, hT, hvol, Nfedraws, Nhorizons,rndStream)
hshocks = hvol * randn(rndStream, 1, Nfedraws, Nhorizons);
h = NaN(1, Nfedraws, Nhorizons);
j = 1;
h(:,:,j) = rho .* hT + hshocks(:,:,j);
for j = 2 : Nhorizons
    h(:,:,j) = rho .* h(:,:,j-1) + hshocks(:,:,j);
end
SV = exp(h * .5);