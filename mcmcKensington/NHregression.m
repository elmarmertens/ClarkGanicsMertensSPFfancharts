function slopes = NHregression(G,PSI,H)
% NHregression - computes the slopes of the Norhdaus ("NH") regression
% slopes = NHregression(G,PSI,H)
% G are VAR slopes, PSI is VCV of forecast updates and H is max horizon

slopes  = NaN(H,1);
for h = 1 : H
    ii          = h + 2;
    xx          = PSI(ii,ii);
    yx          = G(ii-1,:) * PSI(:,ii);
    slopes(h)   = yx ./ xx;
end % h
