function slopes = CGregressionpooledTC(PItilde,SIGMAtilde,sigTrend,Hpools)
% CGregressionpooled - Compute the slopes of pooled CG regressions in trend-cycle model
% slopes = CGregressionpooledTC(PItilde,SIGMAtilde,sigTrend,Hpools)
% PItilde are VAR slopes, SIGMAtlde is VCV of gap updates, sigtrend is ternd-shock variance and Hpools are cells of vectors listing the horizons to be pooled

Ngap    = size(PItilde,1);
I       = eye(Ngap);
Npools  = length(Hpools);
slopes  = NaN(Npools,1);
for nn = 1 : Npools
    poolH = Hpools{nn};
    yx = 0;
    xx = 0;
    for j = 1 : length(poolH)
        ii         = poolH(j) + 2;
        ePItilde   = zeros(1,Ngap);
        PItildeKp1 = I;
        
        for k = 0 : poolH(j)
            PItildeKp1    = PItildeKp1 * PItilde;
            ePItilde      = ePItilde + PItildeKp1(ii-k-1,:);
        end
        yx          = yx + ePItilde * SIGMAtilde(:,ii);
        xx          = xx + SIGMAtilde(ii,ii);
    end % j
    
    xx = xx + length(poolH) * sigTrend;

    slopes(nn) = yx ./ xx;
end % nn