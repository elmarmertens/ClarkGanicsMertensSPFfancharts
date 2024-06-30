function slopes = CGregressionTC(PItilde,SIGMAtilde,sigTrend,H)
% CG regression in population of trend-cycle model
% slopes = CGregressionTC(PItilde,SIGMAtilde,sigtrend,H)
% PItilde are VAR slopes, SIGMAtlde is VCV of gap updates, sigtrend is ternd-shock variance and H is the maximum horizon

Ngap    = size(PItilde,1);
I       = eye(Ngap);
slopes  = NaN(H,1);
for h = 1 : H
    ii          = h + 2;
    ePItilde    = zeros(1,Ngap);
    PItildeKp1  = I;
    
    for k = 0 : h
        PItildeKp1 = PItildeKp1 * PItilde;
        ePItilde   = ePItilde + PItildeKp1(ii-k-1,:);
    end
    yx          = ePItilde * SIGMAtilde(:,ii);
    xx          = SIGMAtilde(ii,ii) + sigTrend;
    slopes(h)   = yx ./ xx;
end % h
