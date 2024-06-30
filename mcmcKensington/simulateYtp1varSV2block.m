function Ytp1var = simulateYtp1varSV2block(ndxSV11, ndxSV22, Nstates, Nfedraws, StateTvar, Bgap, SV, mustarVariance, F)
% Ytp1var = NaN(Nstates, Nstates, Nfedraws);

Ngap = Nstates - 1;

YvarTT   = F * StateTvar * F';
Bgap1    = Bgap(:, ndxSV11) .* permute(SV(1,:,1), [1 3 2]);
Bgap2    = Bgap(:, ndxSV22) .* permute(SV(2,:,1), [1 3 2]);
BGAP     = cat(2, Bgap1, Bgap2);
BBGAP    = pagemtimes(BGAP, "none", BGAP, "transpose");

BB                      = zeros(Nstates, Nstates, Nfedraws);
BB(Nstates, Nstates, :) = mustarVariance;
BB(1:Ngap, 1:Ngap, :)   = BBGAP;

Ytp1var = YvarTT + BB;
end