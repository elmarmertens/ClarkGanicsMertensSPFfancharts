function Ytp1var = simulateYtp1var(Nstates, Nfedraws, StateTvar, Bgap, SV, mustarVariance, F)
% Ytp1var = NaN(Nstates, Nstates, Nfedraws);

Ngap = Nstates - 1;

YvarTT   = F * StateTvar * F';
BBgap    = Bgap * Bgap';
SV2      = permute(SV(:,:,1).^2, [1 3 2]); % permutation prepares vectorized operation below

BB                      = zeros(Nstates, Nstates, Nfedraws);
BB(Nstates, Nstates, :) = mustarVariance;
BB(1:Ngap, 1:Ngap, :)   = SV2 .* BBgap;

Ytp1var = YvarTT + BB;
end