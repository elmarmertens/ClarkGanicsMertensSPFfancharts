function R = plotPSRF(thesedraws, thislabel, modellabel, wrap)
if nargin < 2
    thislabel = '';
end
if nargin < 3
    modellabel = '';
end
if nargin < 4
    wrap = [];
end
R = psrf(thesedraws);
if length(R) > 10 % histogram
    thisfig = figure;
    histogram(R);
    title(sprintf('%s: PSRF (histogram)', thislabel));
    wrapthisfigure(thisfig, strcat('PSRF-', modellabel), wrap);
else % barplot
    thisfig = figure;
    bar(R);
    title(sprintf('%s: PSRF', thislabel));
    wrapthisfigure(thisfig, strcat('PSRF-', modellabel), wrap);
end
