function INFF = plotINFF(thesedraws, thislabel, modellabel, wrap)
if nargin < 2
    thislabel = '';
end
if nargin < 3
    modellabel = '';
end
if nargin < 4
    wrap = [];
end
taperlabels = {'4%', '8%', '15%'};
INFF = ineff(thesedraws);
if size(INFF,1) > 10 % histogram
    thisfig = figure;
    for ii = 1 : size(INFF, 2)
        subplot(size(INFF, 2), 1, ii);
        histogram(INFF(:, ii)); %
        title(sprintf('INEFF (%s tapering)', taperlabels{ii}));
    end
    sgtitle(sprintf('%s (histograms)', thislabel));
    wrapthisfigure(thisfig, strcat('INEFF-', modellabel), wrap);
else % barplots
    thisfig = figure;
    for ii = 1 : size(INFF, 2)
        subplot(size(INFF, 2), 1, ii);
        bar(INFF(:, ii));
        title(sprintf('INEFF (%s tapering)', taperlabels{ii}));
    end
    sgtitle(thislabel);
    wrapthisfigure(thisfig, strcat('INEFF-', modellabel), wrap);
end
