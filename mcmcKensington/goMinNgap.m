%% STATE-TREND-SCALE-SV model -- single run for specific thisT
% uses MDS assumption

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*PFBNS>
%#ok<*UNRCH>
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc
rng('default')

%% some parameters

rndStream  = getDefaultStream;

% beginning and end of eval window (set to empty if to use max)

quicky     = false; % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'};


% select date
thisDate       = datenum(2023,10,1);

fontsize = 18;

% SED-PARAMETERS-HERE


%% some settings

if quicky
    MCMCdraws    = 1e2;
    Nfedraws     = 10;
else
    MCMCdraws    = 1e3;
    Nfedraws     = 100;
end

initwrap

%% loop over datalabels
for d = 1 : length(DATALABELS)

    % close all
    %% load dataset for variable
    datalabel = DATALABELS{d};

    datadir = fullfile('..', 'matdataKensington');
    matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
    load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', ...
        'YBARlabel', 'Nbar', ...
        'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czhat', 'Czbar')

    ndxFixedHorizon     = 1 + 1 + (0:4);


    fprintf('\n')
    fprintf('Processing %s ... \n', datalabel)

    %% record location of missing values (NaN)

    samStart   = find(sum(~Znanny,2) > 1,1); % first availability of SPF

    startDateLabel = datestr(dates(samStart), 'yyyyqq');
    dates      = dates(samStart:end); datesQ = datesQ(samStart:end);
    T          = length(dates);

    Zdata      = Zdata(samStart:end,:);
    Znanny     = Znanny(samStart:end,:);
    Cz         = Cz(:,:,samStart:end);

    %% cloop over Ngaps to find minimal compliant value
    NGAPS = 1 : Ny;

    % working with some cellfuns to avoid loops
    OK     = false(T,length(NGAPS));
    obsndx = num2cell(~Znanny, 2);

    for n = 1 : length(NGAPS)

        Ngap    = NGAPS(n);

        thisC   = Cz(:,1:Ngap,:);
        thisC   = cat(2,thisC,ones(Nz,1,T));

        % convert to cell and remove missing-data rows
        thisC   = squeeze(num2cell(thisC, [1, 2]));
        thisC   = cellfun(@(x,n) x(n,:), thisC, obsndx, 'UniformOutput', false);

        OK(:,n) = cellfun(@(x,y) rank(x) == sum(y), thisC, obsndx);

    end

    %% find minNgap
    % for each row of OK, find the column with the first element equal true
    minimalNgap = NaN(T,1);
    for tt = 1 : T
        minimalNgap(tt) = NGAPS(find(OK(tt,:), 1, 'first'));
    end

    %% plot minNgap (all dates)
    thisfig = figure;
    % plot(dates, minNgap, '-', 'LineWidth', 2)
    bc = colors4plots('blue');
    bar(dates, minimalNgap, 1, 'FaceColor', bc, 'EdgeColor', bc)
    xtickdates(dates, 'keeplimits');
    title(sprintf('minNgap for %s', datalabel))
    wrapthisfigure(thisfig, sprintf('minNgap%s', upper(datalabel)), wrap)

    %% plot minNgap Q1 dates only
    ndxQ1 = quarter(dates) == 1;
    thisfig = figure;
    % plot(dates, minNgap, '-', 'LineWidth', 2)
    bc = colors4plots('blue');
    bar(dates(ndxQ1), minimalNgap(ndxQ1), 1, 'FaceColor', bc, 'EdgeColor', bc)
    xtickdates(dates(ndxQ1), 'keeplimits');
    title(sprintf('minNgap (in Q1) for %s', datalabel))
    wrapthisfigure(thisfig, sprintf('minNgap%sQ1', upper(datalabel)), wrap)

    %% tabulate min values
    datesQ1      = dates(ndxQ1);
    [uval, urow] = unique(minimalNgap(ndxQ1));
    fprintf('MinNgaps (Q1 dates only):\n'); 
    for j = 1 : length(uval)
        fprintf('%d \t %s\n', uval(j), datestr(datesQ1(urow(j))));
    end
end

%% finish / clean up
dockAllFigures
finishwrap
finishscript
