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
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};

% select date

thisDate       = datenum(2022,4,1);

fontsize = 18;


%% some settings

YBARlabel = {'Next year', '2-years ahead', '3-years ahead'};
Nbar      = length(YBARlabel);
barx      = -5 : .01 : 20;


if quicky
    MCMCdraws    = 1e2;
    Nfedraws     = 10;
else
    MCMCdraws    = 1e3;
    Nfedraws     = 100;
end

%% loop over datalabels
for d = 1 : 5

    close all
    %% load dataset for variable
    datalabel = DATALABELS{d};

    datadir = fullfile('..', 'kensingtonDataMatfiles');
    matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
    load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', ...
        'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czhat', 'Czbar')

    ndxFixedHorizon     = 1 + 1 + (0:4);

    thisT    = find(dates <= thisDate, 1, 'last');
    thisDate = dates(thisT);
    thisDateLabel = datestr(thisDate, 'yyyyqq');

    % construct a model-label indicating dataset and important parameters, to
    % append to picture names
    modellabel = datalabel;
    %#ok<*UNRCH>

    modellabel = strcat(modellabel, 'STATEtrendgapSV');

    modellabel = strcat(modellabel, thisDateLabel);


    wrap = [];
    titlename = modellabel;
    if ~quicky
        initwrap
    end

    fprintf('Processing %s ... \n', modellabel)


    %% record location of missing values (NaN)

    samStart   = find(any(~Znanny,2),1); % for CPI/TBILL: leaves early sample with only Yrealized ...


    dates      = dates(samStart:end);
    T          = length(dates);
    Nhorizons  = Ny - 1;
    maxHorizon = Ny - 2;

    Zdata      = Zdata(samStart:end,:);
    Znanny     = Znanny(samStart:end,:);
    Cz         = Cz(:,:,samStart:end);



    %% call MCMC sampler
    [YdensityDraws, Ydraws, ETAdraws, ETASVdraws, ...
        SVscaleDraws, hvardraws, hrhodraws, sqrtsigmaDraws, mustarvardraws, ...
        YstateTp1mean, YstateTp1var] = ...
        mcmcsamplerSTATEtrendgapSV(Zdata, Znanny, Cz, ...
        thisT, MCMCdraws, Nfedraws, rndStream, true);


    %% adjust datevector etc for plotting
    datesT     = dates(1:thisT);



    %% plot hrho and hvar

    thisfig = figure;
    plotpriorposteriordraws(hrhodraws)
    wrapthisfigure(thisfig, sprintf('hrho%s', modellabel), wrap)
    thisfig = figure;
    plotpriorposteriordraws(hvardraws)
    wrapthisfigure(thisfig, sprintf('hvar%s', modellabel), wrap)


    %% plut mustarvol
    mustarvoldraws = sqrt(mustarvardraws);
    thisfig = figure;
    plotpriorposteriordraws(mustarvoldraws)
    wrapthisfigure(thisfig, sprintf('mustarvol%s', modellabel), wrap)


    %% plot SVscale
    thisfig = figure;
    hold on
    set(gca, 'fontsize', 18)
    plot(datesT,median(SVscaleDraws,2), 'r-', 'LineWidth', 2)
    plot(datesT,prctile(SVscaleDraws, [5 95], 2), 'r-', 'LineWidth', 1)
    xtickdates(datesT)
    wrapthisfigure(thisfig, sprintf('SVscale-%s', modellabel), wrap)


    %% collect VCV
    sqrtSIGMAdraws = NaN(Ny,Ny,thisT,MCMCdraws);
    Ngap = Ny - 1;
    parfor mm = 1 : MCMCdraws
        thisSIGMA                  = zeros(Ny,Ny,thisT);
        thisSIGMA(:,Ny,:)          = mustarvoldraws(mm);
        thisSV                     = permute(SVscaleDraws(:,mm), [2 3 1]);
        thisSIGMA(1:Ngap,1:Ngap,:) = sqrtsigmaDraws(:,:,mm) .* thisSV; % exploit array expansion
        sqrtSIGMAdraws(:,:,:,mm)   = thisSIGMA;

    end

    %% construct state space for univariate MA representation
    F          = diag(ones(Ny-1,1),1);
    F(end,end) = 1;
    e1         = zeros(1,Ny); e1(1) = 1;

    yMAprocess   = NaN(Ny,thisT,MCMCdraws);
    yMAprocessSD = NaN(Ny,thisT,MCMCdraws);
    parfor mm = 1 : MCMCdraws
        for t = 1 : thisT
            [SIG, K] = abckalman(F,sqrtSIGMAdraws(:,:,t,mm),e1);
            yMAprocess(:,t,mm)   = K;
            yMAprocessSD(:,t,mm) = K * sqrt(SIG(1));
        end
    end

    %% plot yMA results
    yMAmedian = median(yMAprocess, 3);

    %% plot hi-lo-mid
    SVmid = median(SVscaleDraws,2);
    [~, sortNdx] = sort(SVmid);
    loNdx  = sortNdx(1);
    midNdx = sortNdx(floor(thisT / 2));
    hiNdx  = sortNdx(end);

    maxlags = 20;
    irflags = 1 : Ny;
    extralags = irflags(end):maxlags;
    Nxtra     = length(extralags);
    thisfig = figure;
    hold on
    hlo = plot(irflags, yMAmedian(:,loNdx), 'color', Colors4Plots(1), 'linewidth', 2);
    plot(extralags, repmat(yMAmedian(end,loNdx), Nxtra, 1), 'color', Colors4Plots(1), 'linewidth', 2);
    hhi = plot(irflags, yMAmedian(:,hiNdx), 'color', Colors4Plots(2), 'linewidth', 2);
    plot(extralags, repmat(yMAmedian(end,hiNdx), Nxtra, 1), 'color', Colors4Plots(2), 'linewidth', 2);
    hmid = plot(irflags, yMAmedian(:,midNdx), 'color', Colors4Plots(3), 'linewidth', 2);
    plot(extralags, repmat(yMAmedian(end,midNdx), Nxtra, 1), 'color', Colors4Plots(3), 'linewidth', 2);

    xlim([1 maxlags])
    % xticks([-1 0 : 3 : Ny-2])
    grid on
    set(gca, 'FontSize', 18)
    wrapthisfigure(thisfig, sprintf('yMAprocessHiLoMid-%s', modellabel), wrap)
    legend([hlo, hmid, hhi], ...
        sprintf('%s (lowest SV)', datestr(dates(loNdx), 'yyyyQQ')), ...
        sprintf('%s (medium SV)', datestr(dates(midNdx), 'yyyyQQ')), ...
        sprintf('%s (highest SV)', datestr(dates(hiNdx), 'yyyyQQ')));
    wrapthisfigure(thisfig, sprintf('yMAprocessHiLoMid-%s-WITHLEGEND', modellabel), wrap)


    %% plot predictive density of Y
    Zdata(Znanny) = NaN;

    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end

    ndxB            = 2 + 4 + 1;
    horizonsB       = (4 - datesQ(thisT)) + MAweights;
    ndxC            = 2 + 4 + 2;
    horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
    ndxD            = 2 + 4 + 3;
    horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;

    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    % note: using the mean (to match observed)
    % plot density
    hsim = plot(0:Nhorizons-1,mean(YdensityDraws,2), 'r-', 'LineWidth', 4);
    plot(0:Nhorizons-1,prctile(YdensityDraws, normcdf([-1 1]) * 100, 2), 'r-', 'LineWidth', 2)
    plot(0:Nhorizons-1,prctile(YdensityDraws, [5 95], 2), 'r-.', 'LineWidth', 2)

    lightblue =  [0.5843    0.8157    0.9882];

    % plot SPF obs as available
    thisData = Zdata(thisT, 1);
    hData = plot(-1, thisData , 's-', 'color', [0 0 1], 'linewidth', 2);
    hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', lightblue, 'linewidth', 2);
    if ~Znanny(thisT,ndxB)
        hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', lightblue, 'linewidth', 2);
    end
    if ndxC <= Nz && ~Znanny(thisT,ndxC)
        hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', lightblue, 'linewidth', 2);
    end
    if ndxD <= Nz && ~Znanny(thisT,ndxD)
        hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', lightblue, 'linewidth', 2);
    end
    % hFuture = plot(0:Nhorizons-1, Yfuture(thisT,:), 's:', 'color', [0 .8 0], 'linewidth', 2);
    xticks([-1 0 : 2 : Nhorizons-1])
    xlim([-1 Nhorizons-1])
    % wrapthisfigure(thisfig, sprintf('YvsY2predictivedensity-%s', modellabel), wrap)
    hl = legend([hsim hSPF hData], 'Trend model', ...
        'SPF', 'y(-1)', ...
        'location', 'best');
    % wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGEND', modellabel), wrap)
    title(sprintf('%s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
    wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGENDTITLE', modellabel), wrap)

    %% ABS-ETA vs SV (trend scale only)

    etaMid      = median(ETAdraws,3)';
    etaSV       = median(ETASVdraws,3)';
    etaSVtails  = prctile(ETASVdraws, 100 * normcdf([-100 100]), 3);
    etaSVtails  = permute(etaSVtails, [2 1 3]);

    % pre COVID
    doPreCOVID = false;
    ndx = (1 : T) <= thisT;
    if doPreCOVID
        ndx        = ndx & (dates < datenum(2019,1,1));
        covidlabel = 'preCOVID';
    else
        covidlabel = '';
    end


    samEnd = length(dates(ndx));

    % plots
    for n = 0 : Nhorizons
        hanni = NaN(3,1);

        thisfig = figure;
        set(gca, 'FontSize', 18)


        % orient landscape
        if n <= 4
            barcol = .5 * [1 1 1];
        else
            barcol = .75 * [1 1 1];
        end

        hold on
        hanni(1) = bar(dates(ndx), abs(etaMid(ndx,n+1)), 1, 'FaceColor', barcol, 'EdgeColor', barcol);
        hanni(2) = plot(dates(ndx), etaSV(ndx,n+1), 'k-.', 'linewidth', 3);
        xtickdates(dates(1:samEnd))


        % legend(hanni, '|\eta|', 'SV', 'SV (scaleSV)')
        wrapthisfigure(thisfig,sprintf('etaSV%d-%s%s', n, datalabel, covidlabel), wrap, [], [], [], [], false)
    end
    dockAllFigures
    finishwrap

end

%% finish / clean up
finishscript
