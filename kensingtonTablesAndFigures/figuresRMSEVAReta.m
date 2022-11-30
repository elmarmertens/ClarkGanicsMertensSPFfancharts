%% compute and plot MSE decomposition of SPF-consistent expectations of VAReta

% todo:
% - load predictive means directly (once stored with MCMC)
% - in MCMC: compute predictive means analytically (instead of MC)
% - compute bias directly in MCMC to generate also uncertainty around bias

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*ASGLU>
%#ok<*NASGU>
%#ok<*UNRCH>

%% clear workspace
clear variables

close all
fclose all;
clc


%% some parameters


modeltype   = 'VARprior10STATEtrendgapSV';
modelpretty = 'non-MDS (SV)';


Ndraws      = 3e3;

resultdir = localresultsMCMC;

fontsize = 18;

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};

%% initwrap
wrap = [];
titlename = sprintf('SPFmse-%s', modeltype);

if isempty(wrap) && ~isdesktop
    initwrap
end

%% loop over datalabels

for d = 1 : length(DATALABELS)

    close all
    datalabel   = DATALABELS{d};
    modellabel  = strcat(datalabel, modeltype);

    fprintf('Processing %s ... \n', datalabel)

    %% load data
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));

    Yhat          = mat.fcstYhatRB;
    Ydraws        = mat.Ydraws;
    Yvol          = mat.fcstYvol;


    Yfuture       = mat.Yfuture;
    Nz            = mat.Nz;
    Zdata         = mat.Zdata;
    Znanny        = mat.Znanny;

    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    Tstart        = mat.Tstart;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    horizonLabels = arrayfun(@(x) sprintf('h=%d', x), 0:Nhorizons-1, 'uniformoutput', false);

    T2019 = find(dates == datenum(2019,10,1));

    %% compute bias
    spfYhat  = mean(Ydraws(2:end,:,:),3)';
    spfBias  = Yhat - spfYhat;
    spfBias2 = spfBias.^2;
    spfYvar  = Yvol.^2;
    spfMSE   = spfBias2 + spfYvar;

    for doPreCOVID = true % [true false]
        if doPreCOVID
            covidLabel = '-preCOVID';
            Tend       = T2019;
        else
            covidLabel = '';
            Tend       = T;
        end

        %% plot MSE decomposition

        for hh = 1 : Nhorizons
            thisfig = figure;
            hold on
            hanni  = bar(dates(Tstart:Tend) , [spfYvar(Tstart:Tend, hh) spfBias2(Tstart:Tend, hh)], 1, 'stacked');
            hanni(1).FaceColor = [1 1 1] * .7;
            hanni(1).EdgeColor =  hanni(1).FaceColor;
            hanni(2).FaceColor = [1 1 1] * 0;
            hanni(2).EdgeColor =  hanni(2).FaceColor;
            box off
            set(gca, 'FontSize', fontsize)
            nberlines(dates(Tstart:Tend))
            xticks(datenum(1970:5:2020,1,1))
            xtickdates(dates(Tstart:Tend), 'keepticks')
            wrapthisfigure(thisfig, sprintf('SPFmsebars-h%d-%s-%s%s', hh-1, ...
                datalabel, modeltype, covidLabel), wrap, [], [], [], [], true);
            legend(hanni, 'Forecast Error Variance', 'Squared Bias', 'location', 'northwest')
            wrapthisfigure(thisfig, sprintf('SPFmsebars-h%d-%s-%s%s-WITHLEGEND', hh-1, ...
                datalabel, modeltype, covidLabel), wrap, [], [], [], [], false);
        end

        %% plot MSE decomposition - yyplot
        for hh = 1 : Nhorizons
            thisfig = figure;
            hold on
            [axhanni, hVar, hBias] = plotyy(dates(Tstart:Tend) , spfYvar(Tstart:Tend, hh), ...
                dates(Tstart:Tend) , spfBias2(Tstart:Tend, hh)); %#ok<PLOTYY> 
            hVar.LineStyle  = '-.';
            hVar.LineWidth  = 2;
            hBias.LineWidth = 1;
            
            box off
            set(axhanni, 'FontSize', fontsize)
            %             yyaxis left
            nbershades(dates(Tstart:Tend))
            xticks(datenum(1970:5:2020,1,1))
            xtickdates(dates(Tstart:Tend), 'keepticks')
            legend([hVar, hBias], 'Forecast Error Variance', 'Squared Bias', 'location', 'northwest')
            title(sprintf('MSE components %s h=%d', datalabel, hh-1))
            wrapthisfigure(thisfig, sprintf('SPFmseyy-h%d-%s-%s%s-WITHLEGEND', hh-1, ...
                datalabel, modeltype, covidLabel), wrap, [], [], [], [], false);
        end


        %% plot selected Bias
        horizons      = 0:Nhorizons-1;
        HorizonLabels = arrayfun(@(x) sprintf('h=%d', x), horizons, 'uniformoutput', false);

        ndxHnear = [1 4];
        ndxH     = unique([8 : 4 : Nhorizons, Nhorizons]); % note: unique also does an automatic sort

        % full sample
        thisfig = figure;
        hold on
        set(gca, 'FontSize', fontsize)
        box off
        hanniNear = plot(dates(Tstart:Tend), spfBias(Tstart:Tend,ndxHnear), '-.', 'linewidth', 2);
        hanni = plot(dates(Tstart:Tend), spfBias(Tstart:Tend,ndxH), 'linewidth', 2);
        nbershades(dates(Tstart:Tend))
        xticks(datenum(1970:5:2020,1,1))
        xtickdates(dates(Tstart:Tend), 'keepticks')
        wrapthisfigure(thisfig, sprintf('SPFbias-%s%s', modellabel, covidLabel), wrap, [], [], [], [], true);
        legend([hanniNear; hanni], HorizonLabels([ndxHnear ndxH]), 'Location','best')
        wrapthisfigure(thisfig, sprintf('SPFbias-%s%s-WITHLEGEND', modellabel, covidLabel), wrap)

    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures
