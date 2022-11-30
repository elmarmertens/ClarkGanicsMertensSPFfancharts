%% plot SV figures for Baseline model (ETA SV) vs an alternative

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/

%% clear workspace
clear variables
 
close all
fclose all;
clc

wrap = [];


%% select choice of data

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

modeltype      = 'STATEtrendgapSV';
titlename      = strcat('ETA4SV-', modeltype);
initwrap

doSince1992    = false;

resultdir      = localresultsMCMC;
Ndraws         = 3e3;

%#ok<*UNRCH>
for doPreCOVID = false % [true false]


    if doPreCOVID
        covidlabel = '-preCOVID';
    else
        covidlabel = '';
    end

    %% loop over data lables
    for d = 1 : length(DATALABELS)

        datalabel = DATALABELS{d};


        modellabel  = strcat(datalabel, modeltype);
        %% load data

        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
        mat = matfile(fullfile(resultdir, matfilename));


        Nhorizons  = mat.Nhorizons;
        Ylabel     = mat.Ylabel;
        dates      = mat.dates;

        T1992      = find(dates == datenum(1992,1,1));

        fprintf('Processing %s ... \n', modellabel)


        etaMid      = median(mat.ETAdraws,3)';
        etaFinalMid = median(mat.ETAFINALdraws,3)';

        etaSV       = median(mat.ETASVdraws,3)';
        etaSVtails  = prctile(mat.ETASVdraws, 100 * normcdf([-100 100]), 3);
        etaSVtails  = permute(etaSVtails, [2 1 3]);
        etaSVFINAL  = median(mat.ETASVFINALdraws,3)';

        %% pre COVID
        if doPreCOVID
            ndx = dates < datenum(2019,1,1);
        else
            ndx = true(size(dates));
        end


        samEnd = length(dates(ndx));

        %% abs ETA vs SV
        for n = 0 : Nhorizons
            hanni = NaN(2,1);

            thisfig = figure;
            set(gca, 'FontSize', 18)


            % orient landscape
            hold on
            if n <= 4
                barcol = .5 * [1 1 1];
            else
                barcol = .75 * [1 1 1];
            end
            hanni(1) = bar(dates(ndx), abs(etaFinalMid(ndx,n+1)), 1, 'FaceColor', barcol, 'EdgeColor', barcol);
            hanni(2) = plot(dates(ndx), etaSVFINAL(ndx,n+1), 'k-', 'linewidth', 4); 


            xtickdates(dates(1:samEnd))
            maxy = max(ylim);
            ylim([0 maxy])

            wrapthisfigure(thisfig,sprintf('etaSV%d-%s-%s%s', n, datalabel, modeltype, covidlabel), wrap, [], [], [], [], false)
            hl = legend(hanni, '|{\bf\eta}|', 'SV', 'location', 'northwest');
            wrapthisfigure(thisfig,sprintf('etaSV%d-%s-%s%s-WITHLEGEND', n, datalabel, modeltype, covidlabel), wrap, [], [], [], [], false)

            if doSince1992
                delete(hl)
                xtickdates(dates(T1992:samEnd))
                wrapthisfigure(thisfig,sprintf('etaSV%d-%s-%s%s-since1992', n, datalabel, modeltype, covidlabel), wrap, [], [], [], [], false)
                hl = legend(hanni, '|{\bf\eta}|', 'SV', 'location', 'northwest');
                wrapthisfigure(thisfig,sprintf('etaSV%d-%s-%s%s-since1992-WITHLEGEND', n, datalabel, modeltype, covidlabel), wrap, [], [], [], [], false)
            end

        end


        %% prepare next loop
        if ~isempty(wrap)
            close all
        end
    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures

