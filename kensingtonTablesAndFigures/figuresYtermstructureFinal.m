%% compares terms structure of expectations from two models

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

DATALABELS = {'UNRATE', 'RGDP', 'TBILL', 'CPI', 'PGDP'};

modeltype   = 'STATEtrendgapSV';
modelpretty = 'SV model';
Ndraws      = 3e3;
resultdir   = localresultsMCMC;

% modeltype   = 'STATEtrendgapSVshortlong';
% modelpretty = 'SV model (short/long)';
% Ndraws      = 1e3;
% resultdir   = '~/jam/lager/KENSINGTON/kensingtonresultsShortlong/'; 


quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
fontsize  = 24;

doShowSPFdata = true;

colorSPF1y = Colors4Plots(1);
colorSPF2y = Colors4Plots(5);
colorSPF3y = Colors4Plots(2);
colorSPFquarterly =  [1 0 0];
colorYterm  = .3 * [1 1 1]; % Colors4Plots(7);
colorYterm2 = Colors4Plots(3);

%% loop over datalabels
for d =  1 : length(DATALABELS)

    %% prepare things

    %     tndx = [];

    tndx   = 164 : 200; % since 2009
    % tndx = 194 : 215; % since 2017
    % tndx = 206 : 215; % since 2020Q1

    datalabel   = DATALABELS{d};


    modellabel  = strcat(datalabel, modeltype);
    %% prepare latexwrapper
    wrap = [];
    titlename = sprintf('termstructureFinal-%s-%s', datalabel, modeltype);
    %     initwrap
    if isempty(wrap) && ~isdesktop
        initwrap
    end


    fprintf('Processing %s ... \n', datalabel)

    %% load data

    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
    mat = matfile(fullfile(resultdir, matfilename));



    %     YdensityDraws = mat.YdensityDraws;
    Ydraws        = mat.Ydraws;
    YFINALdraws   = mat.YFINALdraws;

    Yfuture       = mat.Yfuture;
    Nz            = mat.Nz;
    Zdata         = mat.Zdata;
    Znanny        = mat.Znanny;

    Nhorizons     = mat.Nhorizons;
    dates         = mat.dates;
    T             = mat.T;
    datesQ        = mat.datesQ;
    doNIPA        = mat.doNIPA;

    Tstart        = mat.Tstart;



    if doNIPA
        MAweights       = -2 : 4;
    else
        MAweights       = 1 : 4;
    end
    ndxFixedHorizon     = 1 + 1 + (0:4);

    %% plot term structures
    if isempty(tndx)
        tndx = Tstart : T;
    end

    for thisT = tndx

        ndxB            = 2 + 4 + 1;
        horizonsB       = (4 - datesQ(thisT)) + MAweights;
        ndxC            = 2 + 4 + 2;
        horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
        ndxD            = 2 + 4 + 3;
        horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;

        these  = squeeze(Ydraws(2:end,thisT,:));
        mid1   = mean(these, 2);
        tails1 = prctile(these, normcdf([-1 1]) * 100, 2);

        these  = squeeze(YFINALdraws(2:end,thisT,:));
        mid2   = mean(these, 2);
        tails2 = prctile(these, normcdf([-1 1]) * 100, 2);

        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)

        h1 = plot(0:Nhorizons-1, mid1, '-', 'color', colorYterm, 'linewidth', 3);
        %         plot(0:Nhorizons-1, tails1, '-', 'color', colorYterm, 'linewidth', 1);

        h2 = plot(0:Nhorizons-1, mid2, '-', 'color', colorYterm2, 'linewidth', 3);
        %         plot(0:Nhorizons-1, tails2, '-', 'color', colorYterm2, 'linewidth', 1);

        ylimits = ylim;
        switch upper(datalabel)
            case {'CPI'}
                if all(ylimits > 0)
                    ylimits(1) = 0;
                elseif all(ylimits < 0)
                    ylimits(2) = 0;
                end
            case 'UNRATE'
                %                 if all(ylimits > 0)
                %                     ylimits(1) = 0;
                %                 elseif all(ylimits < 0)
                %                     ylimits(2) = 0;
                %                 end
                %                 if year(dates(thisT)) == 2020
                %                     ylimits(2) = 20;
                %                 end
            case 'TBILL'
                if all(ylimits > 0)
                    ylimits(1) = 0;
                elseif all(ylimits < 0)
                    ylimits(2) = 0;
                end
                if ylimits(2) < 4
                    ylimits(2) = 4;
                end
                if year(dates(thisT)) <= 2012
                    ylimits(2) = 4;
                end        
                if year(dates(thisT)) == 2020
                    ylimits(2) = 2.5;
                end                    
                if year(dates(thisT)) == 2021
                    ylimits(2) = 2.5;
                end                    
            case {'RGDP', 'PGDP'}
                if year(dates(thisT)) == 2009
                    ylimits(1) = -6;
                    ylimits(2) = 6;
                elseif year(dates(thisT)) > 2019
                    ylimits(1) = -10;
                    ylimits(2) = 10;
                else
                    ylimits(1) = 1;
                    ylimits(2) = 4.5;
                end

                
                % if (year(dates(thisT)) >= 2017) && (year(dates(thisT)) <= 2019)
                %    ylimits(1) = min(1, ylimits(1));
                %    ylimits(2) = max(3, ylimits(2));
                % end  


        end
        ylim(ylimits)
        xlim([-1 Nhorizons-1])
        wrapthisfigure(thisfig, sprintf('Ytermdraws-%s%s-%s', datestr(dates(thisT), 'yyyyqq'), ...
            datalabel, modeltype), wrap, [], [], [], [], true);
        if doShowSPFdata
           
            hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'd', 'color', colorSPFquarterly, 'linewidth', 3);
            spfLegend = {'SPF quarterly'};
            hl = legend([h1, h2, hSPF], modelpretty, sprintf('%s (full sample)', modelpretty), spfLegend{:}, 'location', 'best');
            ylim(ylimits)
            wrapthisfigure(thisfig, sprintf('Ytermdraws-%s%s-%s-WITHQUARTERLYSPF-WITHLEGEND', datestr(dates(thisT), 'yyyyqq'), ...
                datalabel, modeltype), wrap, [], [], [], [], false);
     
            if ~Znanny(thisT,ndxB)
                hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':', 'color', colorSPF1y, 'linewidth', 3);
                hSPF  = [hSPF, hSPFb]; %#ok<AGROW> 
                spfLegend = cat(1, spfLegend, 'SPF next year');
            end
            if ndxC <= Nz && ~Znanny(thisT,ndxC)
                hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':', 'color', colorSPF2y, 'linewidth', 3);
                hSPF  = [hSPF, hSPFc];  %#ok<AGROW> 
                spfLegend = cat(1, spfLegend, 'SPF 2 years ahead');
            end
            if ndxD <= Nz && ~Znanny(thisT,ndxD)
                hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':', 'color', colorSPF3y, 'linewidth', 3);
                hSPF  = [hSPF, hSPFd];  %#ok<AGROW>
                spfLegend = cat(1, spfLegend, 'SPF 3 years ahead');
            end
            hl = legend([h1, h2, hSPF], modelpretty, sprintf('%s (full sample)', modelpretty), spfLegend{:}, 'location', 'best');
            ylim(ylimits)
            wrapthisfigure(thisfig, sprintf('Ytermdraws-%s%s-%s-WITHSPF-WITHLEGEND', datestr(dates(thisT), 'yyyyqq'), ...
                datalabel, modeltype), wrap, [], [], [], [], false);
            delete(hl)
            ylim(ylimits)
            wrapthisfigure(thisfig, sprintf('Ytermdraws-%s%s-%s-WITHSPF', datestr(dates(thisT), 'yyyyqq'), ...
                datalabel, modeltype), wrap, [], [], [], [], true);
            
        end

        %% clean up
        if ~isempty(wrap)
            close(thisfig);
        end
    end


    %% finish datalabel loop
    finishwrap


end

%% finish / clean up
finishscript
dockAllFigures
