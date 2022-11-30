%% forecast fan charts of STATE-SCALE-SV model
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
MODELTYPES = {'trendgapSV', 'trendgapSVnoise2', 'const'};

Ndraws    = 3e3;

resultdir = localresultsMCMC;

quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
ndx68     = [3 7];
ndx90     = [2 8];
ndx6890   = [2 3 7 8];
fontsize  = 18;

doShow90  = false;

if doShow90
    show90label = '-bands90';
else
    show90label = '';
end

dataColor   = [0.4660 0.6740 0.1880]; % medium green
modelColor  = [1 0 0]; % [0.6350 0.0780 0.1840];  % dark red


%% loop over variables

for m = 1 % : 2 % length(MODELTYPES)

    modeltype = MODELTYPES{m};

    for d = 1 : length(DATALABELS)

        %% prepare things

        %         tndx = 206 : 215; % since 2022Q1
        tndx = 162 : 215; % since 2022Q1

        close all
        datalabel = DATALABELS{d};

        modellabel = strcat(datalabel, 'STATE', modeltype);

        %% load data
        fprintf('Processing %s... \n', modellabel);
        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d', modellabel, Ndraws);
        mat = matfile(fullfile(resultdir, matfilename));

        Yhat          = mat.fcstYhatRB;
        Yquantiles    = mat.fcstYquantiles;
        Yfuture       = mat.Yfuture;
        Nz            = mat.Nz;
        Zdata         = mat.Zdata;
        Znanny        = mat.Znanny;

        Nhorizons     = mat.Nhorizons;
        dates         = mat.dates;

        datesQ        = mat.datesQ;
        doNIPA        = mat.doNIPA;

        % patch up firstT and lastT from previos matfile formats
        Tstart           = mat.Tstart;
        T                = mat.T;

        %% prepare latexwrapper
        wrap = [];
        titlename = sprintf('fancharts-%s', modellabel);
        initwrap

        %% plot predictive density of Y


        if doNIPA
            MAweights       = -2 : 4;
        else
            MAweights       = 1 : 4;
        end

        ndxFixedHorizon     = 1 + 1 + (0:4);


        if isempty(tndx)
            if ismac && isdesktop
                tndx = 166; % T - 12;
            else
                tndx = Tstart : T;
            end
        end

        tndx = tndx(tndx >= Tstart);

        for thisT = tndx


            ndxB            = 2 + 4 + 1;
            horizonsB       = (4 - datesQ(thisT)) + MAweights;
            ndxC            = 2 + 4 + 2;
            horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
            ndxD            = 2 + 4 + 3;
            horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;


            mid     = Yhat(thisT,:);
            tails68 = squeeze(Yquantiles(thisT,:,ndx68));
            tails90 = squeeze(Yquantiles(thisT,:,ndx90));


            thisfig = figure;
            hold on
            set(gca, 'fontsize', fontsize)

            % model density
            hsim = plot(0:Nhorizons-1,mid, ':', 'color', modelColor, 'LineWidth', 5);
            plot(0:Nhorizons-1,tails68, '-', 'color', modelColor, 'LineWidth', 3)
            if doShow90
                plot(0:Nhorizons-1,tails90, '-', 'color', modelColor, 'LineWidth', 2)
            end

            xticks([-1 0 : 2 : Nhorizons-1])
            xlim([0 Nhorizons-1])


            ylimits = ylim;
            switch upper(datalabel)
                case {'CPI'}
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    elseif all(ylimits < 0)
                        ylimits(2) = 0;
                    end
                    if dates(thisT) >= datenum(2020,1,1)
                        ylimits = [-4 10];
                    end
                case 'UNRATE'
                    switch year(dates(thisT))
                        case {2020,2021,2022}
                            ylimits = ylim;
                            if ylimits(1) >= 0
                                ylimits(1) = 0;
                            end
                            if ylimits(2) <= 15
                                ylimits(2) = 15;
                            end
                        otherwise
                            ylimits = ylim;
                    end
                    if all(ylimits > 0)
                        ylimits(1) = 0;
                    end
                case 'RGDP'
                    switch year(dates(thisT))
                        case 2020
                            ylimits = [-50 40];
                        case {2021,2022}
                            ylimits = [-15 20];
                        otherwise
                            ylimits = ylim;
                    end
            end
            ylim(ylimits)

            wrapthisfigure(thisfig, sprintf('Ypredictivedensity%s-%s%s', datestr(dates(thisT), 'yyyyqq'), ...
                modellabel, show90label), wrap, [], [], [], [], true);
            %             ht = title(sprintf('%s', datestr(dates(thisT), 'yyyyqq')));
            %             wrapthisfigure(thisfig, sprintf('Ypredictivedensity%s-%s%s-WITHTITLE', datestr(dates(thisT), 'yyyyqq'), ...
            %                 modellabel, show90label), wrap, [], [], [], [], true);

            % add realized data
            if ~all(isnan(Yfuture(thisT,:)))
                ylim('auto')
                hFuture = plot(0:Nhorizons-1, Yfuture(thisT,:), '-', 'color', dataColor, ...
                    'linewidth', 5);
                ylimits = ylim;
                switch upper(datalabel)
                    case {'CPI'}
                        if all(ylimits > 0)
                            ylimits(1) = 0;
                        elseif all(ylimits < 0)
                            ylimits(2) = 0;
                        end
                        if dates(thisT) >= datenum(2020,1,1)
                            ylimits = [-4 10];
                        end
                    case 'UNRATE'
                        switch year(dates(thisT))
                            case {2020,2021,2022}
                                ylimits = ylim;
                                if ylimits(1) >= 0
                                    ylimits(1) = 0;
                                end
                                if ylimits(2) <= 15
                                    ylimits(2) = 15;
                                end
                            otherwise
                                ylimits = ylim;
                        end
                    case 'RGDP'
                        switch year(dates(thisT))
                            case 2020
                                ylimits = [-50 40];
                            case {2021,2022}
                                ylimits = [-15 20];
                            otherwise
                                ylimits = ylim;
                        end
                end
                ylim(ylimits)
                wrapthisfigure(thisfig, sprintf('Ypredictivedensity%s-%s%s-WITHDATA', datestr(dates(thisT), 'yyyyqq'), ...
                    modellabel, show90label), wrap, [], [], [], [], true);
                hl = legend([hsim hFuture], 'Predictive Density', 'Realized Data', ...
                    'location', 'best');
                wrapthisfigure(thisfig, sprintf('Ypredictivedensity%s-%s%s-WITHLEGENDDATA', datestr(dates(thisT), 'yyyyqq'), ...
                    modellabel, show90label), wrap, [], [], [], [], false);
                delete(hl)
            end

            if ispc
                close(thisfig);
            end

        end % thisT

        finishwrap

    end % datalabel
end % modeltype

%% finish / clean up
finishscript
dockAllFigures
