%% collect and plot moments for CG slopes from QRT runs

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
%#ok<*SAGROW>
%#ok<*DEFNU>
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc


%% parameters
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'}; % , 'TBILLcensored'};
fontsize   = 18;


resultdir = '~/jam/lager/kensingtonStateSpaceDraws/';

%% define list of models
% NGAP = {'Ny', 'BOP'};
NGAP = {'BOP'};

MODELTYPES         = {'VAR0trendHScycleSVt','VAR0trendHScycleSVtnoise','VAR0trendHScycleSVtnoiseHS'};
SAMSTART           = {'1968Q4'};
SPFquarterly       = {'', '-y1q4', '-SPFquarterlyOnly'};
prettySPFquarterly = {'(excl Q4 of y1)','(incl Q4 of y1)', ' (fixed h only)'};

MODELTYPES         = {'VAR0trendHScycleSVtnoiseHS','VAR0trendHScycleSVt2blockNoiseHS'};
SPFquarterly       = {'-EXy1q4', '-y1q4', '-SPFquarterlyOnly'};
prettySPFquarterly = {'(excl Q4 of y1)','(incl Q4 of y1)', ' (fixed h only)'};

resultdir = '~/prod/kensington-buba/mcmcKensington/matfiles';
NGAP = {'BOP'};
MODELTYPES         = {'VAR0trendHScycleSVt2blockNoiseHS'};
SPFquarterly       = {'-y1q4'};
prettySPFquarterly = {'(incl Q4 of y1)'};

for nn = 1 : length(NGAP)
    for mm = 1 : length(MODELTYPES)
        for ss = 1 : length(SAMSTART)
            for qq = 1 : length(SPFquarterly)
                samStartLabel = SAMSTART{ss};
                models(nn,mm,ss,qq).type   = sprintf('%s%s-Ngap%s-samStart%s', MODELTYPES{mm}, SPFquarterly{qq}, NGAP{nn}, samStartLabel);
                models(nn,mm,ss,qq).Ngap   = NGAP{nn};
                models(nn,mm,ss,qq).pretty = sprintf('%s-%s%s', MODELTYPES{mm}, NGAP{nn}, prettySPFquarterly{qq});
                models(nn,mm,ss,qq).Ndraws = 3e3;
            end
        end
    end
end

%% loop over models
for m0 = 1 : numel(models)

    modeltype0   = models(m0).type;
    modelpretty0 = models(m0).pretty;
    Ndraws0      = models(m0).Ndraws;

    %% latexwrapper (per modeltype)
    wrap   = [];
    titlename = sprintf('reportCGslopes-%s', modeltype0);
    initwrap
    if isempty(wrap) && ~isdesktop
        initwrap
    end

    %% loop over datalabel
    for d = 1 : length(DATALABELS)

        datalabel    = DATALABELS{d};

        %% load data
        modellabel0  = strcat(datalabel, '-', modeltype0);
        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
        fullmatfilename = fullfile(resultdir, matfilename);
        if ~exist(fullmatfilename, 'file')
            warning('%s not found', matfilename);
            continue
        end
        mat0 = matfile(fullmatfilename);

        dates    = mat0.dates;
        datesQ   = mat0.datesQ;
        CGdraws  = permute(mat0.CGPOOLEDdraws, [3 1 2]); % Ndates x Npools x Ndraws
        CGpools  = mat0.CGpools;
        Tstart   = mat0.Tstart;

        ndx  = Tstart : length(dates);

        %% first, plot Glambda
        thesedates     = dates(ndx);
        thesedatesQ    = datesQ(ndx);

        thiscolor = colors4plots("darkred");
        thatcolor = colors4plots("black");
        thisfig = figure;
        hold on
        hmid = plot(thesedates, mat0.GlambdaMedian(:,ndx), '-', 'color', thiscolor, 'linewidth', 2);
        htail = plot(thesedates, mat0.GlambdaTails(:,ndx), '-', 'color', thiscolor, 'linewidth', 1);
        hprior = yline(mat0.GpriorlambdaMedian(:,ndx(1)), '-', 'color', thatcolor, 'linewidth', 2);
        yline(mat0.GpriorlambdaTails(:,ndx(1)), '-', 'color', thatcolor, 'linewidth', 2);
        % plot(dates(ndx), mat0.GpriorlambdaMedian(:,ndx), '-', 'color', thatcolor, 'linewidth', 2);
        % plot(dates(ndx), mat0.GpriorlambdaTails(:,ndx), '-', 'color', thatcolor, 'linewidth', 1);
        hQ1 = xline(thesedates(thesedatesQ==1), 'k:', 'LineWidth', .5);
        legend([hmid, htail(1) hprior hQ1(1)], 'posterior median', 'posterior IQR', 'prior', 'Q1 dates', 'location', 'best');
        ylim([0 1])
        xtickdates(dates(ndx))
        title(sprintf('Max root of VAR\n %s %s', datalabel, modelpretty0))
        wrapthisfigure(thisfig, sprintf('Glambda-%s-%s', datalabel, modeltype0), wrap);

        %% collect moments
        mid    = median(CGdraws, 3);
        tail68 = prctile(CGdraws, normcdf11, 3);
        tail90 = prctile(CGdraws, [5 95], 3);

        %% loop over pools
        thiscolor = colors4plots("darkblue");
        for ii = 1 : length(CGpools)
            thispool = CGpools{ii};
            % cut sample
            thesedates     = dates(ndx);
            thesedatesQ    = datesQ(ndx);
            thismid   = mid(ndx,ii);
            thistail68  = squeeze(tail68(ndx,ii,:));
            thistail90  = squeeze(tail90(ndx,ii,:));
            % plot tDof moments over time
            thisfig = figure;
            hold on
            hmid = plot(thesedates, thismid, '-', 'color', thiscolor, 'linewidth', 2);
            htail68 = plot(thesedates, thistail68, '-', 'color', thiscolor, 'linewidth', 1);
            htail90 = plot(thesedates, thistail90, '--', 'color', thiscolor, 'linewidth', 1);
            hQ1 = xline(thesedates(thesedatesQ==1), 'k:', 'LineWidth', .5);
            xtickdates(thesedates)
            legend([hmid, htail68(1) htail90(1) hQ1(1)], 'median', '68% band', '90% band', 'Q1 dates', 'location', 'best')
            title(sprintf('CG pooled (h=%d:%d)\n%s %s', thispool(1), thispool(end), datalabel, modelpretty0))
            wrapthisfigure(thisfig, sprintf('CGPOOL%d-%s-%s', ii, datalabel, modeltype0), wrap);
        end % pools
    end % datalabel

    %% finish latexwrapper
    finishwrap
    if ~isempty(wrap)
        close all
    end
end % models

%% finish
dockAllFigures
finishscript