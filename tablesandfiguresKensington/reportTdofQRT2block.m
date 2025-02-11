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
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'};
fontsize   = 18;


resultdir = '../mcmcKensington/foo/';

%% define list of models
NGAP = {'Ny', 'BOP'};
NGAP = {'BOP'};
MODELTYPES = {'VAR0trendHScycleSVt2blockNoiseHS'};
SAMSTART = {'1968Q4'};
SPFquarterly       = {'-EXy1q4', '-y1q4'};
prettySPFquarterly = {'(excl Q4 of y1)','(incl Q4 of y1)'};

MODELTYPES         = {'VAR0trendHScycleSVt2blockNoiseHS', 'CGpriorVAR0trendHScycleSVt2blockNoiseHS'};
SPFquarterly       = {'-EXy1q4', '-y1q4'};
prettySPFquarterly = {'(excl Q4 of y1)','(incl Q4 of y1)'};

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

%% initwrap
initwrap

%% loop over models
for d = 1 : length(DATALABELS)

    datalabel    = DATALABELS{d};
    for m0 = 1 : numel(models)

        modeltype0   = models(m0).type;
        modelpretty0 = models(m0).pretty;
        Ndraws0      = models(m0).Ndraws;

        % %% latexwrapper (per modeltype)
        % wrap   = [];
        % titlename = sprintf('reportTdof-%s', modeltype0);
        % if isempty(wrap) && ~isdesktop
        %     initwrap
        % end
        %
        %% loop over datalabel

        %% load data
        modellabel0  = strcat(datalabel, '-', modeltype0);
        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
        fullmatfilename = fullfile(resultdir, matfilename);
        if ~exist(fullmatfilename, 'file')
            warning('%s not found', matfilename);
            continue
        end
        mat0 = matfile(fullmatfilename);

        dates  = mat0.dates;
        datesQ = mat0.datesQ;
        tDof   = mat0.tDOFdraws;
        Tstart = mat0.Tstart;

        if ismatrix(tDof)
            error('does not look like 2block output')
        end
        tDof = permute(tDof, [2 1 3]); % T x Nsv x Ndraws
        Nsv  = size(tDof,2);

        %% collect tDof median and 68% interval
        mid  = median(tDof, 3);
        tail = prctile(tDof, normcdf11, 3);

        %% cut sample
        ndx    = Tstart : length(dates);
        dates  = dates(ndx);
        datesQ = datesQ(ndx);
        mid    = mid(ndx,:);
        tail   = tail(ndx,:,:);

        %% plot tDof moments over time
        thiscolor = colors4plots("darkblue");
        thisfig = figure;
        for ii = 1 : Nsv
            subplot(Nsv,1,ii)
            hold on
            hmid = plot(dates, mid(:,ii), '-', 'color', thiscolor, 'linewidth', 2);
            htail = plot(dates, squeeze(tail(:,ii,:)), '-', 'color', thiscolor, 'linewidth', 1);
            hQ1 = xline(dates(datesQ==1), 'k:', 'LineWidth', .5);
            yline(3, ':')
            yline(40, ':')
            xtickdates(dates)
            ylim([0 45])
            yticks([0 : 5 : 45])
            % legend([hmid, htail(1) hQ1(1)], 'median', '68% interval', 'Q1 dates', 'location', 'best')
            title(sprintf('SV Block %d', ii))
        end
        sgtitle(sprintf('tDof %s %s', datalabel, modelpretty0))
        wrapthisfigure(thisfig, sprintf('tDof-%s-%s', datalabel, modelpretty0), wrap);

    end % datalabel

    %% finish latexwrapper
    if ~isempty(wrap)
        close all
    end
end % models

%% finish
dockAllFigures
finishwrap
finishscript