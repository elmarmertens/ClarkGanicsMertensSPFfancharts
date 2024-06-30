%% collect and compare logscores

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

Ndraws   = 3e3;
Npitbins = 10;

fontsize = 24;


samStartLabel = '1968Q4';
%% set modeldir
resultdir = '~/jam/lager/kensingtonStateSpaceDraws/';



%% define models
NGAP               = {'BOP'};
samStartLabel      = '1968Q4';
MDStype            = {'MDS', 'VAR0'};
MDSpretty          = {'MDS', 'VAR'};
SHOCKtype          = {'trendHScycleSVt2blockNoiseHS'};
SPFquarterly       = {'-y1q4'};
prettySPFquarterly = {'incl Q4 of y1'};

for nn = 1 : length(NGAP)
    for mm = 1 : length(MDStype)
        for ss = 1 : length(SHOCKtype)
            for qq = 1 : length(SPFquarterly)
                models(nn,mm,ss,qq).type      = sprintf('%s%s%s-Ngap%s-samStart%s', MDStype{mm}, SHOCKtype{ss}, SPFquarterly{qq}, NGAP{nn}, samStartLabel);
                models(nn,mm,ss,qq).MDStype   = MDStype{mm};
                models(nn,mm,ss,qq).SHOCKtype = SHOCKtype{ss};
                models(nn,mm,ss,qq).Ngap      = NGAP{nn};
                models(nn,mm,ss,qq).datatype  = prettySPFquarterly{qq};
                models(nn,mm,ss,qq).pretty    = MDSpretty{mm}; % sprintf('%s-%s (%s)', MDStype{mm}, SHOCKtype{ss}, prettySPFquarterly{qq});
                models(nn,mm,ss,qq).Ndraws    = 3e3;
            end
        end
    end
end

modelsNaN = true(size(models));

%% define modelgroups for comparison
groups      = cell(0);
grouplabels = cell(0);
groupprettylabels = cell(0);


% group by modeltype (and data input) and compare across shock types
nn = 1; % BOP only
for mm = 1 : 2 : length(MDStype)
    for qq = 1 : length(SPFquarterly)
        ndx = NaN(length(SHOCKtype),2);
        for ss = 1 : length(SHOCKtype)
            ndx(ss,1) = sub2ind(size(models),nn,mm,ss,qq);
            ndx(ss,2) = sub2ind(size(models),nn,mm+1,ss,qq);
        end
        % modeltype = strcat(MDStype{mm}, SPFquarterly{qq});
        groups = cat(1, groups, {ndx(:)});
        if isempty(SPFquarterly{qq})
            grouplabels = cat(1, grouplabels, '-EXy1q4');
        else
            grouplabels = cat(1, grouplabels, SPFquarterly{qq});
        end
        groupprettylabels = cat(1, groupprettylabels, sprintf('Model comparison (%s)', prettySPFquarterly{qq}));
        % group0label = cat(1, group0label, SHOCKtype{1});
    end % ss
end % mm


%% loop over datalabels
wrap = [];
titlename = 'ZlogscoresACROSSmodels';
initwrap

for d = 1 : length(DATALABELS)

    datalabel = DATALABELS{d};
    close all

    %% probe for sample length (Based on m0) and allocate memory
    m = 1;
    modellabel  = strcat(datalabel, '-',  models(m).type);
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel, Ndraws);
    fullmatfilename = fullfile(resultdir, matfilename);
    mat             = matfile(fullmatfilename);

    dates         = mat.dates;
    T             = mat.T;
    Tstart        = mat.Tstart;

    fcstZlogscores = NaN(T,numel(models));


    %% loop over all models

    for m = 1 : numel(models)

        modellabel = strcat(datalabel, '-',  models(m).type);
        fprintf('Processing %s ... \n', modellabel)

        %% load data

        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel, Ndraws);
        fullmatfilename = fullfile(resultdir, matfilename);
        if ~exist(fullmatfilename, 'file')
            modelsNaN(m) = true;
            warning('%s not found', matfilename);
            continue
        else
            modelsNaN(m) = false;
        end
        mat = matfile(fullmatfilename);

        fcstZlogscores(:,m) = mat.fcstZmvlogscore; % 2D indexing needed by matfile


    end % model


    %% process logscores per group
    for gg = 1 : numel(groups)
        thesemodels = groups{gg};
        % if modelsNaN(thesemodels(1)) && ~isempty(thesemodels(~modelsNaN(thesemodels)))
        %     warning('benchmark model in group <<%s>> contains a NaN element', grouplabels{gg})
        %     continue
        % end
        thesemodels = thesemodels(~modelsNaN(thesemodels));
        if numel(thesemodels) < 2
            warning('group <<%s>> contains fewer than 2 elements', grouplabels{gg})
            continue
        end
        thesescores = fcstZlogscores(:,thesemodels);

        sammy          = Tstart:T; % for now: full sample plots
        thesedates     = dates(sammy);
        cumscores      = cumsum(thesescores(sammy,:),1);
        relativescores = cumscores(:,2:end) - cumscores(:,1);
        % take avg scores
        relativescores = relativescores ./ transpose(1 : length(thesedates));

        % plot relativescores
        thisfig = figure;
        hZ = plot(thesedates, relativescores, 'linewidth', 2);
        set(gca, 'fontsize', fontsize)
        yline(0, 'k:', 'linewidth', 2)
        nbershades(thesedates)
        box(gca, 'off')
        % title(sprintf('%s relative to %s', models(thesemodels(2)).pretty, models(thesemodels(1)).pretty)); 
        wrapthisfigure(thisfig, sprintf('avglogscorediffs%s-%s', grouplabels{gg}, datalabel),wrap)

    end
    %% finish latexwrapper (per datalabel)
    % finishwrap
end % datalabel

%% finish / clean up
finishwrap
finishscript
dockAllFigures


