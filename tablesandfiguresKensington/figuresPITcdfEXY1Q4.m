%% collect and compare PITs "CDF" style

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

fontsize = 18;

doTitle = false;

samStartLabel = '1968Q4';
%% set modeldir
resultdir = '~/jam/lager/kensingtonStateSpaceDraws/';



%%  define models
% NGAP               = {'Ny', 'BOP'};
NGAP               = {'BOP'};
samStartLabel      = '1968Q4';

MDStype            = {'MDS','VAR0'};
MDSpretty          = {'MDS', 'VAR'};
SHOCKtype          = {'trendHScycleSVt2block', 'trendHScycleSVt2blockNoiseHS'};
SPFquarterly       = {'-EXy1q4', '-y1q4'};
prettySPFquarterly = {'w/o Noise','w/Noise'};

% MDStype            = {'MDS', 'VAR0'};
% MDSpretty          = {'MDS', 'VAR'};
% SHOCKtype          = {'trendHScycleSVt2blockNoiseHS'};
% SPFquarterly       = {'-y1q4'};
% prettySPFquarterly = {''};

for nn = 1 : length(NGAP)
    for mm = 1 : length(MDStype)
        
        for ss = 1 : length(SHOCKtype)
            for qq = 1 : length(SPFquarterly)
                models(nn,mm,ss,qq).type      = sprintf('%s%s%s-Ngap%s-samStart%s', MDStype{mm}, SHOCKtype{ss}, SPFquarterly{qq}, NGAP{nn}, samStartLabel);
                models(nn,mm,ss,qq).MDStype   = MDStype{mm};
                models(nn,mm,ss,qq).SHOCKtype = SHOCKtype{ss};
                models(nn,mm,ss,qq).Ngap      = NGAP{nn};
                models(nn,mm,ss,qq).datatype  = prettySPFquarterly{qq};
                models(nn,mm,ss,qq).pretty    = prettySPFquarterly{qq}; % sprintf('%s', MDSpretty{mm});
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
for mm = 1 : length(MDStype)
    ndx = NaN(length(SPFquarterly),1);
    for qq = 1 : length(SPFquarterly)
        ss = qq;
        ndx(ss) = sub2ind(size(models),nn,mm,ss,qq);
    end
    groups = cat(1, groups, {ndx(:)});
    grouplabels = cat(1, grouplabels, MDStype{mm});
    groupprettylabels = cat(1, groupprettylabels, sprintf('Model comparison (%s)', MDSpretty{mm}));
    % group0label = cat(1, group0label, SHOCKtype{1});
end % mm

%% define set of eval windows
s = 0;

% s = s + 1;
% sam(s).start  = [];
% sam(s).stop   = [];
% sam(s).label  = 'fullsample';
% sam(s).pretty = 'full sample';

s = s + 1;
sam(s).start  = datenum(1990,1,1);
sam(s).stop   = datenum(2023,10,1);
sam(s).label  = 'fullsampleSince1990';
sam(s).pretty = 'full sample (since 1990)';

% s = s + 1;
% sam(s).start  = [];
% sam(s).stop   = datenum(2019,10,1);
% sam(s).label  = 'preCOVID';
% sam(s).pretty = 'pre COVID (2019)';

s = s + 1;
sam(s).start  = datenum(1990,1,1);
sam(s).stop   = datenum(2019,10,1);
sam(s).label  = 'since1990preCOVID';
sam(s).pretty = '1990 -- 2019';
%
% s = s + 1;
% sam(s).start  = [];
% sam(s).stop   = datenum(2008,10,1);
% sam(s).label  = 'preGFC';
% sam(s).pretty = 'pre GFC (until 2008)';
%
% s = s + 1;
% sam(s).start  = datenum(1990,1,1);
% sam(s).stop   = datenum(2008,10,1);
% sam(s).label  = '1990to2008';
% sam(s).pretty = '1990 -- 2008';
%
%
% s = s + 1;
% sam(s).start  = datenum(2009,1,1);
% sam(s).stop   = datenum(2019,10,1);
% sam(s).label  = 'sinceGFCpreCOVID';
% sam(s).pretty = 'GFC until COVID (2009-2019)';
%
% s = s + 1;
% sam(s).start  = datenum(2009,1,1);
% sam(s).stop   = [];
% sam(s).label  = 'sinceGFC';
% sam(s).pretty = 'since GFC (as of 2009)';

%% for now, purge all samples starting prior to 1990
ndx = arrayfun(@(x) isempty(x.start) || x.start < datenum(1990,1,1), sam);
sam = sam(~ndx);

%% loop over datalabels

for d = 1 : length(DATALABELS)

    datalabel = DATALABELS{d};


    %% loop over samples
    for s = 1 : length(sam)

        %% prepare latexwrapper
        wrap = [];
        titlename = strcat('PITcdf-', datalabel, '-', sam(s).label);
        initwrap


        %% close figures
        if ~isempty(wrap)
            close all
        end

        %% probe for sample length (Based on m0)
        m = 1;
        modellabel  = strcat(datalabel, '-',  models(m).type);
        matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel, Ndraws);
        fullmatfilename = fullfile(resultdir, matfilename);
        mat             = matfile(fullmatfilename);

        dates         = mat.dates;
        T             = mat.T;
        Nhorizons     = mat.Nhorizons;
        if isempty(sam(s).start)
            samStart = Tstart;
        else
            samStart = find(dates >= sam(s).start, 1, 'first');
        end
        if isempty(sam(s).stop)
            samStop = T;
        else
            samStop = find(dates <= sam(s).stop, 1, 'last');
        end
        sammy = (dates >= dates(samStart)) & (dates <= dates(samStop));
        Tsam  = sum(sammy);



        %% loop over all models
        fcstYpits = NaN(Tsam,Nhorizons,numel(models));
        for m = 1 : numel(models)

            modellabel = strcat(datalabel, '-',  models(m).type);

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


            Nz            = mat.Nz;
            Nhorizons     = mat.Nhorizons;
            Nbar          = mat.Nbar;

            dates         = mat.dates;
            T             = mat.T;
            datesQ        = mat.datesQ;
            doNIPA        = mat.doNIPA;

            Tstart        = mat.Tstart;


            fprintf('Processing %s ... \n', modellabel)

            %% cut sample
            if isempty(sam(s).start)
                samStart = Tstart;
            else
                samStart = find(dates >= sam(s).start, 1, 'first');
            end

            if isempty(sam(s).stop)
                samStop = T;
            else
                samStop = find(dates <= sam(s).stop, 1, 'last');
            end

            sammy = find((dates >= dates(samStart)) & (dates <= dates(samStop))); % find needed for matfile access

            fcstYpits(:,:,m) = mat.fcstYpits(sammy,:);


        end % model

        fcstYpits = permute(fcstYpits, [1 3 2]);

        % go over models and horizons, and compute bootstrapped KS and CvM p-values and KS band half width
        [KSp, CvMp, KSBandHalf]  = deal(NaN(numel(models) , Nhorizons));

        for m = 1 : size(KSp, 1)
            if ~modelsNaN(m)
                for hz = 1 : size(KSp, 2)

                    thispit = fcstYpits(:, m,hz);
                    [~,tmpKSp,~,tmpCvMp,tmpKS95cv,tmpTpit,~] = PITtest(thispit);

                    KSp(m,hz)        = tmpKSp;
                    CvMp(m,hz)       = tmpCvMp;
                    KSBandHalf(m,hz) = tmpKS95cv / sqrt(tmpTpit);

                end
            end
        end

        

        %% plot different groups

        for g =  1 : numel(groups)
            thesemodels = groups{g};

            if ~all(modelsNaN(thesemodels))
                %% single panel figures
                for h = 1 : Nhorizons
                    if doTitle
                        thistitle = sprintf('%s, h=%d\n%s', datalabel, h-1, sam(s).pretty);
                    else 
                        thistitle = [];
                    end
                    [thisfig, hl] = plotPITcdfstairs(fcstYpits(:,thesemodels,h), KSBandHalf(thesemodels,h), {models(thesemodels).pretty}, thistitle);
                    wrapthisfigure(thisfig, sprintf('PITs-%s-h%d-%s-EXY1Q4-%s-WITHLEGEND', datalabel, h-1, sam(s).label, grouplabels{g}), wrap);
                    delete(hl)
                    wrapthisfigure(thisfig, sprintf('PITs-%s-h%d-%s-EXY1Q4-%s', datalabel, h-1, sam(s).label, grouplabels{g}), wrap);                    
                end
            end
            
            if ~isempty(wrap)
                close all
            end
        end

        %% compile latexwrapper per sample
        finishwrap

    end % sample
end % datalabel

%% finish / clean up
finishscript
dockAllFigures

%% local functions
function [thisfig, hl] = plotPITcdfstairs(thesepits, thesehalfbands, theselabels, thistitle)

fontsize   = 18;
Nmodels    = length(theselabels);
rvec    = 0 : 0.001 : 1;
legHandles = [];
legLabels  = cell(0);
thisfig    = figure;
set(gca, 'FontSize', fontsize);
hold on;
for m = 1 : Nmodels
    thispit   = thesepits(:,m);
    if ~all(isnan(thispit))
        [f,x]     = ecdf(thispit);
        this      = stairs(x,f, 'LineWidth', 3, 'Color', colors4plots(m));
        thishalfband = thesehalfbands(m);
        legHandles = cat(1, legHandles, this);
        legLabels  = cat(1, legLabels, theselabels(m));
        if  ~isodd(m) % contains(theselabels(m), 'VAR')
            this.LineStyle = '-.';
        end
        % get bands
        [~,~,~,~,KS95cv,Tpit,rvec] = PITtest(thispit);
        bplus  = plot(rvec,rvec + thishalfband,'Color',this.Color,'LineStyle','-','LineWidth',1);
        bminus = plot(rvec,rvec - thishalfband,'Color',this.Color,'LineStyle','-','LineWidth',1);
        if ~isodd(m) % contains(theselabels(m), 'VAR')
            bplus.LineStyle = '-.';
            bminus.LineStyle = '-.';
        end
        legHandles = cat(1, legHandles, bplus);
        legLabels  = cat(1, legLabels, sprintf('Confidence band (%s)', theselabels{m}));
    end
end
plotfortyfive('k:', 'linewidth', 2);
xlim([0 1]); ylim([0 1]);
% legHandles = legHandles([1 4 2]);
% legLabels  = legLabels([1 3 2]);

hl = legend(legHandles, legLabels, 'location', 'northwest');
if ~isempty(thistitle)
    title(thistitle)
end

end % function plotPITcdfstairs



function [Kvstat,Kpval,CVMvstat,CvMpval,KS95cv,T,rvec] = PITtest(PITvec)
%% PIT test of Rossi and Sekhposyan (2019)

% get rid of NaNs - input can contain NaNs but the valid observations must be contiguous
PITvec = PITvec(~isnan(PITvec));
T      = size(PITvec, 1); % PITvec is (T x 1)
el     = floor(T^(1/3)); % block length, must be o(T^0.5), RS used T^1/3 and T^1/4 in MC
bootMC = 1000;
rvec   = 0 : 0.001 : 1; % rvec is (1 x numr)
numr   = size(rvec,2);

cumcumz = (PITvec < rvec) - rvec;
v       = sum(cumcumz,1)/sqrt(T);

Kvstat   = max(abs(v));
CVMvstat = mean((v.^2));
pgrid    = 0.01 : 0.01 :0.99;
tablecv  = bootstrapInoue(el,bootMC,PITvec,rvec,pgrid);
% 95th percentile of the bootstrap distribution of the Kolmogorov-Smirnov statistic, for plotting
KS95cv=tablecv(1,95);
% now look up the p-value
% KS stat
Kpval=1-interp1(tablecv(1,:)',pgrid',Kvstat,'linear','extrap');
if Kpval < 0
    Kpval = 0;
elseif Kpval > 1
    Kpval = 1;
end
% CvM stat
CvMpval=1-interp1(tablecv(2,:)',pgrid',CVMvstat,'linear','extrap');
if CvMpval < 0
    CvMpval = 0;
elseif CvMpval > 1
    CvMpval = 1;
end

end % function PITtest

function result = bootstrapInoue(el, bootMC, pit, rvec, pgrid)
%% bootstrapping PIT stats
KSv  = zeros(bootMC,1);
CVMv = zeros(bootMC,1);

P         = size(pit,1);
invsqrtP  = 1 ./ sqrt(P);

emp_cdf        = pit <= rvec;
emp_cdf        = emp_cdf - mean(emp_cdf,1);
emp_cdf_sumel  = NaN(size(emp_cdf));
for j = 1 : P-el+1
    emp_cdf_sumel(j,:) = sum(emp_cdf(j:j+el-1,:), 1);
end
bigNormMat = 1/sqrt(el) * randn(P-el+1,bootMC);

parfor bootrep = 1 : bootMC

    z      = bigNormMat(:, bootrep);
    K_star = zeros(1,length(rvec));
    
    for j = 1 : P-el+1
         K_star = K_star + invsqrtP .* z(j,1) .* emp_cdf_sumel(j,:); %#ok<PFBNS>
    end

    KSv(bootrep,1) =  max(abs(K_star'));
    CVMv(bootrep,1) = mean(K_star'.^2);
end

ndx  = round(bootMC*pgrid);
KSv  = sort(KSv,'ascend');
cvKv = KSv(ndx);
CVMv = sort(CVMv,'ascend');
cvMv = CVMv(ndx);

result = [cvKv'; cvMv'];
end % bootstrapInoue
