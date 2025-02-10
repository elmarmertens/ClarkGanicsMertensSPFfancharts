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

theseDates = datenum([2019 2023], 10, 1);

resultdir = '../mcmcKensington/foo/';

%% define list of models
NGAP = {'BOP'};
SAMSTART           = {'1968Q4'};
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


%% allocate memory
m0 = 1;
datalabel    = DATALABELS{1};
modeltype0   = models(m0).type;
modelpretty0 = models(m0).pretty;
Ndraws0      = models(m0).Ndraws;

% preload a few things for dd=1
modellabel0  = strcat(datalabel, '-', modeltype0);
matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
fullmatfilename = fullfile(resultdir, matfilename);
mat0 = matfile(fullmatfilename);
dates    = mat0.dates;
datesQ   = mat0.datesQ;
Tstart   = mat0.Tstart;
T        = length(dates);

Ndata = length(DATALABELS);
CGslopesMid    = NaN(T, Ndata);
CGslopesTail68 = NaN(T, 2, Ndata);
CGslopesTail90 = NaN(T, 2, Ndata);


%% loop over datalabel
for dd = 1 : length(DATALABELS)
    
    datalabel    = DATALABELS{dd};
    
    %% load data
    modellabel0  = strcat(datalabel, '-', modeltype0);
    matfilename = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel0, Ndraws0);
    fullmatfilename = fullfile(resultdir, matfilename);
    if ~exist(fullmatfilename, 'file')
        warning('%s not found', matfilename);
        continue
    end
    mat0 = matfile(fullmatfilename);
    
    CGdraws  = permute(mat0.CGPOOLEDdraws, [3 2 1]); % Ndates x Ndraws x Npools with Npools = 1
    
    %% collect moments
    CGslopesMid(:,dd)     = median(CGdraws, 2);
    CGslopesTail68(:,:,dd) = prctile(CGdraws, normcdf11, 2);
    CGslopesTail90(:,:,dd) = prctile(CGdraws, [5 95], 2);
    
end % datalabel

%% select results 
ndxT = ismember(dates, theseDates);

Ndates    = length(theseDates);
Nstats    = 3; % mid and 90% tail
slopesMid = CGslopesMid(ndxT,:);
slopes68  = CGslopesTail68(ndxT,:,:);
slopes90  = CGslopesTail90(ndxT,:,:);

%% tabulate values and output in TeX
initwrap
tabname    = sprintf('table-CGslopes-%s.tex', modeltype0);
tabcaption = sprintf('Slopes of Coibion-Gorodnichenko Regressions');
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
else
    tabdir = fullfile(localtemp, 'foo');
end
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, Nstats * Ndates));
fprintf(fid, '\\toprule\n');
% header: dates row
for tt = 1 : Ndates
    fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, datestr(theseDates(tt), 'yyyyqq'));
end
fprintf(fid, '\\\\\n');
offset = 1;
for dd = 1 : Ndates
    fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats);
    offset = offset + Nstats;
end
fprintf(fid, '\n');
% header stats row
fprintf(fid, ' Variable ');
for dd = 1 : Ndates
    fprintf(fid, ' & \\multicolumn{1}{c}{5\\%%} & \\multicolumn{1}{c}{50\\%%} & \\multicolumn{1}{c}{95\\%%} ');
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for dd = 1 : Ndata
    fprintf(fid, '%s ', DATALABELS{dd});
    for tt = 1 : Ndates
        fprintf(fid, ' & %6.2f & %6.2f & %6.2f ', slopes90(tt,1,dd), slopesMid(tt,dd), slopes90(tt,2,dd));
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fclose(fid);
type(fullfile(tabdir, tabname))

%% finish
dockAllFigures
finishwrap
finishscript