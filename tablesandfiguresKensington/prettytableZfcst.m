%% collect predicted SPF data (Z)

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
% addpath(genpath('~/matlab/library/jplv7'), '-end')

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

resultdir = '../mcmcKensington/foo/';
datadir   = fullfile('..', 'matdataKensington');

for doEXY1Q4 = [false true]

    %% define models
    SAMSTART           = {'1968Q4'};
    NGAP               = {'BOP'};
    MODELTYPESpretty   = {'MDS','VAR'};

    if ~doEXY1Q4
        MODELTYPES         = {'MDStrendHScycleSVt2blockNoiseHS','VAR0trendHScycleSVt2blockNoiseHS'};
        SPFquarterly       = {'-y1q4'};
        prettySPFquarterly = {'(incl Q4 of y1)'};
    else
        MODELTYPES         = {'MDStrendHScycleSVt2block','VAR0trendHScycleSVt2block'};
        SPFquarterly       = {'-EXy1q4'};
        prettySPFquarterly = {'(excl Q4 of y1)'};
    end

    for nn = 1 : length(NGAP)
        for mm = 1 : length(MODELTYPES)
            for ss = 1 : length(SAMSTART)
                for qq = 1 : length(SPFquarterly)
                    samStartLabel = SAMSTART{ss};
                    models(nn,mm,ss,qq).type   = sprintf('%s%s-Ngap%s-samStart%s', MODELTYPES{mm}, SPFquarterly{qq}, NGAP{nn}, samStartLabel);
                    models(nn,mm,ss,qq).Ngap   = NGAP{nn};
                    models(nn,mm,ss,qq).pretty = MODELTYPESpretty{mm};
                    models(nn,mm,ss,qq).Ndraws = 3e3;
                end
            end
        end
    end


    %% define set of eval windows
    s = 0;

    % s = s + 1;
    % sam(s).start  = [];
    % sam(s).stop   = [];
    % sam(s).label  = 'fullsample';
    % sam(s).pretty = 'full sample';

    s = s + 1;
    sam(s).start  = datenum(1990,1,1);
    sam(s).stop   = [];
    sam(s).label  = 'fullsampleSince1990';
    sam(s).pretty = 'full sample (since 1990)';

    s = s + 1;
    sam(s).start  = datenum(1990,1,1);
    sam(s).stop   = datenum(2019,10,1);
    sam(s).label  = 'since1990preCOVID';
    sam(s).pretty = '1990 -- 2019';

    Nsam = numel(sam);

    %% prepare latexwrapper
    wrap   = [];
    titlename = 'prettytableZfcst';
    if doEXY1Q4
        titlename = strcat(titlename, '-EXY1Q4')
    end
    initwrap
    if isempty(wrap) && ~isdesktop
        initwrap
    end

    %% loop over models and collect results

    maxNz = 9;
    % allocate memory to collect stats for tables
    [MZslopes, MZslopesSE, MZconst, MZconstSE] = deal(NaN(maxNz,Nsam,length(DATALABELS), numel(models)));
    matarray = cell(length(DATALABELS),numel(models));

    for mm = 1 : numel(models)

        modeltype   = models(mm).type;
        modelpretty = models(mm).pretty;
        Ndraws      = models(mm).Ndraws;


        %% loop over datalabel
        for dd = 1 : length(DATALABELS)

            close all
            datalabel    = DATALABELS{dd};


            %% load data
            modellabel   = strcat(datalabel, '-', modeltype);
            matfilename  = sprintf('slimCGMmcmc-%s-Ndraws%d.mat', modellabel, Ndraws);
            fullmatfilename = fullfile(resultdir, matfilename);
            if ~exist(fullmatfilename, 'file')
                % warning('%s not found', matfilename);
                continue
            end
            mat          = matfile(fullfile(resultdir, matfilename));
            matarray{dd,mm} = mat;
            matfilename  = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
            datmat       = matfile(matfilename);

            dates    = mat.dates;
            datesQ   = mat.datesQ;
            Tstart   = mat.Tstart;
            Nz       = mat.Nz;
            ndx      = Tstart : length(dates);
            Zlabel   = datmat.Zlabel;

            Zhat     = mat.fcstZhat;
            Zerror   = mat.fcstZhaterror;
            Ztp1     = Zerror + Zhat;

            %% collect regressions for each sample
            for s = 1 : length(sam)

                % hrulefill
                % fprintf('Sample: %s\n', sam(s).pretty);
                % hrulefill

                thisstart = sam(s).start;
                thisstop  = sam(s).stop;
                thislabel = sam(s).label;
                thispretty = sam(s).pretty;
                if isempty(thisstart)
                    thisstart = dates(1);
                end
                if isempty(thisstop)
                    thisstop = dates(end);
                end
                ndx = dates >= thisstart & dates <= thisstop;
                thisZhat   = Zhat(ndx,:);
                thisZtp1   = Ztp1(ndx,:);
                thisZlabel = Zlabel;
                thisDatesQ = datesQ(ndx);
                thisDates  = dates(ndx);

                % loop over elements of Z
                for z = 1 : Nz

                    ndx       = ~any(isnan([thisZhat(:,z), thisZtp1(:,z)]),2);
                    MZreg     = nwest(thisZtp1(ndx,z), [ones(sum(ndx),1) thisZhat(ndx,z)], 0);

                    % collect stats for table
                    MZslopes(z,s,dd,mm)   = MZreg.beta(2);
                    MZslopesSE(z,s,dd,mm) = sqrt(MZreg.Vbeta(2,2));
                    MZconst(z,s,dd,mm)    = MZreg.beta(1);
                    MZconstSE(z,s,dd,mm)  = sqrt(MZreg.Vbeta(1,1));

                end % z

                clear MZreg

                if ~isempty(wrap)
                    close all
                end
            end % s
        end % datalabel
    end % models


    %% tabulate results in one big table
    % Zlabel  = matarray{1,1}.Zlabel; % assuming it contains maxNz elements
    % if length(Zlabel) ~= maxNz
    %     error('Zlabel mismatch')
    % end
    % for thisSam = 1 : Nsam
    %     Ndata   = length(DATALABELS);
    %     Nmodels = numel(models);
    %     Nstats  = 2;
    %     if ~doEXY1Q4
    %         tabname = sprintf('table-ZhatZtp1-%s.tex', sam(thisSam).label);
    %     else
    %         tabname = sprintf('table-ZhatZtp1-EXY1Q4-%s.tex', sam(thisSam).label);
    %     end
    %     tabcaption = sprintf('Zhat predictability regressions (%s)', sam(thisSam).pretty);
    %     if ~isempty(wrap)
    %         tabdir = wrap.dir;
    %         latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
    %     else
    %         tabdir = fullfile(localtemp, 'foo');
    %     end
    %     fid = fopen(fullfile(tabdir, tabname), 'wt');
    %     fprintf(fid, '\\begin{center}\n');
    %     fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.2', 1, Nstats * Nmodels * Ndata));
    %     fprintf(fid, '\\toprule\n');
    %     % header: datalabel row
    %     for dd = 1 : Ndata
    %         datalabel = DATALABELS{dd};
    %         fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats * Nmodels, datalabel);
    %     end
    %     fprintf(fid, '\\\\');
    %     offset = 1;
    %     for dd = 1 : Ndata
    %         fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats*Nmodels);
    %         offset = offset + Nstats*Nmodels;
    %     end
    %     fprintf(fid, '\n');
    %     % header: models
    %     for dd = 1 : Ndata
    %         for mm = 1 : Nmodels
    %             modeltype = models(mm).pretty;
    %             fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, modeltype);
    %         end
    %     end
    %     fprintf(fid, '\\\\');
    %     offset = 1;
    %     for dd = 1 : Ndata*Nmodels
    %         fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats);
    %         offset = offset + Nstats;
    %     end
    %     fprintf(fid, '\n');
    %     % header const/slopes
    %     fprintf(fid, 'Forecast ');
    %     for s = 1 : Nmodels * Ndata
    %         fprintf(fid, ' & \\multicolumn{1}{c}{const.}  & \\multicolumn{1}{c}{slope} ');
    %     end
    %     fprintf(fid, '\\\\ ');
    %     fprintf(fid, '\\midrule\n');
    %
    %     fprintf(fid, '\\\\');
    %     fprintf(fid, '\\midrule\n');
    %     for zz = 2 : maxNz % skip lagged realized value
    %         fprintf(fid, '%s ', Zlabel{zz});
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     slopestar = abs(MZslopes(zz,thisSam,dd,mm) - 1) > 1.96 * MZslopesSE(zz,thisSam,dd,mm);
    %                     conststar = abs(MZconst(zz,thisSam,dd,mm)) > 1.96 * MZconstSE(zz,thisSam,dd,mm);
    %                     fprintf(fid, ' & %s ', dcolbf(MZconst(zz,thisSam,dd,mm), '%5.2f', conststar));
    %                     fprintf(fid, ' & %s ', dcolbf(MZslopes(zz,thisSam,dd,mm), '%5.2f', slopestar));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\ \n');
    %         % print se
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     fprintf(fid, ' & (%5.2f) & (%5.2f) ', MZconstSE(zz,thisSam,dd,mm), MZslopesSE(zz,thisSam,dd,mm));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\\n');
    %         if zz < maxNz
    %             fprintf(fid, '\\midrule\n');
    %         end
    %     end
    %     fprintf(fid, '\\bottomrule\n');
    %     fprintf(fid, '\\end{tabular}\n');
    %     fprintf(fid, '\\end{center}\n');
    %     fclose(fid);
    %     type(fullfile(tabdir, tabname))
    % end % sam

    %% tabulate results in one big table -- FIRST INTERCEPTS THEN SLOPES
    Zlabel  = matarray{1,1}.Zlabel; % assuming it contains maxNz elements
    if length(Zlabel) ~= maxNz
        error('Zlabel mismatch')
    end
    for thisSam = 1 : Nsam
        Ndata   = length(DATALABELS);
        Nmodels = numel(models);
        Nstats  = 2;
        if ~doEXY1Q4
            tabname = sprintf('table-ZhatZtp1-%s.tex', sam(thisSam).label);
        else
            tabname = sprintf('table-ZhatZtp1-EXY1Q4-%s.tex', sam(thisSam).label);
        end
        tabcaption = sprintf('Zhat predictability regressions (%s)', sam(thisSam).pretty);
        if ~isempty(wrap)
            tabdir = wrap.dir;
            latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
        else
            tabdir = fullfile(localtemp, 'foo');
        end
        fid = fopen(fullfile(tabdir, tabname), 'wt');
        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.2', 1, Nstats * Nmodels * Ndata));
        fprintf(fid, '\\toprule\n');
        % header: intercept vs slopes row
        fprintf(fid, ' & \\multicolumn{%d}{c}{\\bf %s} ', Ndata * Nmodels, 'intercept');
        fprintf(fid, ' & \\multicolumn{%d}{c}{\\bf %s} ', Ndata * Nmodels, 'slope');
        fprintf(fid, '\\\\');
        offset = 1;
        fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Ndata*Nmodels);
        offset = offset + Ndata*Nmodels;
        fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Ndata*Nmodels);
        fprintf(fid, '\n');

        % header: datalabel row
        for dd = 1 : Ndata
            datalabel = DATALABELS{dd};
            fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nmodels, datalabel);
        end
        for dd = 1 : Ndata
            datalabel = DATALABELS{dd};
            fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nmodels, datalabel);
        end
        fprintf(fid, '\\\\');
        offset = 1;
        for dd = 1 : Ndata*2
            fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nmodels);
            offset = offset + Nmodels;
        end
        fprintf(fid, '\n');
        % header: models
        fprintf(fid, 'Forecast ');
        for dd = 1 : Ndata*2
            for mm = 1 : Nmodels
                modeltype = models(mm).pretty;
                fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', 1, modeltype);
            end
        end
        fprintf(fid, '\\\\');
        fprintf(fid, '\\midrule\n');
        for zz = 2 : maxNz % skip lagged realized value
            fprintf(fid, '%s ', Zlabel{zz});
            %% intercepts
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        conststar = abs(MZconst(zz,thisSam,dd,mm)) > 1.96 * MZconstSE(zz,thisSam,dd,mm);
                        fprintf(fid, ' & %s ', dcolbf(MZconst(zz,thisSam,dd,mm), '%5.2f', conststar));
                    else
                        fprintf(fid, ' &  \\multicolumn{1}{c}{---}');
                    end % thisNz
                end % models
            end % datalabel
            %% slopes
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        slopestar = abs(MZslopes(zz,thisSam,dd,mm) - 1) > 1.96 * MZslopesSE(zz,thisSam,dd,mm);
                        fprintf(fid, ' & %s ', dcolbf(MZslopes(zz,thisSam,dd,mm), '%5.2f', slopestar));
                    else
                        fprintf(fid, ' &  \\multicolumn{1}{c}{---}');
                    end % thisNz
                end % models
            end % datalabel
            fprintf(fid, '\\\\ \n');
            % print se: intercepts
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        fprintf(fid, ' & (%5.2f)  ', MZconstSE(zz,thisSam,dd,mm));
                    else
                        fprintf(fid, ' &  ');
                    end % thisNz
                end % models
            end % datalabel
            % print se: slopes
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        fprintf(fid, ' & (%5.2f) ', MZslopesSE(zz,thisSam,dd,mm));
                    else
                        fprintf(fid, ' &  ');
                    end % thisNz
                end % models
            end % datalabel
            fprintf(fid, '\\\\\n');
            if zz < maxNz
                fprintf(fid, '\\midrule\n');
            end
        end
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        fclose(fid);
        type(fullfile(tabdir, tabname))
    end % sam

    %% tabulate results in one big table -- SLOPES ONLY
    Zlabel  = matarray{1,1}.Zlabel; % assuming it contains maxNz elements
    if length(Zlabel) ~= maxNz
        error('Zlabel mismatch')
    end
    for thisSam = 1 : Nsam
        Ndata   = length(DATALABELS);
        Nmodels = numel(models);
        Nstats  = 1;
        if ~doEXY1Q4
            tabname = sprintf('table-ZhatZtp1-slopes-%s.tex', sam(thisSam).label);
        else
            tabname = sprintf('table-ZhatZtp1-slopes-EXY1Q4-%s.tex', sam(thisSam).label);
        end
        tabcaption = sprintf('Zhat predictability regressions (slopes, %s)', sam(thisSam).pretty);
        if ~isempty(wrap)
            tabdir = wrap.dir;
            latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
        else
            tabdir = fullfile(localtemp, 'foo');
        end
        fid = fopen(fullfile(tabdir, tabname), 'wt');
        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.2', 1, Nstats * Nmodels * Ndata));
        fprintf(fid, '\\toprule\n');
        % header: datalabel row
        for dd = 1 : Ndata
            datalabel = DATALABELS{dd};
            fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats * Nmodels, datalabel);
        end
        fprintf(fid, '\\\\');
        offset = 1;
        for dd = 1 : Ndata
            fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats*Nmodels);
            offset = offset + Nstats*Nmodels;
        end
        fprintf(fid, '\n');
        % header: models
        fprintf(fid, 'Forecast ');
        for dd = 1 : Ndata
            for mm = 1 : Nmodels
                modeltype = models(mm).pretty;
                fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, modeltype);
            end
        end
        % fprintf(fid, '\\\\');
        % offset = 1;
        % for dd = 1 : Ndata*Nmodels
        %     fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats);
        %     offset = offset + Nstats;
        % end
        fprintf(fid, '\n');
        % header const/slopes
        % fprintf(fid, 'Forecast ');
        % for s = 1 : Nmodels * Ndata
        %     fprintf(fid, ' & \\multicolumn{1}{c}{const.}  & \\multicolumn{1}{c}{slope} ');
        % end
        % fprintf(fid, '\\\\ ');
        % fprintf(fid, '\\midrule\n');

        fprintf(fid, '\\\\');
        fprintf(fid, '\\midrule\n');
        for zz = 2 : maxNz % skip lagged realized value
            fprintf(fid, '%s ', Zlabel{zz});
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        slopestar = abs(MZslopes(zz,thisSam,dd,mm) - 1) > 1.96 * MZslopesSE(zz,thisSam,dd,mm);
                        % conststar = abs(MZconst(zz,thisSam,dd,mm)) > 1.96 * MZconstSE(zz,thisSam,dd,mm);
                        % fprintf(fid, ' & %s ', dcolbf(MZconst(zz,thisSam,dd,mm), '%5.2f', conststar));
                        fprintf(fid, ' & %s ', dcolbf(MZslopes(zz,thisSam,dd,mm), '%5.2f', slopestar));
                    else
                        fprintf(fid, ' &  \\multicolumn{1}{c}{---}');
                    end % thisNz

                end % models
            end % datalabel
            fprintf(fid, '\\\\ \n');
            % print se
            for dd = 1 : Ndata
                thisNz = matarray{dd,1}.Nz;
                for mm = 1 : Nmodels
                    if zz <= thisNz
                        % fprintf(fid, ' & (%5.2f) & (%5.2f) ', MZconstSE(zz,thisSam,dd,mm), MZslopesSE(zz,thisSam,dd,mm));
                        fprintf(fid, ' & (%5.2f) ', MZslopesSE(zz,thisSam,dd,mm));
                    else
                        fprintf(fid, ' &  \\multicolumn{1}{c}{---}');
                    end % thisNz
                end % models
            end % datalabel
            fprintf(fid, '\\\\\n');
            if zz < maxNz
                fprintf(fid, '\\midrule\n');
            end
        end
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        fclose(fid);
        type(fullfile(tabdir, tabname))
    end % sam


    %% tabulate results in one big table -- SLOPES AND INTERCEPTS IN SEPARATE PANELS
    % Zlabel  = matarray{1,1}.Zlabel; % assuming it contains maxNz elements
    % if length(Zlabel) ~= maxNz
    %     error('Zlabel mismatch')
    % end
    % for thisSam = 1 : Nsam
    %     Ndata   = length(DATALABELS);
    %     Nmodels = numel(models);
    %     Nstats  = 1;
    %     if ~doEXY1Q4
    %         tabname = sprintf('table-ZhatZtp1-twopanel-%s.tex', sam(thisSam).label);
    %     else
    %         tabname = sprintf('table-ZhatZtp1-twopanel-EXY1Q4-%s.tex', sam(thisSam).label);
    %     end
    %     tabcaption = sprintf('Zhat predictability regressions (two-panel version, %s)', sam(thisSam).pretty);
    %     if ~isempty(wrap)
    %         tabdir = wrap.dir;
    %         latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
    %     else
    %         tabdir = fullfile(localtemp, 'foo');
    %     end
    %     fid = fopen(fullfile(tabdir, tabname), 'wt');
    %     fprintf(fid, '\\begin{center}\n');
    %     fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.2', 1, Nstats * Nmodels * Ndata));
    %     fprintf(fid, '\\toprule\n');
    %     % header: datalabel row
    %     for dd = 1 : Ndata
    %         datalabel = DATALABELS{dd};
    %         fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats * Nmodels, datalabel);
    %     end
    %     fprintf(fid, '\\\\');
    %     offset = 1;
    %     for dd = 1 : Ndata
    %         fprintf(fid, '\\cmidrule(lr){%d-%d}', offset+1, offset+Nstats*Nmodels);
    %         offset = offset + Nstats*Nmodels;
    %     end
    %     fprintf(fid, '\n');
    %     % header: models
    %     fprintf(fid, 'Forecast ');
    %     for dd = 1 : Ndata
    %         for mm = 1 : Nmodels
    %             modeltype = models(mm).pretty;
    %             fprintf(fid, ' & \\multicolumn{%d}{c}{%s} ', Nstats, modeltype);
    %         end
    %     end
    %     fprintf(fid, '\n');
    %     fprintf(fid, '\\\\ \n');
    %     fprintf(fid, '\\midrule\n');
    %     fprintf(fid, '\\multicolumn{%d}{c}{PANEL A: slope $\\beta$} ', Nstats * Nmodels * Ndata);
    %     fprintf(fid, '\\\\ \n');
    %     fprintf(fid, '\\midrule\n');
    %     for zz = 2 : maxNz % skip lagged realized value
    %         fprintf(fid, '%s ', Zlabel{zz});
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     slopestar = abs(MZslopes(zz,thisSam,dd,mm) - 1) > 1.96 * MZslopesSE(zz,thisSam,dd,mm);
    %                     % conststar = abs(MZconst(zz,thisSam,dd,mm)) > 1.96 * MZconstSE(zz,thisSam,dd,mm);
    %                     % fprintf(fid, ' & %s ', dcolbf(MZconst(zz,thisSam,dd,mm), '%5.2f', conststar));
    %                     fprintf(fid, ' & %s ', dcolbf(MZslopes(zz,thisSam,dd,mm), '%5.2f', slopestar));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\ \n');
    %         % print se
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     % fprintf(fid, ' & (%5.2f) & (%5.2f) ', MZconstSE(zz,thisSam,dd,mm), MZslopesSE(zz,thisSam,dd,mm));
    %                     fprintf(fid, ' & (%5.2f) ', MZslopesSE(zz,thisSam,dd,mm));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\\n');
    %         % if zz < maxNz
    %         fprintf(fid, '\\midrule\n');
    %         % end
    %     end
    %     % fprintf(fid, '\\\\ \n');
    %     % fprintf(fid, '\\midrule\n');
    %     fprintf(fid, '\\multicolumn{%d}{c}{PANEL B: intercept $\\alpha$} ', Nstats * Nmodels * Ndata);
    %     fprintf(fid, '\\\\ \n');
    %     fprintf(fid, '\\midrule\n');
    %     for zz = 2 : maxNz % skip lagged realized value
    %         fprintf(fid, '%s ', Zlabel{zz});
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     % slopestar = abs(MZslopes(zz,thisSam,dd,mm) - 1) > 1.96 * MZslopesSE(zz,thisSam,dd,mm);
    %                     conststar = abs(MZconst(zz,thisSam,dd,mm)) > 1.96 * MZconstSE(zz,thisSam,dd,mm);
    %                     fprintf(fid, ' & %s ', dcolbf(MZconst(zz,thisSam,dd,mm), '%5.2f', conststar));
    %                     % fprintf(fid, ' & %s ', dcolbf(MZslopes(zz,thisSam,dd,mm), '%5.2f', slopestar));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\ \n');
    %         % print se
    %         for dd = 1 : Ndata
    %             thisNz = matarray{dd,1}.Nz;
    %             if zz <= thisNz
    %                 for mm = 1 : Nmodels
    %                     % fprintf(fid, ' & (%5.2f) & (%5.2f) ', MZconstSE(zz,thisSam,dd,mm), MZslopesSE(zz,thisSam,dd,mm));
    %                     fprintf(fid, ' & (%5.2f) ', MZconstSE(zz,thisSam,dd,mm));
    %                 end % models
    %             end % thisNz
    %         end % datalabel
    %         fprintf(fid, '\\\\\n');
    %         if zz < maxNz
    %             fprintf(fid, '\\midrule\n');
    %         end
    %     end
    %     fprintf(fid, '\\bottomrule\n');
    %     fprintf(fid, '\\end{tabular}\n');
    %     fprintf(fid, '\\end{center}\n');
    %     fclose(fid);
    %     type(fullfile(tabdir, tabname))
    % end % sam


    %%  wrap up loop
    dockAllFigures
    finishwrap
end

%% finish
finishscript