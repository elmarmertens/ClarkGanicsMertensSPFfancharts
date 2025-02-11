%% Estimate MDS and VAR0 models with and without noise -- single run for specific thisT

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
addpath ../matlabtoolbox/empbsbox/

%#ok<*NANMEAN>
%#ok<*NANVAR>
%#ok<*PFBNS>
%#ok<*UNRCH>
%#ok<*DATNM>
%#ok<*DATST>

%% clear workspace
clear variables

close all
fclose all;
clc
rng('default')

%% some parameters

rndStream  = getDefaultStream;

quicky     = false;  % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI'};

% select parameters
thisDate       = datenum(2024,1,1);
doSamStartSPF  = false;

NGAP     = 'BOP';
LABEL1   = 'WithNoiseY1Q4';
LABEL2   = 'noNoiseExY1Q4';

PRETTY1  = 'With Noise';
PRETTY2  = 'Without Noise';

MODELTYPES = {'MDS', 'VAR0'};

thinSTEPS = 1;

doStore = true;

%% general parameters
samStart = []; % 52 is 1981Q3

%% batch parameters
% SED-PARAMETERS-HERE
thisDate       = datenum(2019,10,1);

%% some settings


fontsize = 18;

color1   = colors4plots("darkblue");
color2   = colors4plots("darkred");

if quicky
    MCMCdraws    = 1e2;
    Nfedraws     = 10;
else
    MCMCdraws    = 3e3;
    Nfedraws     = 100;
end

%% loop over modeltypes
for mtp = MODELTYPES

    modeltype = mtp{:};
    switch modeltype
        case {'VAR0'}
            doVAR = true;
        case {'MDS'}
            doVAR = false;
        otherwise
            error('modeltype <<%s>> not recognized', modeltype);
    end

    %% loop over datalabels
    for d = 1 : length(DATALABELS)

        close all
        %% load dataset for variable
        datalabel = DATALABELS{d};

        datadir = fullfile('..', 'matdataKensington');
        matfilename = fullfile(datadir, sprintf('kensington%sdata', upper(datalabel)));
        load(matfilename, 'Ny', 'Nz', 'dates', 'datesQ', 'doNIPA', ...
            'Yfuture', 'Ylabel', 'Nbar', ...
            'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czhat', 'Czbar')

        ndxFixedHorizon     = 1 + 1 + (0:4);


        % select start of sample and cut sample
        if isempty(samStart)
            if doSamStartSPF
                samStart   = find(sum(~Znanny,2) > 1,1); % first availability of SPF
            else
                samStart   = find(any(~Znanny,2),1); % for CPI: leaves early sample with only Yrealized ...
            end
        end
        startDateLabel = datestr(dates(samStart), 'yyyyqq');
        dates          = dates(samStart:end); datesQ = datesQ(samStart:end);
        T              = length(dates);
        Zdata      = Zdata(samStart:end,:);
        Znanny     = Znanny(samStart:end,:);
        Cz         = Cz(:,:,samStart:end);
        ndxZbar    = cat(1, false(Nz-Nbar,1), true(Nbar,1));

        thisT    = find(dates <= thisDate, 1, 'last');
        thisDate = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');

        % process NGAP
        Ngap = switchNGAP(NGAP, Ny, datalabel, thisDate);


        modellabel = sprintf('%s-NoiseVsNoNoise-%strendHScycleSV2blockNoiseHS-Ngap%s-samStart%s-%s-thin%d', ...
            datalabel, modeltype, NGAP, startDateLabel, thisDateLabel, thinSTEPS);


        wrap = [];
        titlename = modellabel;
        if ~quicky
            initwrap
        end

        fprintf('Processing %s ... \n', modellabel)

        %% call MCMC sampler with Noise
        modelname1 = LABEL1; % sprintf('Ngap%d', Ngap1);
        fprintf('\n');
        OK = false;
        while ~OK
            try
                fprintf('STARTED mcmc1  ...\n')
                tic
                warning('off', 'backtrace');
                if doVAR
                    [YdensityDraws1, Ydraws1, ~, ETAdraws1, RESIDSVdraws1, ~, ...
                        Gdraws1, Glambda1, ~, YgapCONSTdraws1, ...
                        SVdraws1, SVtdraws1, hvcvdraws1, hrhodraws1, tdofdraws1, tScaleDraws1, sqrtsigmaDraws1, ~, mustarvardraws1, ...
                        NOISEdraws1, NOISEvarstate1, NOISEsig1, ...
                        ] = ...
                        mcmcsamplerVAR0trendHScycleSVt2blockNoiseHS(Zdata, Znanny, Cz, ...,
                        Ngap, ...
                        [], ...
                        ndxZbar, dates, ...
                        thisT, MCMCdraws, Nfedraws, rndStream, true, thinSTEPS);
                else
                    [YdensityDraws1, Ydraws1, ETAdraws1, RESIDSVdraws1, ...
                        SVdraws1, SVtdraws1, hvcvdraws1, hrhodraws1, tdofdraws1, tScaleDraws1, sqrtsigmaDraws1, ~, mustarvardraws1, ...
                        NOISEdraws1, NOISEvarstate1, NOISEsig1, ...
                        ] = ...
                        mcmcsamplerMDStrendHScycleSVt2blockNoiseHS(Zdata, Znanny, Cz, ...,
                        Ngap, ...
                        ndxZbar, dates, ...
                        thisT, MCMCdraws, Nfedraws, rndStream, true, thinSTEPS);
                end
                toc
                warning('on', 'backtrace');
                fprintf('DONE mcmc1.\n')
                OK = true;
            catch mcmcME
                warning('mcmc1 aborted at thisT=%d, message: %s', thisT, getReport(mcmcME, 'basic'))
            end
        end


        Nhorizons1 = size(YdensityDraws1,1);
        Nstates1   = Ngap + 1;

        %% collect PSRF and INFF for MCMC1
        % Ydraws1
        thesedraws          = Ydraws1;
        % NaN out the known values
        for ii = 1 : 6
            ndx = Znanny(1:thisT,ii);
            thesedraws(ii,~ndx,:) = NaN;
        end
        thesedraws = transpose(reshape(thesedraws, Nstates1 * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: Ydraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: Ydraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        % SV
        Nsv = 2;
        thesedraws = transpose(reshape(SVdraws1, Nsv * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: SVdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: SVdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        % SVt
        thesedraws = transpose(reshape(SVtdraws1, Nsv * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: SVtdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: SVtdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        % tdofdraws
        thesedraws = transpose(tdofdraws1);
        plotPSRF(thesedraws, sprintf('%s: tdofdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: tdofdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        if doVAR
            % Gdraws1
            thesedraws = transpose(reshape(Gdraws1, Ngap * Ngap, MCMCdraws));
            plotPSRF(thesedraws, sprintf('%s: Gdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
            plotINFF(thesedraws, sprintf('%s: Gdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
            % YgapCONSTdraws1
            thesedraws = transpose(YgapCONSTdraws1);
            plotPSRF(thesedraws, sprintf('%s: YgapCONSTdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
            plotINFF(thesedraws, sprintf('%s: YgapCONSTdraws1', LABEL1), strcat(modelname1, '-', modellabel), wrap);
        end

        %% call MCMC sampler without noise
        modelname2 = LABEL2; % sprintf('Ngap%d', Ngap1);
        % drop Y1Q4 from Zdata
        ndxQ4                = quarter(dates) == 4;
        ndxY1                 = 7; % y(-1), h=0:4 means six elements in Z before Y1
        Zdata2                = Zdata;
        Znanny2               = Znanny;
        Zdata2(ndxQ4, ndxY1)  = NaN;
        Znanny2(ndxQ4, ndxY1) = true;

        fprintf('\n');

        OK = false;
        while ~OK
            try
                fprintf('STARTED mcmc2 ...\n')
                tic
                warning('off', 'backtrace');
                if doVAR
                    [YdensityDraws2, Ydraws2, YhatRB2, ETAdraws2, RESIDSVdraws2, ZRESIDdraws2, ...
                        Gdraws2, Glambda2, ~, YgapCONSTdraws2, ...
                        SVdraws2, SVtdraws2, hvcvdraws2, hrhodraws2, tdofdraws2, tScaleDraws2, sqrtsigmaDraws2, ~, mustarvardraws2, ...
                        ] = ...
                        mcmcsamplerVAR0trendHScycleSVt2block(Zdata2, Znanny2, Cz, ...,
                        Ngap, ...
                        [], ...
                        thisT, MCMCdraws, Nfedraws, rndStream, true);
                else
                    [YdensityDraws2, Ydraws2, ETAdraws2, ETASVdraws2, ...
                        SVdraws2, SVtdraws2, hvcvdraws2, hrhodraws2, tdofdraws2, tScaleDraws2, sqrtsigmaDraws2, ~, mustarvardraws2, ...
                        ] = ...
                        mcmcsamplerMDStrendHScycleSVt2block(Zdata2, Znanny2, Cz, ...,
                        Ngap, ...
                        thisT, MCMCdraws, Nfedraws, rndStream, true);
                end
                toc
                warning('on', 'backtrace');
                fprintf('DONE mcmc2.\n')
                OK = true;
            catch mcmcME
                warning('mcmc2 aborted at thisT=%d, message: %s', thisT, getReport(mcmcME, 'basic'))
            end
        end
        Nhorizons2 = size(YdensityDraws2,1);
        Nstates2   = Ngap + 1;


        %% collect PSRF and INFF for MCMC2
        % Ydraws2
        thesedraws          = Ydraws2;
        % NaN out the known values
        for ii = 1 : 6
            ndx = Znanny(1:thisT,ii);
            thesedraws(ii,~ndx,:) = NaN;
        end
        thesedraws = transpose(reshape(thesedraws, Nstates2 * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: Ydraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: Ydraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        % SV
        Nsv = 2;
        thesedraws = transpose(reshape(SVdraws2, Nsv * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: SVdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: SVdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        % SVt
        thesedraws = transpose(reshape(SVtdraws2, Nsv * thisT, MCMCdraws));
        plotPSRF(thesedraws, sprintf('%s: SVtdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: SVtdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        % tdofdraws
        thesedraws = transpose(tdofdraws2);
        plotPSRF(thesedraws, sprintf('%s: tdofdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        plotINFF(thesedraws, sprintf('%s: tdofdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        if doVAR
            % Gdraws2
            thesedraws = transpose(reshape(Gdraws2, Ngap * Ngap, MCMCdraws));
            plotPSRF(thesedraws, sprintf('%s: Gdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
            plotINFF(thesedraws, sprintf('%s: Gdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
            % YgapCONSTdraws2
            thesedraws = transpose(YgapCONSTdraws2);
            plotPSRF(thesedraws, sprintf('%s: YgapCONSTdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
            plotINFF(thesedraws, sprintf('%s: YgapCONSTdraws2', LABEL2), strcat(modelname2, '-', modellabel), wrap);
        end

        %% adjust datevector etc for plotting
        datesT     = dates(1:thisT);

        %% plot Ydraws over time
        Nstates = max(Nstates1,Nstates2);
        for n = 1 : Nstates
            these  = squeeze(Ydraws1(min(n,Nstates1),:,:));
            mid1   = median(these, 2);
            tails1 = prctile(these, normcdf([-1 1]) * 100, 2);

            these  = squeeze(Ydraws2(min(n,Nstates2),:,:));
            mid2   = median(these, 2);
            tails2 = prctile(these, normcdf([-1 1]) * 100, 2);

            thisfig = figure;
            hold on
            h0 = plot(datesT, mid1, '-', 'color', color1, 'LineWidth', 2);
            plot(datesT, tails1, '--', 'color', color1, 'LineWidth', 1);
            h2 = plot(datesT, mid2, '-', 'color', color2, 'LineWidth', 2);
            plot(datesT, tails2, '--', 'color', color2, 'LineWidth', 1);
            if any(ismember(ndxFixedHorizon, n))
                thisZ = Zdata(1:thisT, n);
                thisZ(Znanny(1:thisT,n)) = NaN;
                hSPF = plot(datesT, thisZ, 'rx', 'LineWidth', .5);
            else
                hSPF = [];
            end

            xtickdates(datesT)
            if n == 1
                title(sprintf('data (lag 1)'))
            elseif ~isempty(hSPF)
                title(sprintf('forecast h=%d', n-2))
                legend([h0 h2 hSPF], modelname1, modelname2, 'SPF', 'location', 'best')
            else
                title(sprintf('forecast h=%d', n-2))
                legend([h0 h2], modelname1, modelname2, 'location', 'best')
            end
            wrapthisfigure(thisfig, sprintf('Ydraws%d-%s', n, modellabel), wrap)
        end
        clear Nstates


        %% plot SV
        Nsv = 2;
        % SVtdraws1 = SVtdraws1 ./ SVtdraws1(:,1,:);
        % SVtdraws2 = SVtdraws2 ./ SVtdraws2(:,1,:);
        % SVdraws1  = SVdraws1 ./ SVdraws1(:,1,:);
        % SVdraws2  = SVdraws2 ./ SVdraws2(:,1,:);

        midt1   = median(SVtdraws1,3);
        % tailst1 = prctile(SVtdraws1, normcdf([-1 1]) * 100, 3);
        midt2   = median(SVtdraws2,3);
        % tailst2 = prctile(SVtdraws2, normcdf([-1 1]) * 100, 3);

        mid1   = median(SVdraws1,3);
        % tails1 = prctile(SVdraws1, normcdf([-1 1]) * 100, 3);
        mid2   = median(SVdraws2,3);
        % tails2 = prctile(SVdraws2, normcdf([-1 1]) * 100, 3);

        for ii = 1 : Nsv
            thisfig = figure;
            subplot(2,1,1)
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT,midt1(ii,:), '-', 'color', color1, 'LineWidth', 2);
            % plot(datesT,squeeze(tailst1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
            h2 = plot(datesT,midt2(ii,:), '-', 'color', color2, 'LineWidth', 2);
            % plot(datesT,squeeze(tailst2(ii,:,:)), '-', 'color', color2, 'LineWidth', 1)
            xtickdates(datesT)
            legend([h1, h2], modelname1, modelname2, 'location', 'best')
            title('SV-t')
            subplot(2,1,2)
            hold on
            set(gca, 'fontsize', 18)
            h1 = plot(datesT,mid1(ii,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
            % plot(datesT,squeeze(tails1(ii,:,:)), '-', 'color', color1, 'LineWidth', 1)
            h2 = plot(datesT,mid2(ii,:), '-', 'color', color2, 'LineWidth', 2); %#ok<NASGU>
            % plot(datesT,squeeze(tails2(ii,:,:)), '-', 'color', color2, 'LineWidth', 1)
            xtickdates(datesT)
            % legend([h1, h2], modelname1, modelname2, 'location', 'best')
            title('SV')
            wrapthisfigure(thisfig, sprintf('SVblock%d-%s', ii, modellabel), wrap)
        end

        %% compare histograms of tdofdraws
        for ii = 1 : Nsv
            thisfig = figure;
            hold on
            set(gca, 'fontsize', 18)
            histogram(tdofdraws1(ii,:), 'Normalization', 'pdf', 'FaceAlpha', .5);
            histogram(tdofdraws2(ii,:), 'Normalization', 'pdf', 'FaceAlpha', .5);
            legend(modelname1, modelname2, 'location', 'best')
            title(sprintf('tdof Block %d', ii))
            wrapthisfigure(thisfig, sprintf('tdof-block%d-%s', ii, modellabel), wrap)
        end


        if doVAR
            Nstates   = Ngap + 1;
            H         = Nstates - 2 - 1;
            Hpools    = {0:3};
            ndxSV11   = 1:5;
            ndxSV22   = ndxSV11(end)+1:Ngap; % OK, let's hard code this. Could use cell arrays for easier extension to N blocks, but leave that for later

            %% Compute CG regressions model1
            CGbeta1   = NaN(H,MCMCdraws);
            CGpooled1 = NaN(length(Hpools),MCMCdraws);

            for m = 1 : MCMCdraws

                thisrho     = min(hrhodraws1(:,m), .99);
                ESV         = exp(0.5 .* diag(hvcvdraws1(:,:,m)) ./ (1 - thisrho.^2));
                ESV         = ESV .* tdofdraws1(:,m) ./ (tdofdraws1(:,m) - 2); % adjustment for t
                sigTrend    = mustarvardraws1(m);
                sqrtSigma   = sqrtsigmaDraws1(:,:,m);
                sqrtSigma(:,ndxSV11) = sqrtSigma(:,ndxSV11) .* ESV(1);
                sqrtSigma(:,ndxSV22) = sqrtSigma(:,ndxSV22) .* ESV(2);
                Sigma       = sqrtSigma * sqrtSigma';

                G           = Gdraws1(:,:,m);
                VCV         = dlyapdoubling(G, Sigma);

                CGpooled1(:,m) = CGregressionpooledTC(G,VCV,sigTrend,Hpools);
                % CG regression
                CGbeta1(:,m) = CGregressionTC(G,VCV,sigTrend,H);

            end
            %% Compute CG regressions model2
            CGbeta2   = NaN(H,MCMCdraws);
            CGpooled2 = NaN(length(Hpools),MCMCdraws);

            for m = 1 : MCMCdraws

                thisrho     = min(hrhodraws2(:,m), .99);
                ESV         = exp(0.5 .* diag(hvcvdraws2(:,:,m)) ./ (1 - thisrho.^2));
                ESV         = ESV .* tdofdraws2(:,m) ./ (tdofdraws2(:,m) - 2); % adjustment for t
                sigTrend    = mustarvardraws2(m);
                sqrtSigma   = sqrtsigmaDraws2(:,:,m);
                sqrtSigma(:,ndxSV11) = sqrtSigma(:,ndxSV11) .* ESV(1);
                sqrtSigma(:,ndxSV22) = sqrtSigma(:,ndxSV22) .* ESV(2);
                Sigma       = sqrtSigma * sqrtSigma';

                G           = Gdraws2(:,:,m);
                VCV         = dlyapdoubling(G, Sigma);

                CGpooled2(:,m) = CGregressionpooledTC(G,VCV,sigTrend,Hpools);
                % CG regression
                CGbeta2(:,m) = CGregressionTC(G,VCV,sigTrend,H);

            end
            clear Ngap Nstates

            %% plot densities of CG slopes
            CGbeta1 = transpose(CGbeta1);
            CGbeta2 = transpose(CGbeta2);
            thisfig = figure;
            for h = 1 : H
                subplot(ceil(H/2),2,h)
                [~, h1, h2] = plotpriorposteriordraws(CGbeta1(:,h),CGbeta2(:,h));
                h1.Color=color1;
                h2.Color=color2;
                xline(0, '-', 'LineWidth',1)
                title(sprintf('h = %d', h))
            end
            sgtitle(sprintf('%s \n CG slopes', datalabel))
            wrapthisfigure(thisfig, sprintf('CGslopes%d-%s', H, modellabel), wrap)

            %% plot pooled CG slopes
            for ii = 1 : length(Hpools)
                thisfig = figure;
                [~,h1,h2] = plotpriorposteriordraws(CGpooled1(ii,:),CGpooled2(ii,:));
                h1.Color=color1;
                h2.Color=color2;
                xline(0, '-', 'LineWidth',1)
                thisPool = Hpools{ii};
                title(sprintf('%s -- Pooled CG slope (h=%d:%d)', datalabel, thisPool(1), thisPool(end)))
                wrapthisfigure(thisfig, sprintf('CGpooled%d-%s', ii, modellabel), wrap)
            end
        end

        %% plot medians of noise levels
        Zbarlabel = Zlabel(ndxZbar);
        ndxQ4 = quarter(datesT) == 4;
        mid = median(NOISEdraws1,3);
        % mid2 = median(NOISEdraws2,3);
        % with Q4 diamond
        for nn = 1 : Nbar
            firstObs = find(~isnan(mid(nn,:)), 1);
            if ~isempty(firstObs)
                thisfig = figure;
                hold on
                set(gca, 'fontsize', 18)
                h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2);
                h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
                xtickdates(datesT(firstObs:end))
                plotOrigin;
                wrapthisfigure(thisfig, sprintf('NOISElevelsWITHQ4-y%d-%s', nn, modellabel), wrap, [], [], [], [], true);
                legend(h1q4, 'Q4', 'location', 'best')
                wrapthisfigure(thisfig, sprintf('NOISElevelsWITHQ4-y%d-%s-WITHLEGEND', nn, modellabel), wrap, [], [], [], [], true);
                % title(sprintf('NOISE levels %s', Zbarlabel{nn}))
                % wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s-WITHLEGENDTITLE', nn, modellabel), wrap)
            end
        end
        % without Q4 diamond
        for nn = 1 : Nbar
            firstObs = find(~isnan(mid(nn,:)), 1);
            if ~isempty(firstObs)
                thisfig = figure;
                hold on
                set(gca, 'fontsize', 18)
                h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
                % h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);
                xtickdates(datesT(firstObs:end))
                plotOrigin;
                wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s', nn, modellabel), wrap, [], [], [], [], true);
                % legend([h1, h1q4, h2, h2q4], modelname1, 'Q4', modelname2, 'Q4', 'location', 'best')
                % wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s-WITHLEGEND', nn, modellabel), wrap, [], [], [], [], true);
                % title(sprintf('NOISE levels %s', Zbarlabel{nn}))
                % wrapthisfigure(thisfig, sprintf('NOISElevels-y%d-%s-WITHLEGENDTITLE', nn, modellabel), wrap)
            end
        end

        %% plot medians of noise variances
        mid = median(NOISEvarstate1,3);
        for nn = 1 : Nbar
            firstObs = find(~isnan(mid(nn,:)), 1);
            if ~isempty(firstObs)
                thisfig = figure;
                hold on
                set(gca, 'fontsize', 18)
                h1 = plot(datesT, mid(nn,:), '-', 'color', color1, 'LineWidth', 2); %#ok<NASGU>
                h1q4 = plot(datesT(ndxQ4), mid(nn,ndxQ4), 'd', 'color', color1, 'LineWidth', 2);

                % h2 = plot(datesT, mid2(nn,:), ':', 'color', color2, 'LineWidth', 2);
                % h2q4 = plot(datesT(ndxQ4), mid2(nn,ndxQ4), 'd', 'color', color2, 'LineWidth', 2);
                xtickdates(datesT(firstObs:end))
                plotOrigin;
                % legend([h1, h1q4, h2, h2q4], modelname1, 'Q4', modelname2, 'Q4', 'location', 'best')
                title(sprintf('NOISE variances %s', Zbarlabel{nn}))
                wrapthisfigure(thisfig, sprintf('NOISEvariances-y%d-%s', nn, modellabel), wrap)
            end
        end

        %% plot predictive density of Y
        SPFcolor = colors4plots('yellow');
        Zdata(Znanny) = NaN;
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density
        mid1   = mean(YdensityDraws1,2);
        tails1 = prctile(YdensityDraws1, [5 normcdf([-1 1]) * 100 95], 2);
        h1 = plotCI(mid1, tails1, 0:Nhorizons1-1,[], '-', 'color', 'black', 'LineWidth', 4);

        h2 = plot(0:Nhorizons2-1,mean(YdensityDraws2,2), '-', 'color', colors4plots('orange'), 'LineWidth', 4);
        plot(0:Nhorizons2-1,prctile(YdensityDraws2, normcdf([-1 1]) * 100, 2), '-', 'color', colors4plots('orange'), 'LineWidth', 2)
        plot(0:Nhorizons2-1,prctile(YdensityDraws2, [5 95], 2), '-.', 'color', colors4plots('orange'), 'LineWidth', 2)

        % redo mid1 as dashes
        plot(0:Nhorizons1-1, mid1, 'k--', 'linewidth', 3);

        % finish plot
        xmax = max(Nhorizons1,Nhorizons2) - 1;
        xticks(0 : 2 : xmax)
        xlim([0 xmax])
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s', modellabel), wrap)
        hl = legend([h1 h2], PRETTY1, PRETTY2, 'location', 'best');
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHLEGEND', modellabel), wrap)
        delete(hl)


        % add SPF plots
        if doNIPA
            MAweights       = -2 : 4;
        else
            MAweights       = 1 : 4;
        end
        ndxB            = 2 + 4 + 1;
        horizonsB       = (4 - datesQ(thisT)) + MAweights;
        ndxC            = 2 + 4 + 2;
        horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
        ndxD            = 2 + 4 + 3;
        horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
        hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', SPFcolor, 'linewidth', 2);
        if ~Znanny(thisT,ndxB)
            hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', SPFcolor, 'linewidth', 2);
        end
        if ndxC <= Nz && ~Znanny(thisT,ndxC)
            hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', SPFcolor, 'linewidth', 2);
        end
        if ndxD <= Nz && ~Znanny(thisT,ndxD)
            hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', SPFcolor, 'linewidth', 2);
        end
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHSPF', modellabel), wrap)

        hl = legend([h1 h2 hSPF], PRETTY1, PRETTY2,'SPF', 'location', 'best');
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHSPFLEGEND', modellabel), wrap)

        title(sprintf('%s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('Ypredictivedensity-%s-WITHSPFLEGENDTITLE', modellabel), wrap)
        YtermLIM = ylim; % for Yterm plot

        %% plot term structures Ydraws (at end of sample)

     
        thisH = max(Nstates1,Nstates2)-2;

        Zdata(Znanny) = NaN;
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % add plot of term structure estimates
        h1 = plot(0:Nstates1-2,mean(Ydraws1(2:Nstates1,end,:),3), '-', 'color', color1, 'Linewidth', 3);
        theseTtails = squeeze(prctile(Ydraws1(2:Nstates1,end,:), normcdf([-1 1]) * 100, 3));
        plot(0:Nstates1-2,theseTtails, '-', 'color', color1, 'Linewidth', 1);

        h2 = plot(0:Nstates2-2,mean(Ydraws2(2:Nstates2,end,:),3), '--', 'color', color2, 'Linewidth', 3);
        theseTtails = squeeze(prctile(Ydraws2(2:Nstates2,end,:), normcdf([-1 1]) * 100, 3));
        plot(0:Nstates2-2,theseTtails, '--', 'color', color2, 'Linewidth', 1);
        % % add SPF plots
        % if doNIPA
        %     MAweights       = -2 : 4;
        % else
        %     MAweights       = 1 : 4;
        % end
        % ndxB            = 2 + 4 + 1;
        % horizonsB       = (4 - datesQ(thisT)) + MAweights;
        % ndxC            = 2 + 4 + 2;
        % horizonsC       = (4 - datesQ(thisT)) + 4 + MAweights;
        % ndxD            = 2 + 4 + 3;
        % horizonsD       = (4 - datesQ(thisT)) + 8 + MAweights;
        % SPFcolor =  [0 0 0];
        % hSPF  = plot(0:4, Zdata(thisT, ndxFixedHorizon), 'o', 'color', SPFcolor, 'linewidth', 2);
        % if ~Znanny(thisT,ndxB)
        %     hSPFb = plot(horizonsB, Zdata(thisT, repmat(ndxB, 1, length(horizonsB))), ':o', 'color', SPFcolor, 'linewidth', 2);
        % end
        % if ndxC <= Nz && ~Znanny(thisT,ndxC)
        %     hSPFc = plot(horizonsC, Zdata(thisT, repmat(ndxC, 1, length(horizonsC))), ':o', 'color', SPFcolor, 'linewidth', 2);
        % end
        % if ndxD <= Nz && ~Znanny(thisT,ndxD)
        %     hSPFd = plot(horizonsD, Zdata(thisT, repmat(ndxD, 1, length(horizonsD))), ':o', 'color', SPFcolor, 'linewidth', 2);
        % end

        xticks(0 : 2 : thisH)
        xlim([0 thisH])
        wrapthisfigure(thisfig, sprintf('Yterm-%s', modellabel), wrap, [], [], [], [], true);
        legend([h1 h2], PRETTY1, PRETTY2,'location', 'best');
        wrapthisfigure(thisfig, sprintf('Yterm-%s-WITHLEGEND', modellabel), wrap, [], [], [], [], true);
        title(sprintf('Term structure: %s per %s', datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('Yterm-%s-WITHLEGENDTITLE', modellabel), wrap)


        %% for model1: plot predictive density together with term structure of Ydraws (at end of sample)
        colorTerm = colors4plots('black');
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density
        hDense = plot(0:Nhorizons1-1,mean(YdensityDraws1,2), '-', 'color', color1, 'LineWidth', 4);
        plot(0:Nhorizons1-1,prctile(YdensityDraws1, normcdf([-1 1]) * 100, 2), '-', 'color', color1, 'LineWidth', 2)
        plot(0:Nhorizons1-1,prctile(YdensityDraws1, [5 95], 2), '-.', 'color', color1, 'LineWidth', 2)
        % plot term structure estimates
        hTerm = plot(0:Nstates1-2,mean(Ydraws1(2:Nstates1,end,:),3), '--', 'color', colorTerm, 'Linewidth', 4);
        theseTtails = squeeze(prctile(Ydraws1(2:Nstates1,end,:), normcdf([-1 1]) * 100, 3));
        plot(0:Nstates1-2,theseTtails, '--', 'color', colorTerm, 'Linewidth', 1);
        % axis settings
        xmax = max(Nhorizons1,Nhorizons2) - 1;
        xticks(0 : 2 : xmax)
        xlim([0 xmax])
        legend([hDense hTerm], 'Predictive Density', 'SPF-consistent term structure','location', 'best');
        title(sprintf('%s: %s per %s', modelname1, datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGENDTITLE', modelname1, modellabel), wrap)

        %% for model2: plot predictive density together with term structure of Ydraws (at end of sample)
        colorTerm = colors4plots('black');
        thisfig = figure;
        hold on
        set(gca, 'fontsize', fontsize)
        % plot density
        hDense = plot(0:Nhorizons2-1,mean(YdensityDraws2,2), '-', 'color', color2, 'LineWidth', 4);
        plot(0:Nhorizons2-1,prctile(YdensityDraws2, normcdf([-1 1]) * 100, 2), '-', 'color', color2, 'LineWidth', 2)
        plot(0:Nhorizons2-1,prctile(YdensityDraws2, [5 95], 2), '-.', 'color', color2, 'LineWidth', 2)
        % plot term structure estimates
        hTerm = plot(0:Nstates2-2,mean(Ydraws2(2:Nstates2,end,:),3), '--', 'color', colorTerm, 'Linewidth', 4);
        theseTtails = squeeze(prctile(Ydraws2(2:Nstates2,end,:), normcdf([-1 1]) * 100, 3));
        plot(0:Nstates2-2,theseTtails, '--', 'color', colorTerm, 'Linewidth', 1);
        % axis settings
        xmax = max(Nhorizons1,Nhorizons2) - 1;
        xticks(0 : 2 : xmax)
        xlim([0 xmax])
        legend([hDense hTerm], 'Predictive Density', 'SPF-consistent term structure','location', 'best');
        title(sprintf('%s: %s per %s', modelname2, datalabel, datestr(datesT(end), 'yyyyqq')))
        wrapthisfigure(thisfig, sprintf('YpredictivedensityTerm-%s-%s-WITHLEGENDTITLE', modelname2, modellabel), wrap)

        %% wrap up
        dockAllFigures
        finishwrap

        if ~isempty(wrap) && doStore
            close all
            filename = fullfile(wrap.dir, modellabel);
            save(filename, '-v7.3');
        end

    end % datalabel
end % modeltype

%% finish / clean up
finishscript
