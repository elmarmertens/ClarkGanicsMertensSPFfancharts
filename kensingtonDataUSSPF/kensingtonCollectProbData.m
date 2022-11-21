%% CGM US SPF probability forecasts
%% NOTE 1: in the SPF tables, the responses start from the RIGHTMOST bin, while the bins below start from the LEFTMOST ones
%% NOTE 2: kensingtonCollectData.m must be run before running this script (constructing the structure objects with bin data requires the output of that script)
clear; close all; clc;

%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/

%% bin endpoints, these are valid as of the 2022Q2 SPF round
% if SPF bins change after that date, this need to be extended

% PGDP
PRPGDP1_bins=[-3,-2.05:1:9.95]; % 1968Q4-1973Q1
PRPGDP2_bins=[-1,-0.05:1:11.95]; % 1973Q2-1974Q3
PRPGDP3_bins=[3,3.95:1:15.95]; % 1974Q4-1981Q2
PRPGDP4_bins=[4,5.95:2:11.95]; % 1981Q3-1985Q1  next year got included
PRPGDP5_bins=[2,3.95:2:9.95]; % 1985Q2-1991Q4
PRPGDP6_bins=[0,0.95:1:7.95]; % 1992Q1-2013Q4
PRPGDP7_bins=[0,0.45:0.5:3.95]; % 2014Q1-present
% RGDP
PRGDP1_bins=[-3,-2.05:1:9.95]; % 1968Q4-1973Q1
PRGDP2_bins=[-1,-0.05:1:11.95]; % 1973Q2-1974Q3
PRGDP3_bins=[3,3.95:1:15.95]; % 1974Q4-1981Q2
PRGDP4_bins=[-2,-0.05:2:5.95]; % 1981Q3-1991Q4  next year got included
PRGDP5_bins=[-2,-1.05:1:5.95]; % 1992Q1-2009Q1
PRGDP6_bins=[-3,-2.05:1:5.95]; % 2009Q2-2020Q1  years 3 and 4 got included
PRGDP7_bins=[-12,-6.05,-3.05,-0.5,1.45,2.45,3.95,6.95,9.95,15.95]; % 2020Q2-present
% UNEMP
PRUNEMP1_bins=[6,6.95,7.45,7.95,8.45,8.95,9.45,9.95,10.95]; % 2009Q2-2013Q4     next, years 3 and 4 included
PRUNEMP2_bins=[4,4.95,5.45,5.95,6.45,6.95,7.45,7.95,8.95]; % 2014Q1-2020Q1
PRUNEMP3_bins=[3,3.95,4.95,5.95,6.95,7.95,9.95,11.95,14.95]; % 2020Q2-present

%% select choice of data
DATALABELS = {'PRGDP','PRPGDP','PRUNEMP'};
% prob.xlsx needs to be downloaded from the Philadelphia Fed's Real-Time Data Research Center website: 
% https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/probability-variables
SPFimport=importdata('prob.xlsx');

% ADAPT DATES HERE
dates    = genrQdates(1968,2022); % only adjust the second argument (if needed)
dates    = dates(4:end-2); % only adjust the second argument (if needed)
% the vector "dates" must be the dates as stated in the SPF file (1968:Q4 through 2022:Q2)

for d=1:length(DATALABELS)

    data_act = SPFimport.data.(DATALABELS{d});

    % check dates
    Y = data_act(:,1);
    if ~isequal(Y, year(dates))
        error('inconsistent date vector in the SPF file (years)');
    end
    Q = data_act(:,2);
    if ~isequal(Q, quarter(dates))
        error('inconsistent date vector in the SPF file (quarters)');
    end

    data_probs=data_act(:,3:end);
    allgrids=NaN(size(dates,1),60); % first column: number of horizons,
                                    % after that repeating the bin edges (as many times as the number of horizons - adjust this if new horizons are added after 2022Q2)

    for t=1:size(dates,1)

        if strcmp(DATALABELS{d},'PRGDP')==1

            if     dates(t)>=datenum('1968Q4','yyyyQQ') && dates(t)<=datenum('1973Q1','yyyyQQ')
                numhz=1;
                thisbin=PRGDP1_bins;
            elseif dates(t)>=datenum('1973Q2','yyyyQQ') && dates(t)<=datenum('1974Q3','yyyyQQ')
                numhz=1;
                thisbin=PRGDP2_bins;
            elseif dates(t)>=datenum('1974Q4','yyyyQQ') && dates(t)<=datenum('1981Q2','yyyyQQ')
                numhz=1;
                thisbin=PRGDP3_bins;
            elseif dates(t)>=datenum('1981Q3','yyyyQQ') && dates(t)<=datenum('1991Q4','yyyyQQ')
                numhz=2;
                thisbin=PRGDP4_bins;
            elseif dates(t)>=datenum('1992Q1','yyyyQQ') && dates(t)<=datenum('2009Q1','yyyyQQ')
                numhz=2;
                thisbin=PRGDP5_bins;
            elseif dates(t)>=datenum('2009Q2','yyyyQQ') && dates(t)<=datenum('2020Q1','yyyyQQ')
                numhz=4;
                thisbin=PRGDP6_bins;
            elseif dates(t)>=datenum('2020Q2','yyyyQQ') % if SPF bins for PRGDP change after 2022Q2, this needs to be adjusted
                numhz=4;
                thisbin=PRGDP7_bins;
            end

        elseif strcmp(DATALABELS{d},'PRPGDP')==1

            if     dates(t)>=datenum('1968Q4','yyyyQQ') && dates(t)<=datenum('1973Q1','yyyyQQ')
                numhz=1;
                thisbin=PRPGDP1_bins;
            elseif dates(t)>=datenum('1973Q2','yyyyQQ') && dates(t)<=datenum('1974Q3','yyyyQQ')
                numhz=1;
                thisbin=PRPGDP2_bins;
            elseif dates(t)>=datenum('1974Q4','yyyyQQ') && dates(t)<=datenum('1981Q2','yyyyQQ')
                numhz=1;
                thisbin=PRPGDP3_bins;
            elseif dates(t)>=datenum('1981Q3','yyyyQQ') && dates(t)<=datenum('1985Q1','yyyyQQ')
                numhz=2;
                thisbin=PRPGDP4_bins;
            elseif dates(t)>=datenum('1985Q2','yyyyQQ') && dates(t)<=datenum('1991Q4','yyyyQQ')
                numhz=2;
                thisbin=PRPGDP5_bins;
            elseif dates(t)>=datenum('1992Q1','yyyyQQ') && dates(t)<=datenum('2013Q4','yyyyQQ')
                numhz=2;
                thisbin=PRPGDP6_bins;
            elseif dates(t)>=datenum('2014Q1','yyyyQQ') % if SPF bins for PRPGDP change after 2022Q2, this needs to be adjusted
                numhz=2;
                thisbin=PRPGDP7_bins;
            end

        elseif strcmp(DATALABELS{d},'PRUNEMP')==1

            if     dates(t)>=datenum('2009Q2','yyyyQQ') && dates(t)<=datenum('2013Q4','yyyyQQ')
                numhz=4;
                thisbin=PRUNEMP1_bins;
            elseif dates(t)>=datenum('2014Q1','yyyyQQ') && dates(t)<=datenum('2020Q1','yyyyQQ')
                numhz=4;
                thisbin=PRUNEMP2_bins;
            elseif dates(t)>=datenum('2020Q2','yyyyQQ') % if SPF bins for PRUNEMP change after 2022Q2, this needs to be adjusted
                numhz=4;
                thisbin=PRUNEMP3_bins;
            end

        end

        if ~isnan(data_act(t,3)) % only save data if it exists
            allgrids(t,1)=numhz;
            allgrids(t,2:numhz*length(thisbin)+1)=repmat(thisbin,1,numhz);
        end

    end

    tmp = sum(any(~isnan(allgrids),1)); % keep only up to the last non-NaN column 
    allgrids=allgrids(:,1:tmp);

    if strcmp(DATALABELS{d},'PRPGDP')==1
        PGDP.bins=allgrids;
        PGDP.probs=data_probs;
        PGDP.dates=dates;
    elseif strcmp(DATALABELS{d},'PRGDP')==1
        RGDP.bins=allgrids;
        RGDP.probs=data_probs;
        RGDP.dates=dates;
    elseif strcmp(DATALABELS{d},'PRUNEMP')==1
        UNRATE.bins=allgrids;
        UNRATE.probs=data_probs;
        UNRATE.dates=dates;
    end


end

save('kensingtonPROB.mat','PGDP','RGDP','UNRATE', '-v7.3');

%% save data --- separate files for each variable

% DATALABELS = {'RGDP','PGDP','UNRATE'};
%
% for d = 1 : length(DATALABELS)
%     datalabel = DATALABELS{d};
%     save(sprintf('kensington%sdataPROB.mat', datalabel), '-struct', datalabel, '-v7.3');
% end

%% construct structure objects with bin data
DATALABELS = {'RGDP','PGDP','UNRATE'};

probmat = load('kensingtonPROB.mat');

SPF_diffs_allvars = repmat(struct,length(DATALABELS),1);

for d = 1 : length(DATALABELS)
    datalabel = DATALABELS{d};
    this = probmat.(datalabel);

    dates = this.dates;
    T     = length(dates);

    tndx = find(any(~isnan(this.bins(:,2:end)), 2))';

    SPFhistograms = repmat(struct, T, 1);

    SPF_diffs     = NaN(T,4);

    % RPS starting with next calendar year's forecasts, loading outcomes
    FUTUREdata           = load(sprintf('kensington%sdata', upper(datalabel)),'YBARfuture', 'YCURRENTfuture');
    Nbar                 = size(FUTUREdata.YBARfuture, 2);
    fcstSPFdrps          = NaN(T,Nbar);
    fcstSPFdrpsCURRENT   = NaN(T,1);

    for t = tndx
        Nhz         = this.bins(t,1);
        that        = this.bins(t,2:end);
        Nbins       = find(~isnan(that),1,'last') / Nhz + 1;

        gridEdges   = that(1:Nbins-1);

        % collect parameters
        SPFhistograms(t).Nbins  = Nbins;
        SPFhistograms(t).Nhz    = Nhz;

        % collect bin data
        SPFhistograms(t).binEdges  = [gridEdges, Inf]';

        % compute center point (assuming closed bins)
        binEdgesLeftClosed         = [gridEdges(1)-2*(gridEdges(2)-gridEdges(1)), gridEdges]';
        binEdgesRightClosed        = [gridEdges, gridEdges(end)+2*(gridEdges(end)-gridEdges(end-1))]';
        SPFhistograms(t).binCenter = 0.5 * (binEdgesLeftClosed + binEdgesRightClosed);

        % collect probabilities
        histProbs = NaN(Nbins, Nhz);
        for hz = 1 : Nhz
            that            = this.probs(t,(hz-1)*Nbins+(1:Nbins)) / 100;
            histProbs(:,hz) = fliplr(that); % since SPF stores highest-to-lowest bins from left to right
        end
        SPFhistograms(t).probs = histProbs;
        SPFhistograms(t).cdf   = cumsum(histProbs, 1);

        SPF_diffs(t,1:Nhz) = 100*(sum(histProbs,1,'omitnan') - 1);
        %checkdiff(sum(histProbs) - 1); % TODO: renormalize to add to one exactly?

        % construct upper/lower contours of CDF
        [ upperEdges, lowerEdges, upperCDF, lowerCDF ] = deal(NaN(Nbins + Nbins - 1, Nhz));
        for hz = 1 : Nhz

            % upper contour
            theseEdges              = SPFhistograms(t).binEdges;
            thatCDF                 = SPFhistograms(t).cdf(:,hz);
            theseEdges              = cat(1, theseEdges, theseEdges(1:end-1));
            thatCDF                 = cat(1, thatCDF, thatCDF(2:end));
            [upperEdges(:,hz), ndx] = sort(theseEdges);
            upperCDF(:,hz)          = thatCDF(ndx);

            % lower contour
            theseEdges              = SPFhistograms(t).binEdges;
            thatCDF                 = SPFhistograms(t).cdf(:,hz);
            theseEdges              = cat(1, theseEdges(2:end), theseEdges); % note: order in which things are stacked matters here
            thatCDF                 = cat(1, thatCDF(1:end-1), thatCDF);
            [lowerEdges(:,hz), ndx] = sort(theseEdges);
            lowerCDF(:,hz)          = thatCDF(ndx);

            % calculate RPS based on Clements (2019), p.72, eq. 5.9
            % it starts from next year's forecasts
            if hz == 1
                thisBinRPS = SPFhistograms(t).binEdges;
                thisFuture = FUTUREdata.YCURRENTfuture(t); 
                if ~isnan(thisFuture)
                    thisHit               = double(thisFuture < thisBinRPS);
                    fcstSPFdrpsCURRENT(t) = sum((SPFhistograms(t).cdf(:,hz) - thisHit).^2);
                end
            else
                thisBinRPS     = SPFhistograms(t).binEdges;
                thisFuture = FUTUREdata.YBARfuture(t,hz-1); % YBAR starts with next calendar year's outcome
                if ~isnan(thisFuture)
                    thisHit             = double(thisFuture < thisBinRPS);
                    fcstSPFdrps(t,hz-1) = sum((SPFhistograms(t).cdf(:,hz) - thisHit).^2);
                end
            end

        end

        SPFhistograms(t).upperEdges = upperEdges;
        SPFhistograms(t).lowerEdges = lowerEdges;
        SPFhistograms(t).upperCDF   = upperCDF;
        SPFhistograms(t).lowerCDF   = lowerCDF;

    end

    if Nbar ~= max(cell2mat(arrayfun(@(x) x.Nhz, SPFhistograms, 'UniformOutput', false))) - 1
        error('Nbar mismatch')
    end
    save(sprintf('kensington%sdataHISTOGRAMS.mat', datalabel), 'SPFhistograms', 'fcstSPFdrps', 'fcstSPFdrpsCURRENT', 'dates', '-v7.3');

    SPF_diffs_allvars(d).SPF_diffs = SPF_diffs;
end


%% plot differences from (in %)

wrap = [];
titlename = sprintf('SPFdiffs');
initwrap

legtext = {'Current year','Next year','2 years ahead','3 years ahead'};
close all;

marker_vec={'none','s','d','o'};

for d = 1 : length(DATALABELS)

    datalabel = DATALABELS{d};
    this_diffs = SPF_diffs_allvars(d).SPF_diffs;
    SPF_histogram_tmp = load(sprintf('kensington%sdataHISTOGRAMS.mat', datalabel), 'SPFhistograms');
    thisdates = probmat.(datalabel).dates;
    T = length(thisdates);
    Nhz_max = 0;
    for t = 1: T
       if SPF_histogram_tmp.SPFhistograms(t).Nhz > Nhz_max
            Nhz_max = SPF_histogram_tmp.SPFhistograms(t).Nhz;
       end
    end

    thisfig = figure;
    hold on
    for hz = 1:Nhz_max
        plot(thisdates,this_diffs(:,hz),'Marker',marker_vec{hz});
    end

    ylabel('%');

    xline(thisdates(thisdates==datenum('1992Q1','yyyyQQ')),'r','LineWidth',1.5);
    
    xlim([thisdates(find(~isnan(this_diffs(:,1)),1,"first")) thisdates(T)]);
    datetick('x','yyyy','keeplimits');

    title(datalabel);
    legend(legtext{1:Nhz_max},'Box','on','Location','northeast');

    wrapthisfigure(thisfig, sprintf('SPFdiffs-%s',datalabel), wrap, [], [], [], [], false);

    hold off

    % clean up
        if ~isempty(wrap)
            close(thisfig);
        end

end

% finish datalabel loop
finishwrap
