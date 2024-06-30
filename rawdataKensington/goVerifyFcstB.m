%% construct data set for project "kensington" and check fcst B (RGDP and PGDP only)


% USE getKensingtonData.ipynb to collect latest raw data
% from SPF, RTDSM and FRED

% run collectRTData first

% PLEASE ADAPT start and end dates of sample

%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtexbox
addpath ../matlabtoolbox/emtools/


%#ok<*NANVAR>
%#ok<*UNRCH>
%#ok<*DATST>
%#ok<*DATNM>

%% clear workspace
clear
close all
fclose all;
clc

%% ADAPT END OF SAMPLE HERE
SPFsamstop  = datenum(2024,4,1);  % 2024Q2
% dates follow the FRED convention of placing timestamp at beginning of period

%% initwrap
initwrap
fontsize = 18;

%% select choice of data
DATALABELS = {'UNRATE', 'CPI', 'PGDP', 'RGDP'};

for d = 1 : length(DATALABELS)
    
    datalabel = DATALABELS{d};
    
    doNIPA = any(ismember({'PGDP', 'RGDP'}, datalabel));
    
    %% load RealTime Data
    RTdata = load(fullfile('..', 'matdataKensington', sprintf('kensington%sdataRT', datalabel)));
    
    %% prepare dates
    % quarters are time-stamped on the first day of the quarter (FRED convention)
    dates     = genrQdates(1968,year(SPFsamstop));
    dates     = dates(dates <= SPFsamstop);
    
    
    % realizedDataDates points to date vector starting at first lag of realized Y till last observation (T+1 observations)
    if doNIPA
        realizedDataDates = dates(3:end);
    else
        switch upper(datalabel)
            case {'UNRATE'}
                realizedDataDates = dates(3:end);     % include 68Q3 to grab lagged values
            case {'CPI'}
                realizedDataDates = dates(2:end);     % include 68Q2 to grab lagged *CPI* level for construction of lagged inflation data
            otherwise
                error('datalabel %s not yet supported', upper(datalabel))
        end
    end
    dates      = dates(4:end); % SPF dates begin only in 1968Q4
    T          = length(dates);
    
    
    %% load SPF responses for the level of PGDP and convert into inflation rates
    % the results of this cell are to be stored in Ydata(:,2:6)
    
    switch upper(datalabel)
        case 'UNRATE'
            SPFimport = importdata('mean_UNEMP_level.xlsx'); % this is the original xls file converted into csv
        case {'CPI'}
            SPFimport = importdata('mean_CPI_level.xlsx'); % this is the original xls file converted into csv
        case 'PGDP'
            SPFimport = importdata('mean_PGDP_level.xlsx'); % this is the original xls file converted into csv
        case 'RGDP'
            SPFimport = importdata('mean_RGDP_level.xlsx'); % this is the original xls file converted into csv
        otherwise
            error('datalabel %s not yet supported', upper(datalabel))
    end
    
    % col 1: year
    % col 2: quarter
    % col 3: "forecasts" for previous quarter
    % col 4: nowcast
    % col 5: forecast for next quarter
    % col 6: 2-quarter forecast
    % col 7: 3-quarter forecast
    % col 8: 4-quarter forecast
    % col 9:  "A" Annual forecast for current year
    % col 10: "B" Annual forecast for next year
    % col 11: "C" Annual forecast for year after next (not available for all variables)
    % col 12: "D" Annual forecast for three years ahead (not available for all variables)
    % convert -999 into NaN
    SPFimport.data(SPFimport.data == -999) = NaN;
    
    % check dates
    Y = SPFimport.data(:,1);
    if ~isequal(Y, year(dates))
        error('inconsistent date vector in the SPF file (years)');
    end
    Q = SPFimport.data(:,2);
    if ~isequal(Q, quarter(dates))
        error('inconsistent date vector in the SPF file (quarters)');
    end
    
    %% prepare Ybar parameters
    maxHorizon    = 4; % fixed horizon data
    Nbar          = size(SPFimport.data, 2) - 9; % only using next year and the years thereafter (columns B, C, and D)
    maxBarHorizon = 3 + Nbar * 4; % maximal forecast horizon in Q1 has three quarters in current year, and four quarters in subsequent years
    
    %% collect forecast data
    if doNIPA
        
        % convert SPF level forecasts for NIPA variables into forecasts of
        % *quarterly* changes
        Yforecast            = 400 * (log(SPFimport.data(:,4:8)) - log(SPFimport.data(:,3:7)));
        
        Ycurrentlevelforecast = SPFimport.data(:,9); % note: for now just the level, pulling rest later
        
        Ybarforecast = 100 * (log(SPFimport.data(:,10:end)) - log(SPFimport.data(:,9:end-1))); % for B: using A as basis (not the quarterly forecasts)
        checkdiff(Ybarforecast, 100 * (log(SPFimport.data(:,10:end) ./ SPFimport.data(:,9:end-1))));
        
    else
        % pull out SPF forecasts for quarter 0 through 4
        Yforecast        = SPFimport.data(:,4:8);
        Ybarforecast     = SPFimport.data(:,10:end);
        Ycurrentforecast = SPFimport.data(:,9);
    end
    
    
    if size(Ybarforecast) ~= Nbar
        error houston
    end
    
    %% collect measures of realized data
    Yrealized  = NaN(T+1,1); % includes one lag
    Yvintlags  = NaN(T+1,4);
    Yadvlags   = NaN(T+1,4);
    
    if doNIPA
        % starting with the 1968Q4 vintage (first SPF), collect 1st releases (so first one is 1968Q3)
        first_SPF_date_idx = find(RTdata.vindates==datetime('1968Q4','InputFormat','yyyyQQQ'),1,'first');
        last_SPF_date_idx  = find(RTdata.vindates==datetime(SPFsamstop,'ConvertFrom','datenum'),1,'first');
        first_obs_idx      = find(RTdata.obsdates==datetime(realizedDataDates(1),'ConvertFrom','datenum'),1,'first');
        
        
        % walk along diagonal of RTDSM sheet to collect first-release data
        NobsRT = size(RTdata.data,1);
        ii = 0;
        for vin = first_SPF_date_idx : last_SPF_date_idx
            % check we are on the diagonal of the RTDSM sheet
            if (first_obs_idx+ii < NobsRT) && ~isnan(RTdata.data(first_obs_idx+ii+1,vin))
                warning('we are not collecting Yrealized on the RTDSM diagonal')
            end
            % collect data
            Yrealized(ii+1)   = RTdata.data(first_obs_idx+ii,vin);
            Yvintlags(ii+1,:) = RTdata.data(first_obs_idx+ii+(-3:0),vin);
            % prepare next iteration
            ii              = ii + 1;
        end
        
        for lag = 1 : 4
            for t = lag : T
                Yadvlags(t,4-lag+1) = Yrealized(t-lag+1);
            end
        end
    else
        %% load realized data from current FRED (not using RTDSM data, since release of quarterly vintage may lag our data updates)
        
        switch upper(datalabel)
            case 'UNRATE'
                FRED      = fredreadcsv('UNRATE',[], 'm', 'q', 'avg');
                Yrealized = fredexpand(FRED, realizedDataDates);
            case {'CPI'}
                FRED      = fredreadcsv('CPIAUCSL',[], 'm', 'q', 'avg');
                cpiLevel  = fredexpand(FRED, realizedDataDates);
                Yrealized = ((cpiLevel(2:end) ./ cpiLevel(1:end-1)).^4 - 1) * 100;
            otherwise
                error('datalabel %s not yet supported', upper(datalabel))
        end
        
        
        for lag = 1 : 4
            for t = lag : T
                Yadvlags(t,4-lag+1) = Yrealized(t-lag+1);
            end
        end
        Yvintlags = Yadvlags;
    end
    
    %% prepare SPFpath
    % reconstruct forecastB
    yy = 1;
    fcstB      = Ybarforecast(:,yy);
    fcstBpred  = NaN(T,1); % "pred" for predetermined
    fcstBpred1 = NaN(T,1); % computed based on first-release lags
    fcstBwSum    = NaN(T,1);
    fcstBw0num = NaN(T,1);
    % SPF path collects nowcast, the last four lags of data, plus all quarterly forecasts, and three empty horizons
    SPFpath    = cat(2, Yvintlags(1:end-1,:), Yforecast, zeros(T,3));
    SPFpath1   = cat(2, Yadvlags(1:end-1,:), Yforecast, zeros(T,3));
    
    Nspf       = size(SPFpath,2);
    if doNIPA
        tentshape  = [1 2 3 4 3 2 1] / 16; % used only if doNIPA is true
        for t = 1 : T
            ndxTent          = (4 - Q(t)) + yy * 4 + (-6:0);
            ndxTent          = ndxTent + 1 + 4; % for nowcast and 4 lagged realized values
            weights          = zeros(Nspf,1);
            weights(ndxTent) = tentshape;
            fcstBpred(t)     = SPFpath(t,:) * weights;
            fcstBpred1(t)    = SPFpath1(t,:) * weights;
            fcstBwSum(t)       = sum(weights(10:end)); % sum over the weights attached to values beyond quarterly SPF
            fcstBw0num(t)    = sum(weights(10:end) > 0);
        end
    else
        for t = 1 : T
            ndxMA            = (4 - Q(t)) + yy * 4 + (-3:0);
            ndxMA            = ndxMA + 1 + 4; % for nowcast and 4 lagged realized values
            weights          = zeros(Nspf,1);
            weights(ndxMA)   = 1/4;
            fcstBpred(t)     = SPFpath(t,:) * weights;
            fcstBpred1(t)    = SPFpath1(t,:) * weights;
            fcstBwSum(t)     = sum(weights(10:end)); % sum over the weights attached to values beyond quarterly SPF
        end
    end
    
    firstobsB = find(~isnan(fcstB), 1);
    ndxBobs   = transpose(ismember(1:T,firstobsB-1:T)); % logical col index for future use

    %% time series plots of fcstBpred and fcstB
    % thisfig = figure;
    % hold on
    % hPred = plot(dates(ndxBobs), fcstBpred(ndxBobs), 'color', colors4plots("lightblue"), 'LineWidth', 2);
    % hB    = plot(dates(ndxBobs), fcstB(ndxBobs), 'color', colors4plots("darkblue"), 'LineWidth', 2);
    % hQ4 = xline(dates(quarter(dates) == 4), 'k:');
    % legend([hB hPred hQ4(1)], 'SPF', 'predetermined', 'Q4', 'Location', 'northwest')
    % xtickdates(dates(ndxBobs));
    % title(sprintf('SPF and Predetermined part of FcstB - %s', datalabel))
    % wrapthisfigure(thisfig, sprintf('fcstB-%s', datalabel), wrap)
    
    %% compute deviation between actual and predetermined fcstB
    fcstBdelta = fcstB - fcstBpred;
    
    %% plot fcstBdelta separately for Q4
    thisQ = 4;
    thisfig = figure;
    hold on
    set(gca, 'FontSize', fontsize)
    jack = fcstBdelta;
    jack(quarter(dates) ~= thisQ) = NaN;
    hDelta = bar(dates, jack, 'facecolor', colors4plots("darkblue"), 'edgecolor', colors4plots("darkblue"), 'LineWidth', 2);
    xtickdates(dates(firstobsB-1:end));
    if ~strcmpi(datalabel, 'CPI')
        ylim([-.3 .3])
    else
        ylim([-.4 1.2])
    end
    wrapthisfigure(thisfig, sprintf('fcstBdeltaQ%d-%s', thisQ, datalabel), wrap)
    % title(sprintf('%s: Discrepancy between SPF next year and quarterly forecasts (Q%d only)', datalabel, q))
    % wrapthisfigure(thisfig, sprintf('fcstBdeltaQ%d-%s-WITHTITLE', q, datalabel), wrap)
    
    %% plot imputed values and data
    fcstBweighteddelta = fcstBdelta ./ fcstBwSum;
    fcstBweighteddeltaQ = NaN(T, 4);
    for q = 1:4
        fcstBweighteddeltaQ(quarter(dates) == q, q) = fcstBweighteddelta(quarter(dates) == q);
    end
    markertypes = {'o', 's', 'd'};
    thisfig = figure;
    set(gca, 'FontSize', fontsize)
    hold on
    hanni = NaN(3,1);
    for q=1:3
        ndxQ = quarter(dates) == q;
        hanni(q) = plot(dates(ndxQ), fcstBweighteddeltaQ(ndxQ,q), markertypes{q}, 'color', colors4plots(q), ...
            'MarkerSize', 10, 'LineWidth', 2);
    end
    hRealized   = plot(dates, Yrealized(2:end), ':', 'color', colors4plots("black"), 'LineWidth', 1);
    hfcstB      = plot(dates, fcstB, '-', 'color', colors4plots("black"), 'LineWidth', 1);
    xtickdates(dates(firstobsB:end));
    if strcmpi(datalabel, 'RGDP')
        ylim([-10 10]) % hide COVID outliers
    end
    wrapthisfigure(thisfig, sprintf('fcstBImputation-%s', datalabel), wrap, [], [], [], [], true);
    legend([hanni; hfcstB; hRealized], 'Q1', 'Q2', 'Q3', 'SPF next-year', 'Realized', ...
        'NumColumns', 2, 'Location', 'best')
    wrapthisfigure(thisfig, sprintf('fcstBImputation-%s-WITHLEGEND', datalabel), wrap)

    % redo previous plot, but show only data pre-2020
    ndxPreCOVID = dates < datenum(2020,1,1);
    T2020       = find(dates == datenum(2020,1,1));
    thisfig = figure;
    set(gca, 'FontSize', fontsize)
    hold on
    hanni = NaN(3,1);
    for q=1:3
        ndxQ = quarter(dates) == q & ndxPreCOVID;
        hanni(q) = plot(dates(ndxQ), fcstBweighteddeltaQ(ndxQ,q), markertypes{q}, 'color', colors4plots(q), ...
            'MarkerSize', 10, 'LineWidth', 2);
    end
    thisRealized = Yrealized(2:end);
    hRealized   = plot(dates(ndxPreCOVID), thisRealized(ndxPreCOVID), ':', 'color', colors4plots("black"), 'LineWidth', 1);
    hfcstB      = plot(dates(ndxPreCOVID), fcstB(ndxPreCOVID), '-', 'color', colors4plots("black"), 'LineWidth', 1);
    xtickdates(dates(firstobsB:T2020));
    wrapthisfigure(thisfig, sprintf('fcstBImputation-%s-PreCOVID', datalabel), wrap, [], [], [], [], true);
    legend([hanni; hfcstB; hRealized], 'Q1', 'Q2', 'Q3', 'SPF next-year', 'Realized', ...
        'NumColumns', 2, 'Location', 'best')
    wrapthisfigure(thisfig, sprintf('fcstBImputation-%s-PreCOVID-WITHLEGEND', datalabel), wrap)

    %% plot bar chart of fcstBweighteddelta separately for Q1, Q2, Q3
    fcstBweighteddelta = fcstBdelta ./ fcstBwSum;
    % Collect fcstBweighteddelta in separate columns for each quarter
    fcstBweighteddeltaQ = NaN(T, 4);
    for q = 1:4
        fcstBweighteddeltaQ(quarter(dates) == q, q) = fcstBweighteddelta(quarter(dates) == q);
    end
    
    thisfig = figure;
    set(gca, 'FontSize', fontsize)
    ndxBobsX4 = ndxBobs & quarter(dates) ~= 4;
    hold on
    hanni = NaN(3,1);
    for q = 1 : 3
        hanni(q) = bar(dates(ndxBobsX4), fcstBweighteddeltaQ(ndxBobsX4,q), 1, 'EdgeColor','flat');
    end
    xticks(datenum(1980:5:2025,1,1))
    datetick('x', 'keepticks')
    orient landscape
    wrapthisfigure(thisfig, sprintf('fcstBnonPred-%s', datalabel), wrap, [], [], [], [], true);
    legend(hanni, 'Q1', 'Q2', 'Q3', 'Location', 'northeast')
    wrapthisfigure(thisfig, sprintf('fcstBnonPred-%s-WITHLEGEND', datalabel), wrap)
    % title(sprintf('%s: Non-predetermined part of FcstB (weighted)', datalabel))
    % wrapthisfigure(thisfig, sprintf('nonPredfcstBdeltad-%s-WITHTITLE', datalabel), wrap)
    
    
end

%% finish
dockAllFigures
finishwrap
