%% construct data set for project "kensington"

% USE getKensingtonData.ipynb to collect latest raw data
% from SPF, RTDSM and FRED

% run collectRTData first

% PLEASE ADAPT start and end dates of sample (three cells below)

%% load toolboxes
path(pathdef)
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

Kavg = 4; % for definition of YAVG

%% select choice of data
DATALABELS = {'UNRATE', 'CPI', 'PGDP', 'RGDP'}; 

for d = 1 : length(DATALABELS)

    datalabel = DATALABELS{d};

    doNIPA = any(ismember({'PGDP', 'RGDP'}, datalabel));

    %% load RealTime Data
    RTdata = load(sprintf('kensington%sdataRT', datalabel));

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
    Yvintlags  = NaN(T+1,4); % track four lags of current vintage
    Yadvlags   = NaN(T+1,4); % track four lags of advance-release vintages

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


    end

    %% construct Zdata and loading C
    Ny         = maxBarHorizon + 1 + 1; % adding nowcast and lagged data

    % PART1: fixed horizon
    Nz1        = maxHorizon + 1 + 1; % adding nowcast and lagged data
    Zdata1     = cat(2, Yrealized(1:end-1), Yforecast);
    C1         = repmat(eye(Nz1,Ny), [1 1 T]);
    Z1label    = cell(Nz1,1);
    % Z1label{1} = 'data(-1)';
    for h = -1 : maxHorizon
        Z1label{h+2} = sprintf('h =%2d', h);
    end

    % PART2: fixed event
    Nz2        = Nbar;
    Zdata2     = Ybarforecast;
    Zdata2V1   = Ybarforecast; % earlier version: Q4 obs of next year set to NaN
    C2         = zeros(Nz2,Ny,T);
    tentshape  = [1 2 3 4 3 2 1] / 16; % used only if doNIPA is true
    for t = 1 : T
        for ii = 1 : Nz2
            % ndx jump off counts quarters left in current
            if doNIPA
                ndxTent          = (4 - Q(t)) + ii * 4 + (-6:0);
                ndxTent          = ndxTent + 2; % for nowcast and lagged realized
                if any(ndxTent < 1)
                    if ii == 1 && Q(t) == 4 % amend Bforecast in Q4 to account for Q2 lagged value (which is outside state space)
                        Zdata2V1(t,ii) = NaN; % C2 can be left unaffected for this case, will be set to 0 anyway on account of NaN value
                        if find(ndxTent < 1) ~= 1
                            error('ndxTent off')
                        end
                        C2(ii,ndxTent(ndxTent > 0),t) = tentshape(ndxTent > 0) / sum(tentshape(ndxTent > 0));
                        if ~isnan(Zdata2(t,ii))
                            % take out the Q2 value, and rescale the observation to match (remaining) weights that sum to one
                            Q2realized   = Yvintlags(t-1,3);
                            Zdata2(t,ii) = (Zdata2(t,ii) - tentshape(1) * Q2realized) / sum(tentshape(ndxTent > 0)); 
                        end
                    else
                        error('ndxTent off')
                    end
                else
                    C2(ii,ndxTent,t) = tentshape;
                end
            else
                ndxMA4          = (4 - Q(t)) + ii * 4 + (-3:0);
                ndxMA4          = ndxMA4 + 2; % for nowcast and lagged realized
                C2(ii,ndxMA4,t) = 0.25; % .25 = 1/4
            end
        end
    end
    Z2label    = cell(Nz2,1);
    for n = 1 : Nz2
        Z2label{n} = sprintf('y =%2d', n);
    end

    % PART 3: put all together
    Nz         = Nz1 + Nz2;
    Zdata      = cat(2, Zdata1, Zdata2);
    ZdataV1    = cat(2, Zdata1, Zdata2V1);
    ZnannyV1   = isnan(ZdataV1);
    Cz         = cat(1, C1, C2);
    Zlabel     = cat(1, Z1label, Z2label);

    % construct matrix of "actual forecasts" (includes *All* next year readings, and omits log transformation for NIPA
    YactualForecasts = cat(2, Yforecast, Ybarforecast);
    if doNIPA
        YactualForecasts = (exp(YactualForecasts / 100) - 1) * 100;
    end

    % PART 4: handle NaN
    Znanny        = isnan(Zdata);
    %     Zdata(Znanny) = 0; leave as NaN
    for t = 1 : T
        Cz(Znanny(t,:),:,t) = 0;
    end

    % sanity check 2
    for t = 1 : T
        if rank(Cz(~Znanny(t,:),:,t)) == 0
            warning('zero Rank C at t=%d, d=%d', t, d)
        end
        % if abs(det(Cz(~Znanny(t,:),:,t)*Cz(~Znanny(t,:),:,t)')) < 1e-6
        %     warning('singular C at t=%d, d=%d', t, d)
        % end
    end

    % sanity check 3: all rows sum to one
    sumCz      = sum(Cz,2);
    if ~all(sumCz(~transpose(Znanny)) == 1)
        error('Cz gap loadings should sum to one')
    end

    %% construct future realized values (of quarterly fixed horizon forecasts)
    Nhorizons = Ny; % was: Ny-1, but setting it to Ny allows comparison with trend horizon
    Nhorizons = 17; % was: Ny-1, but setting it to Ny allows comparison with trend horizon
    Yfuture   = NaN(T,Nhorizons); % realized values for all Ny-1 forecast horizons

    % we use 2nd release as realizations for PGDP, RGDP, while latest vintage for UNRATE, CPI


    if doNIPA

        release_idx  = 2;
        first_RT_vin = find(RTdata.vindates==datetime('1968Q4','InputFormat','yyyyQQQ'),1,'first')+release_idx;

        for t=1:T % looping over SPF rounds
            thisSPF_idx = find(RTdata.obsdates==datetime('1968Q4','InputFormat','yyyyQQQ'),1,'first')+t-1;

            % when t=1, thisSPF_idx is in the row of 1968Q4 data
            for n=1:Nhorizons % looping over horizons, n=1 is nowcast
                % Yfuture(t,n) contains the 2nd release of data in quarter t+n-1
                % when t=1 and n=1, we load the 1968Q4 RGDP/PGDP growth according to 1969Q2 vintage (2nd release)
                act_vin_idx = first_RT_vin+t-1+n-1;
                if act_vin_idx <= size(RTdata.data,2)
                    Yfuture(t,n) = RTdata.data(thisSPF_idx+n-1,act_vin_idx);
                end
            end
        end


    else

        % drop lag from Yrealized
        yrealized = Yrealized(2:end);

        for t = 1 : T
            for n = 1 : Nhorizons
                h = n - 1; % first horizon is the nowcast
                if t + h <= T
                    Yfuture(t,n) = yrealized(t+h);
                end
            end
        end

    end

    %% YCUMfuture
    YCUMfuture = cumsum(Yfuture, 2);

    %% construct future realized values (of calendar-year fixed-event forecasts)
    YBARfuture     = NaN(T,Nbar);
    YCURRENTfuture = NaN(T,1);


    % first SPF round is 1968Q4, in YBARfuture we have outcomes from 1969, 1970,...
    % first SPF round is 1968Q4, in YCURRENTfuture we have outcomes from 1968, 1969,...
    leads_vec = 0:-1:-Nbar+1;


    if strcmpi(datalabel, 'CPI')
        % TODO: need to construct Q4/Q4 RTdata for CPI
        warning('RTdata for Q4/Q4 CPI yet to be implemented')
    else
        if doNIPA

            firstYBAR_idx     = find(RTdata.calY_gr_dates==(year(dates(1))+1));
            firstYCURRENT_idx = find(RTdata.calY_gr_dates==(year(dates(1))));

            tmp_BAR_1 = RTdata.calY_gr(firstYBAR_idx:end);
            tmp_BAR   = NaN(size(tmp_BAR_1,1),length(leads_vec));
            l_idx     = 1;

            for l = leads_vec

                tmp_BAR(1:end-abs(l),l_idx) = tmp_BAR_1(l_idx:end,1);

                l_idx = l_idx + 1;
            end

            tmp_CURRENT = RTdata.calY_gr(firstYCURRENT_idx:end);

        else % UNRATE or CPI



            firstYBAR_idx     = find(RTdata.calY_dates==(year(dates(1))+1));
            firstYCURRENT_idx = find(RTdata.calY_dates==(year(dates(1))));


            tmp_BAR_1 = RTdata.calY_avg(firstYBAR_idx:end);
            tmp_BAR   = NaN(size(tmp_BAR_1,1),length(leads_vec));
            l_idx     = 1;

            for l = leads_vec

                tmp_BAR(1:end-abs(l),l_idx) = tmp_BAR_1(l_idx:end,1);

                l_idx = l_idx + 1;
            end


            tmp_CURRENT = RTdata.calY_avg(firstYCURRENT_idx:end);
        end

        tmp2 = repelem(tmp_BAR,4,1);
        tmp2 = tmp2(4:end,:); % in 1968, there was only one SPF round
        YBARfuture(1:size(tmp2,1),:) = tmp2;
        if size(YBARfuture, 1) > T % can happen if more realized values are available than surveys
            YBARfuture = YBARfuture(1:T,:);
        end

        tmp2 = repelem(tmp_CURRENT(~isnan(tmp_CURRENT)),4,1);
        tmp2 = tmp2(4:end,:); % in 1968, there was only one SPF round
        YCURRENTfuture(1:size(tmp2,1),1) = tmp2;

        if size(YCURRENTfuture, 1) > T % can happen if more realized values are available than surveys
            YCURRENTfuture = YCURRENTfuture(1:T,:);
        end
    end

    clear tmp*

    %% construct future realized value of YAVG
    % collect sequence of historical realizations to complete the current year

    ylags = NaN(T,Kavg-1);
    for t = Kavg + 1 : T
        thisDate      = dates(t);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI'}
                % use latest vintage, not real-time data
                RT_vintage_idx = length(RTdata.vindates); % pick last vintage
            case {'RGDP', 'PGDP'}
                RT_vintage_idx = find(datenum(RTdata.vindates) == thisDate,1,'first'); % vintage, time of SPF round
            otherwise
                RT_vintage_idx = []; %#ok<NASGU> % to avoid parfor warning
                error('datalabel <<%s>> not supported', datalabel)
        end
        % last observation available at time of SPF round
        RT_obsdate_idx = find(datenum(RTdata.obsdates)==datenum(dateshift(datetime(thisDateLabel,'InputFormat','yyyyQQQ'),'start','quarter',-1)),1,'first');
        % vector of vintage data
        RT_vec = RTdata.data(1:RT_obsdate_idx,RT_vintage_idx);

        % collect three lags, note that RT_obsdate_idx is one lag realtive to t
        ylags(t,:) = RT_vec(end-(Kavg-1)+1:end); % store lag three first
    end

    % construct YAVGyfuture
    YAVGfuture = NaN(size(Yfuture));
    % h < Kavg
    for h = 1 : (Kavg - 1)
        YAVGfuture(:,h) = (sum(Yfuture(:,1:h), 2) + sum(ylags(:,end-(Kavg-1)+h:end), 2)) / Kavg;
    end
    % h >= Kavg
    for h = Kavg : size(YAVGfuture, 2)
        YAVGfuture(:,h) = sum(Yfuture(:,h+(-(Kavg-1):0)), 2) / Kavg;
    end

    %% patch Ycurrentforecast
    if doNIPA
        % collect sequence of historical realizations to construct previous year's average level
        ylags = NaN(T,1);
        for t = 1 : T
            thisDate      = dates(t);
            thisDateLabel = datestr(thisDate, 'yyyyqq');
            RT_vintage_idx = find(datenum(RTdata.vindates) == thisDate,1,'first'); % vintage, time of SPF round

            prevYear      = year(thisDate) - 1;
            RT_prevYear_idx = year(RTdata.obsdates) == prevYear;
            if sum(RT_prevYear_idx) == 4 % take only full years

                % vector of vintage data
                RT_vec   = RTdata.leveldata(RT_prevYear_idx,RT_vintage_idx);
                ylags(t) = mean(RT_vec); % store lag three first

            end
        end

        Ycurrentforecast = 100 * (Ycurrentlevelforecast ./ ylags - 1);

    end

    %% construct Ylabel
    Ylabel = cell(Ny,1);
    Ylabel{1} = 'y(-1)';
    for n = 2 : Ny
        thisHorizon = n - 2;
        Ylabel{n} = sprintf('h=%d', thisHorizon);
    end

    YBARlabel = cell(Nbar,1);
    for n = 1 : Nbar
        YBARlabel{n} = sprintf('year=%d', n);
    end

    Nforecasts     = 1 + maxHorizon + Nbar;
    if ~isequal(Nforecasts, Nz-1)
        error('something is off')
    end
    YforecastLabel = cell(Nforecasts,1);
    for n = 1 : maxHorizon + 1
        thisHorizon = n - 1;
        YforecastLabel{n} = sprintf('quarter=%d', thisHorizon);
    end
    YforecastLabel(maxHorizon+2:end) = YBARlabel;

    %% store data
    matfile = sprintf('kensington%sdata', upper(datalabel));

    datesQ = Q;
    if doNIPA
        Czhat   = C2;
        Czbar   = [];
    else
        Czhat   = [];
        Czbar   = C2;
    end
    save(matfile, 'Ny', 'Nz', 'T', 'dates', 'datesQ', 'doNIPA', ...
        'Yfuture', 'Ylabel', 'Nbar', 'YBARlabel', ...
        'YactualForecasts', 'Nforecasts', 'YforecastLabel', ...
        'Ycurrentforecast', 'Ybarforecast', ...
        'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czbar', 'Czhat', ...
        'ZdataV1', 'ZnannyV1', ... % legacy code
        'YBARfuture', 'YCURRENTfuture', ...
        'YCUMfuture', ...
        'YAVGfuture', 'Kavg', ...
        '-v7.3')

end