%% construct data set for project "kensington"

% USE getKensingtonData.ipynb to collect latest raw data
% from SPF, RTDSM and FRED

% PLEASE ADAPT start and end dates of sample a few cells below

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/

%#ok<*NANVAR>
%#ok<*UNRCH>

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% ADAPT END OF SAMPLE HERE


SPFsamstop  = datenum(2022,4,1);  % 2022Q1

% dates follow the FRED convention of placing timestamp at beginning of period

%% select choice of data

datalabel   = 'RGDP10';
datalabelRT = 'RGDP';

doNIPA = true;

%% prepare variables

% quarters are time-stamped on the first day of the quarter (FRED convention)
% note: the following could also be automated, but might be good to set a
% few things manually (also forces us to check things when updating the data)

dates     = genrQdates(1968,year(SPFsamstop));
dates     = dates(dates <= SPFsamstop);
% prepare dates for realized data
% realizedDataDates points to date vector starting at first lag of realized Y till last observation (T+1 observations)
realizedDataDates = dates(3:end);
dates      = dates(4:end); % SPF dates begin only in 1968Q4
T          = length(dates);
maxHorizon = 4; % fixed horizon data


%% load SPF responses for the level of PGDP and convert into inflation rates
% the results of this cell are to be stored in Ydata(:,2:6)

SPFimport   = importdata('Mean_RGDP_Level.xlsx'); % this is the original xls file converted into csv

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
NbarAnnual          = size(SPFimport.data, 2) - 9; % only using next year and the years thereafter (columns B, C, and D)
maxBarAnnualHorizon = 3 + NbarAnnual * 4; % maximal forecast horizon in Q1 has three quarters in current year, and four quarters in subsequent years

%% process SPFimport10
SPF10import = importdata('Mean_RGDP10_Level.xlsx'); % this is the original xls file converted into csv

SPF10import.data(SPF10import.data == -999) = NaN;

% check dates
Y = SPF10import.data(:,1);
if ~isequal(Y, year(dates))
    error('inconsistent date vector in the SPF file (years)');
end
Q = SPF10import.data(:,2);
if ~isequal(Q, quarter(dates))
    error('inconsistent date vector in the SPF file (quarters)');
end

Nbar          = NbarAnnual + 1;
maxBarHorizon = 39; % nowcast 0 plus 39

%% collect forecast data

% convert SPF level forecasts for NIPA variables into forecasts of
% *quarterly* changes
Yforecast = 400 * (log(SPFimport.data(:,4:8)) - log(SPFimport.data(:,3:7)));

YbarforecastAnnual = 100 * (log(SPFimport.data(:,10:end)) - log(SPFimport.data(:,9:end-1))); % for B: using A as basis (not the quarterly forecasts)
checkdiff(YbarforecastAnnual, 100 * (log(SPFimport.data(:,10:end) ./ SPFimport.data(:,9:end-1))));

Ybarforecast10 = log(1 + SPF10import.data(:,3) ./ 100) .* 100;
% Ybarforecast   = [YbarforecastAnnual Ybarforecast10];

% if size(Ybarforecast) ~= Nbar
%     error houston
% end



%% collect measures of realized data
Yrealized  = NaN(T+1,1); % includes one lag

%%% load second revision data *of growth rates* and use it as measure of realized outcome

RTdata = load(sprintf('kensington%sdataRT', datalabelRT));

% starting with the 1968Q4 vintage (first SPF), collect 1st releases (so first one is 1968Q3)
first_SPF_date_idx = find(RTdata.vindates==datetime('1968Q4','InputFormat','yyyyQQQ'),1,'first');
last_SPF_date_idx  = find(RTdata.vindates==datetime(SPFsamstop,'ConvertFrom','datenum'),1,'first');
first_obs_idx      = find(RTdata.obsdates==datetime(realizedDataDates(1),'ConvertFrom','datenum'),1,'first');


idx_tmp    = 1;

for vin = first_SPF_date_idx : last_SPF_date_idx
    Yrealized(idx_tmp)   = RTdata.data(first_obs_idx+idx_tmp-1,vin);
    idx_tmp              = idx_tmp+1;
end

%% construct nowcast for current year
YbarNowcastAnnualLevel  = SPFimport.data(:,9);
YbarLastYearAnnualLevel = NaN(T,1);
T1 = find(~isnan(Ybarforecast10), 1, 'first');
for t = T1 : T
    if ~isnan(Ybarforecast10(t))
        if Q(t) ~= 1
            error('RGDPD10 reading found outside Q1')
        end
        RT_vintage_idx = find(datenum(RTdata.vindates) == dates(t),1,'first'); % vintage, time of SPF round

        lastYearQ4                 = find(datenum(RTdata.obsdates) == dates(t)) - 1;
        YbarLastYearAnnualLevel(t) = mean(RTdata.leveldata(lastYearQ4 + (-3:0), RT_vintage_idx)); % taking mean to correspond to SPF nowcast
    end
end

YHATcurrentyear = 100 * (log(YbarNowcastAnnualLevel) - log(YbarLastYearAnnualLevel));
Ybarforecast9   = (10 * Ybarforecast10 - YHATcurrentyear) / 9;

figure
hold on 
plot(dates(T1:T), Ybarforecast10(T1:T), 'x')
plot(dates(T1:T), Ybarforecast9(T1:T), 'o')
xtickdates(dates(T1:T))
legend('RGDP10', 'RGDP10 - 10/9 Ycurrent')

% collection of RAW data is now complete

%% construct Zdata and loading C
Ny         = maxBarHorizon + 1 + 1; % adding nowcast and lagged data

% PART1: fixed horizon
Nz1        = maxHorizon + 1 + 1; % adding nowcast and lagged data
Zdata1     = cat(2, Yrealized(1:end-1), Yforecast);
C1         = repmat(eye(Nz1,Ny), [1 1 T]);
Z1label    = cell(Nz1,1);
Z1label{1} = 'y(-1)';
for n = 1 : maxHorizon + 1
    Z1label{n+1} = sprintf('h%d', n-1);
end

% PART2: fixed event calendar years
Nz2        = NbarAnnual;
Zdata2     = YbarforecastAnnual;
C2         = zeros(Nz2,Ny,T);
tentshape  = [1 2 3 4 3 2 1] / 16;
for t = 1 : T
    for ii = 1 : Nz2
        if ii == 1 && Q(t) == 4 % ignore B foreast at Q4 (since spanned by the quarterly forecasts)
            Zdata2(t,ii) = NaN; % C2 left as zero
        else
            % ndx jump off counts quarters left in current
            ndxTent          = (4 - Q(t)) + ii * 4 + (-6:0);
            ndxTent          = ndxTent + 2; % for nowcast and lagged realized
            if any(ndxTent < 1)
                error('ndxTent off')
            end
            C2(ii,ndxTent,t) = tentshape;
        end
    end
end
Z2label    = cell(Nz2,1);
for n = 1 : Nz2
    Z2label{n} = sprintf('hbar-%d', n);
end

% PART3: add 10-year forecast
Nz3        = 1;
Zdata3     = Ybarforecast9;
Z3label    = 'hbar-40';
C3         = zeros(Nz3,Ny,T);
tentshape            = repmat(4, 1, 39);
tentshape(1:3)       = 1 : 3;
tentshape(end-2:end) = 3: -1 : 1;
tentshape            = tentshape ./ (4 * 4 * 9);
for t = 1 : T
    if Q(t) == 1 && ~isnan(Zdata3(t))
        % ndx jump off counts quarters left in current
        ndxTent          = 1 : 39;
        ndxTent          = ndxTent + 2; % for nowcast and lagged realized
        C3(1,ndxTent,t)  = tentshape;
    end
end

% PART 4: put all together
Nz         = Nz1 + Nz2 + Nz3;
Zdata      = cat(2, Zdata1, Zdata2, Zdata3);
Cz         = cat(1, C1, C2, C3);
Zlabel     = cat(1, Z1label, Z2label, Z3label);

% construct matrix of "actual forecasts" (includes *All* next year
% readigns, and omits log transformation for NIPA
YactualForecasts = cat(2, Yforecast, YbarforecastAnnual);
YactualForecasts = (exp(YactualForecasts / 100) - 1) * 100;

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
    if abs(det(Cz(~Znanny(t,:),:,t)*Cz(~Znanny(t,:),:,t)')) < 1e-6
        warning('singular C at t=%d, d=%d', t, d)
    end
end


%% construct future realized values (of quarterly fixed horizon forecasts)
Nhorizons = Ny-1;
Yfuture   = NaN(T,Nhorizons); % realized values for all Ny-1 forecast horizons

% we use 2nd release as realizations for PGDP, RGDP, while latest vintage for UNRATE, CPI, TBILL



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



%% construct future realized values (of calendar-year fixed-event forecasts)
YBARfuture = NaN(T,NbarAnnual);
YCURRENTfuture = NaN(T,1);


RTdata = load(sprintf('kensington%sdataRT', datalabelRT));

% first SPF round is 1968Q4, in YBARfuture we have outcomes from 1969, 1970,...
% first SPF round is 1968Q4, in YCURRENTfuture we have outcomes from 1968, 1969,...
leads_vec = 0:-1:-NbarAnnual+1;


firstYBAR_idx = find(RTdata.calY_gr_dates==(year(dates(1))+1));
firstYCURRENT_idx = find(RTdata.calY_gr_dates==(year(dates(1))));

tmp_BAR_1 = RTdata.calY_gr(firstYBAR_idx:end);
tmp_BAR = NaN(size(tmp_BAR_1,1),length(leads_vec));
l_idx = 1;

for l = leads_vec

    tmp_BAR(1:end-abs(l),l_idx) = tmp_BAR_1(l_idx:end,1);

    l_idx = l_idx + 1;
end

tmp_CURRENT = RTdata.calY_gr(firstYCURRENT_idx:end);



tmp2 = repelem(tmp_BAR,4,1);
tmp2 = tmp2(4:end,:); % in 1968, there was only one SPF round
YBARfuture(1:size(tmp2,1),:) = tmp2;

tmp2 = repelem(tmp_CURRENT(~isnan(tmp_CURRENT)),4,1);
tmp2 = tmp2(4:end,:); % in 1968, there was only one SPF round
YCURRENTfuture(1:size(tmp2,1),1) = tmp2;




%% construct Ylabel
Ylabel = cell(Ny,1);
Ylabel{1} = 'y(-1)';
for n = 2 : Ny
    thisHorizon = n - 2;
    Ylabel{n} = sprintf('h=%d', thisHorizon);
end

YBARlabel = cell(NbarAnnual,1);
for n = 1 : NbarAnnual
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
YforecastLabel(maxHorizon+2:end-1) = YBARlabel;
YforecastLabel{end} = '10-year';

%% store data
matfile = sprintf('kensington%sdata', upper(datalabel));

datesQ  = Q;
Czhat   = C2; % Add C3?
Czbar   = [];
save(matfile, 'Ny', 'Nz', 'T', 'dates', 'datesQ', 'doNIPA', ...
    'Yfuture', 'Ylabel', 'Nbar',  'YBARlabel', ...
    'YactualForecasts', 'Nforecasts', 'YforecastLabel', ...
    'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czbar', 'Czhat', ...
    'YBARfuture', 'YCURRENTfuture', ...
    '-v7.3')
