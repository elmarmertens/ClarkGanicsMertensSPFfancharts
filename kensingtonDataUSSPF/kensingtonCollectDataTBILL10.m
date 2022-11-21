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

datalabel   = 'TBILL10';
datalabel0  = 'TBILL';
doNIPA      = false;

%% prepare variables

% quarters are time-stamped on the first day of the quarter (FRED convention)
% note: the following could also be automated, but might be good to set a
% few things manually (also forces us to check things when updating the data)

dates     = genrQdates(1968,year(SPFsamstop));
dates     = dates(dates <= SPFsamstop);
% prepare dates for realized data
% realizedDataDates points to date vector starting at first lag of realized Y till last observation (T+1 observations)
switch upper(datalabel0)
    case {'TBILL'}
        realizedDataDates = dates(3:end);     % include 68Q3 to grab lagged values
    case {'CPI'}
        realizedDataDates = dates(2:end);     % include 68Q2 to grab lagged *CPI* level for construction of lagged inflation data
    otherwise
        error('datalabel %s not yet supported', upper(datalabel))
end
dates      = dates(4:end); % SPF dates begin only in 1968Q4
T          = length(dates);
maxHorizon = 4; % fixed horizon data


%% load SPF responses for the level of PGDP and convert into inflation rates
% the results of this cell are to be stored in Ydata(:,2:6)

switch upper(datalabel0)
    case 'TBILL'
        SPFimport = importdata('Mean_TBILL_Level.xlsx'); % this is the original xls file converted into csv
    case {'CPI'}
        SPFimport = importdata('Mean_CPI_Level.xlsx'); % this is the original xls file converted into csv
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
NbarAnnual          = size(SPFimport.data, 2) - 9; % only using next year and the years thereafter (columns B, C, and D)
maxBarAnnualHorizon = 3 + NbarAnnual * 4; % maximal forecast horizon in Q1 has three quarters in current year, and four quarters in subsequent years

%% process SPFimport10
SPF10import = importdata('Mean_BILL10_Level.xlsx'); % this is the original xls file converted into csv

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

% pull out SPF forecasts for quarter 0 through 4
Yforecast          = SPFimport.data(:,4:8);
YbarforecastAnnual = SPFimport.data(:,10:end);
% YbarActualforecast = Ybarforecast;

Ybarforecast10     = SPF10import.data(:,3);


% if size(Ybarforecast) ~= Nbar
%     error houston
% end

%% load realized data from current FRED (not using RTDSM data, since release of quarterly vintage may lag our data updates)

switch upper(datalabel0)
    case 'TBILL'
        FRED      = fredreadcsv('TB3MS',[], 'm', 'q', 'avg');
        Yrealized = fredexpand(FRED, realizedDataDates);
    case {'CPI'}
        FRED      = fredreadcsv('CPIAUCSL',[], 'm', 'q', 'avg');
        cpiLevel  = fredexpand(FRED, realizedDataDates);
        Yrealized = ((cpiLevel(2:end) ./ cpiLevel(1:end-1)).^4 - 1) * 100;
    otherwise
        error('datalabel %s not yet supported', upper(datalabel))
end



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
for t = 1 : T
    for ii = 1 : Nz2
        if ii == 1 && Q(t) == 4 % ignore B foreast at Q4 (since spanned by the quarterly forecasts)
            Zdata2(t,ii) = NaN; % C2 left as zero
        else
            ndxMA4          = (4 - Q(t)) + ii * 4 + (-3:0);
            ndxMA4          = ndxMA4 + 2; % add 2 for nowcast and lagged realized
            C2(ii,ndxMA4,t) = 0.25; % .25 = 1/4
        end
    end
end
Z2label    = cell(Nz2,1);
for n = 1 : Nz2
    Z2label{n} = sprintf('hbar-%d', n);
end

% PART3: add 10-year forecast
Nz3        = 1;
Zdata3     = Ybarforecast10;
Z3label    = 'hbar-40';
C3         = zeros(Nz3,Ny,T);
for t = 1 : T
    if Q(t) == 1 && ~isnan(Zdata3(t))
        % ndx jump off counts quarters left in current
        ndxMA40          = (0:39) + 2;  % add 2 for nowcast and lagged realized
        C3(1,ndxMA40,t)  = 0.025; % = 1 / 40;
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



%% construct future realized values (of calendar-year fixed-event forecasts)
YBARfuture     = NaN(T,NbarAnnual);
YCURRENTfuture = NaN(T,1);



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
Czhat   = [];
Czbar   = C2; % Add C3?
save(matfile, 'Ny', 'Nz', 'T', 'dates', 'datesQ', 'doNIPA', ...
    'Yfuture', 'Ylabel', 'Nbar',  'YBARlabel', ...
    'YactualForecasts', 'Nforecasts', 'YforecastLabel', ...
    'Zlabel', 'Zdata', 'Znanny', 'Cz', 'Czbar', 'Czhat', ...
    'YBARfuture', 'YCURRENTfuture', ...
    '-v7.3')
