%% CGM collecting real-time real GDP growth (annualized logdiff), GDP deflator (annualized logdiff) and unemployment rate (quarterly average of monthly readings) data
% also collecting calendar-year growth rates (RGDP, PGDP) and annual averages (UNRATE, CPI)
clear; close all; clc;
addpath ../matlabtoolbox/emtools/

%% INSTRUCTIONS
% - download xlsx fils via getKensingtonRawData.ipynb
% - run this script to collect the data

%% real GDP, first vintage is 65Q4, first observation is 47Q1
% Import the data
data_tmp = importdata("ROUTPUTQvQd.xlsx");
ROUTPUTQvQd = data_tmp.data;


% set observation and vintage dates
obsdates_RGDP=dateshift(datetime('1947Q1','InputFormat','yyyyQQQ'),'start','quarter',(0:size(ROUTPUTQvQd,1)-1))';
vindates_RGDP=dateshift(datetime('1965Q4','InputFormat','yyyyQQQ'),'start','quarter',(0:size(ROUTPUTQvQd,2)-1));

%% patch in 1996:Feb vintage from ALFRED (to replace missing obs in RTDSM for 1996Q1)
alfred    = readtable('GDPC1_2.xls');
vin1996Q1 = vindates_RGDP == datetime(1996,1,1);

ROUTPUTQvQd(:,vin1996Q1) = NaN; 
[~,ndx1,ndx2] = intersect(obsdates_RGDP,alfred.observation_date);
ROUTPUTQvQd(ndx1,vin1996Q1) = alfred.GDPC1_19960223(ndx2);

%% continue with real GDP collection
% calculate 400*logdiffs
annualized_QoQ_RGDP=400*diff(log(ROUTPUTQvQd),1);
obsdates_RGDP_gr=obsdates_RGDP(2:end);
RGDP.data=annualized_QoQ_RGDP;
RGDP.leveldata=ROUTPUTQvQd(2:end,:); % dropping initial level to conform with obsdata
RGDP.obsdates=obsdates_RGDP_gr;
RGDP.vindates=vindates_RGDP;

% calculate calendar-year growth rates using the 2nd readings, e.g. for 1964 -> 1965, use 1966Q2 vintage

% first calendar-year growth for which we have RT data is 1964 -> 1965, 2nd reading of 1965Q4 was available in 1966Q2
first_calY_vin = find(vindates_RGDP==datetime(1966,04,01));
first_calY_obsQ4 = find(obsdates_RGDP==datetime(1965,10,01)); % last quarter of second annual observation, goes into numerator

last_vin_Q = quarter(vindates_RGDP(1,end));
if last_vin_Q >= 2 % OK, we can use it
    last_vin_Y = year(vindates_RGDP(1,end));
else
    last_vin_Y = year(vindates_RGDP(1,end)) - 1;
end

last_calY_vin = find(vindates_RGDP==datetime(strcat(num2str(last_vin_Y),'Q2'),'InputFormat','yyyyQQQ'));

numcalY = (last_calY_vin - first_calY_vin)/4+1;
RGDP_calY_gr=NaN(numcalY,1);
t_idx = 1;
obs_idx = first_calY_obsQ4;
for y_idx = first_calY_vin : 4 : last_calY_vin
    
    RGDP_calY_gr(t_idx,1) = 100 * ( sum(ROUTPUTQvQd(obs_idx-3:obs_idx,y_idx)) / sum(ROUTPUTQvQd(obs_idx-7:obs_idx-4,y_idx)) - 1);

    t_idx = t_idx + 1;
    obs_idx = obs_idx + 4;

end

RGDP_calY_gr_dates = 1965 + (0:numcalY-1)';

RGDP.calY_gr = RGDP_calY_gr;
RGDP.calY_gr_dates = RGDP_calY_gr_dates;

%% GDP deflator, first vintage is 65Q4, first observation is 47Q1
% Import the data
data_tmp = importdata("PQvQd.xlsx");
PQvQd = data_tmp.data;


% set observation and vintage dates
obsdates_PGDP=dateshift(datetime('1947Q1','InputFormat','yyyyQQQ'),'start','quarter',(0:size(PQvQd,1)-1))';
vindates_PGDP=dateshift(datetime('1965Q4','InputFormat','yyyyQQQ'),'start','quarter',(0:size(PQvQd,2)-1));

%% patch in 1996:Feb vintage from ALFRED (to replace missing obs in RTDSM for 1996Q1)
alfred    = readtable('GDPCTPI_2.xls');
vin1996Q1 = vindates_PGDP == datetime(1996,1,1);

PQvQd(:,vin1996Q1)    = NaN; 
[~,ndx1,ndx2]         = intersect(obsdates_PGDP,alfred.observation_date);
PQvQd(ndx1,vin1996Q1) = alfred.GDPCTPI_19960223(ndx2);

%% continue PGDP
% calculate 400*logdiffs
annualized_QoQ_PGDP=400*diff(log(PQvQd),1);
obsdates_PGDP_gr=obsdates_PGDP(2:end);
PGDP.data=annualized_QoQ_PGDP;
PGDP.leveldata=PQvQd(2:end,:); % dropping initial level to conform with obsdata
PGDP.obsdates=obsdates_PGDP_gr;
PGDP.vindates=vindates_PGDP;

% calculate calendar-year growth rates using the 2nd readings, e.g. for 1964 -> 1965, use 1966Q2 vintage

% first calendar-year growth for which we have RT data is 1964 -> 1965, 2nd reading of 1965Q4 was available in 1966Q2
first_calY_vin = find(vindates_PGDP==datetime(1966,04,01));
first_calY_obsQ4 = find(obsdates_PGDP==datetime(1965,10,01)); % last quarter of second annual observation, goes into numerator

last_vin_Q = quarter(vindates_PGDP(1,end));
if last_vin_Q >= 2 % OK, we can use it
    last_vin_Y = year(vindates_PGDP(1,end));
else
    last_vin_Y = year(vindates_PGDP(1,end)) - 1;
end

last_calY_vin = find(vindates_PGDP==datetime(strcat(num2str(last_vin_Y),'Q2'),'InputFormat','yyyyQQQ'));

numcalY = (last_calY_vin - first_calY_vin)/4+1;
PGDP_calY_gr=NaN(numcalY,1);
t_idx = 1;
obs_idx = first_calY_obsQ4;
for y_idx = first_calY_vin : 4 : last_calY_vin
    
    PGDP_calY_gr(t_idx,1) = 100 * ( sum(PQvQd(obs_idx-3:obs_idx,y_idx)) / sum(PQvQd(obs_idx-7:obs_idx-4,y_idx)) - 1);

    t_idx = t_idx + 1;
    obs_idx = obs_idx + 4;

end

PGDP_calY_gr_dates = 1965 + (0:numcalY-1)';

PGDP.calY_gr = PGDP_calY_gr;
PGDP.calY_gr_dates = PGDP_calY_gr_dates;

%% Unemployment rate, first vintage is 65Q4, first observation is 47M1 in the table, we keep it that way

% Import the data
data_tmp = importdata("rucQvMd.xlsx");
UNRATEQvMd = data_tmp.data;
UNRATEQvMd = [NaN(12,size(UNRATEQvMd,2));UNRATEQvMd];
% important: importdata chops off the first 12 observations, because they are all NaN - so we put them back


% need to convert monthly levels to quarterly levels by averaging
% we keep the monthly values in quarters when only 1 monthly observation was available
numM=size(UNRATEQvMd,1);
numQ=ceil(numM/3);
numvin=size(UNRATEQvMd,2);
obsdates_UNRATE=dateshift(datetime('1947Q1','InputFormat','yyyyQQQ'),'start','quarter',(0:numQ-1))';
vindates_UNRATE=dateshift(datetime('1965Q4','InputFormat','yyyyQQQ'),'start','quarter',(0:numvin-1));

UNRATEQvQd=NaN(numQ,numvin);
for v=1:length(vindates_UNRATE)
    startM=1;

    for d=1:length(obsdates_UNRATE)
        if d*3<numM
            endM=d*3;
        else
            endM=numM;
        end
        if (endM-startM)==2 % only use full quarters
            UNRATEQvQd(d,v)=mean(UNRATEQvMd(startM:endM,v));
        end

        startM=startM+3;
    end
end

UNRATE.data=UNRATEQvQd;
UNRATE.obsdates=obsdates_UNRATE;
UNRATE.vindates=vindates_UNRATE;

% calculate calendar-year annual AVERAGE levels
UNRATE_timetable = array2timetable(UNRATE.data(:,end),"RowTimes",UNRATE.obsdates);
UNRATE_calY_avg_timetable = retime(UNRATE_timetable,'yearly','mean');
UNRATE_calY_avg = UNRATE_calY_avg_timetable.Var1;
UNRATE_calY_avg_dates = year(UNRATE_calY_avg_timetable.Time);

UNRATE.calY_avg = UNRATE_calY_avg;
UNRATE.calY_dates = UNRATE_calY_avg_dates;
%% CPI, first vintage is 65Q4, first observation is 47M1 in the table, we keep it that way
% but actual vintages start in 94Q3

data_tmp = importdata("cpiQvMd.xlsx");
CPIQvMd = data_tmp.data;
CPIQvMd = [NaN(size(CPIQvMd,1),115),CPIQvMd];
% important: importdata chops off first 115 columns, which are all NaNs, as there are no vintages from before 1994Q3
% so we put it back 

% need to convert monthly levels to quarterly levels by averaging
% we keep the monthly values in quarters when only 1 monthly observation was available
numM=size(CPIQvMd,1);
numQ=ceil(numM/3);
numvin=size(CPIQvMd,2);
obsdates_CPI=dateshift(datetime('1947Q1','InputFormat','yyyyQQQ'),'start','quarter',(0:numQ-1))';
vindates_CPI=dateshift(datetime('1965Q4','InputFormat','yyyyQQQ'),'start','quarter',(0:numvin-1));

CPIQvQd=NaN(numQ,numvin);
for v=1:length(vindates_CPI)
    startM=1;

    for d=1:length(obsdates_CPI)
        if d*3<numM
            endM=d*3;
        else
            endM=numM;
        end
        if (endM-startM)==2 % only use full quarters
            CPIQvQd(d,v)=mean(CPIQvMd(startM:endM,v));
        end

        startM=startM+3;
    end
end

% calculate annualized QoQ growth rate, ((X_{t}/X_{t-1})^4-1)*100

annualized_QoQ_CPI=((CPIQvQd(2:end,:)./(CPIQvQd(1:end-1,:))).^4-1)*100;

CPI.data=annualized_QoQ_CPI;
CPI.leveldata=CPIQvQd(2:end,:);
CPI.obsdates=obsdates_CPI(2:end);
CPI.vindates=vindates_CPI;

% TODO: compute calendar-year Q4/Q4 CPI inflation rates

%% save RT data --- separate files for each variable

DATALABELS = {'RGDP','PGDP','UNRATE','CPI'}; % ,'TBILL'};

for d = 1 : length(DATALABELS)
    datalabel = DATALABELS{d};
    save(sprintf('kensington%sdataRT.mat', datalabel), '-struct', datalabel, '-v7.3');
end

