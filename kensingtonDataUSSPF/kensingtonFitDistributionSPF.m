%% CGM: fitting normal and generalized beta distributions to aggregate SPF histograms
clear; close all; clc;
%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emeconometrics/
addpath ../kensingtonToolbox/

%% following Engelberg et al. (2009, JBES) we treat bins as closed from the left and open from the right
% bin right endpoints (and infinity as right endpoint of rightmost bin), these are valid as of the 2022Q2 SPF round
% if SPF bins change after that date, this needs to be extended

% % PGDP
% PRPGDP1_bins = -3 : 1 : 10; % 1968Q4-1973Q1
% PRPGDP2_bins = -1 : 1 : 12; % 1973Q2-1974Q3
% PRPGDP3_bins =  3 : 1 : 16; % 1974Q4-1981Q2
% PRPGDP4_bins =  4 : 2 : 12; % 1981Q3-1985Q1  next year got included
% PRPGDP5_bins =  2 : 2 : 10; % 1985Q2-1991Q4
% PRPGDP6_bins =  0 : 1 : 8;  % 1992Q1-2013Q4
% PRPGDP7_bins =  0 : 0.5 : 4; % 2014Q1-present
% 
% % RGDP
% PRGDP1_bins = PRPGDP1_bins; % 1968Q4-1973Q1
% PRGDP2_bins = PRPGDP2_bins; % 1973Q2-1974Q3
% PRGDP3_bins = PRPGDP3_bins; % 1974Q4-1981Q2
% PRGDP4_bins = -2 : 2 : 6; % 1981Q3-1991Q4  next year got included
% PRGDP5_bins = -2 : 1 : 6; % 1992Q1-2009Q1
% PRGDP6_bins = -3 : 1 : 6; % 2009Q2-2020Q1  years 3 and 4 got included
% PRGDP7_bins = [-12,-6,-3,0,1.5,2.5,4,7,10,16]; % 2020Q2-present
% 
% % UNEMP
% PRUNEMP1_bins = [6,7,7.5,8,8.5,9,9.5,10,11]; % 2009Q2-2013Q4     next, years 3 and 4 included
% PRUNEMP2_bins = [4,5,5.5,6,6.5,7,7.5,8,9]; % 2014Q1-2020Q1
% PRUNEMP3_bins = [3,4,5,6,7,8,10,12,15]; % 2020Q2-present

%% below we use the same points on the CDF as for ET - option 2 (Aug 15)

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

% set optimization options
lsqnl_optimoptions = optimoptions('lsqnonlin','Display','off');
tol_N = 0.01; % tolerance for maximum absolute deviation of fitted normal CDF from cumulative histograms
tol_B = 0.01; % tolerance for maximum absolute deviation of fitted beta CDF from cumulative histograms
tol_SSR_N = 0.01; % tolerance for sum of squared residuals computed as difference of fitted normal CDF from cumulative histograms
tol_SSR_B = 0.01; % tolerance for sum of squared residuals computed as difference of fitted beta CDF from cumulative histograms

% to follow Engelberg et al. (2009, JBES), we need to set the bounds of the generalized beta distribution
RGDP_tmp = importdata('ROUTPUTQvQd.xlsx');
% find 1991Q1 vintage
vin91Q1_idx = find(strcmpi(RGDP_tmp.textdata(1,2:end),'ROUTPUT91Q1'));
RGDP_Q_realiz = RGDP_tmp.data(:,vin91Q1_idx);
RGDP_Q_dates = datenum(cell2mat(RGDP_tmp.textdata(2:end,1)),'yyyy:QQ');
RGDP_timetab = timetable(datetime(RGDP_Q_dates,'ConvertFrom','datenum'),RGDP_Q_realiz);
RGDP_A_timetab = retime(RGDP_timetab,'yearly','mean');
RGDP_A_tmp = RGDP_A_timetab.RGDP_Q_realiz;
RGDP_A_growth = (RGDP_A_tmp(2:end,1)./RGDP_A_tmp(1:end-1,1)-1)*100;
RGDP_A_growth_minmax = [min(RGDP_A_growth);max(RGDP_A_growth)];
RGDP_A_growth_minmax = [-5; 8];
% 
PGDP_tmp = importdata('PQvQd.xlsx');
% find 1991Q1 vintage
vin91Q1_idx = find(strcmpi(PGDP_tmp.textdata(1,2:end),'P91Q1'));
PGDP_Q_realiz = PGDP_tmp.data(:,vin91Q1_idx);
PGDP_Q_dates = datenum(cell2mat(PGDP_tmp.textdata(2:end,1)),'yyyy:QQ');
PGDP_timetab = timetable(datetime(PGDP_Q_dates,'ConvertFrom','datenum'),PGDP_Q_realiz);
PGDP_A_timetab = retime(PGDP_timetab,'yearly','mean');
PGDP_A_tmp = PGDP_A_timetab.PGDP_Q_realiz;
PGDP_A_growth = (PGDP_A_tmp(2:end,1)./PGDP_A_tmp(1:end-1,1)-1)*100;
PGDP_A_growth_minmax = [min(PGDP_A_growth);max(PGDP_A_growth)];
PGDP_A_growth_minmax = [-1; 10];
% 
UNRATE_tmp = importdata('UNRATE.csv');
UNRATE_M_realiz = UNRATE_tmp.data(:,end);
UNRATE_M_dates = datenum(cell2mat(UNRATE_tmp.textdata(2:end,1)),'yyyy-mm-dd');
UNRATE_timetab = timetable(datetime(UNRATE_M_dates,'ConvertFrom','datenum'),UNRATE_M_realiz);
UNRATE_A_timetab = retime(UNRATE_timetab,'yearly','mean');
UNRATE_A_tmp = UNRATE_A_timetab.UNRATE_M_realiz;
UNRATE_A_tmp_90idx = find(UNRATE_A_timetab.Time == datetime(1990,01,01));
UNRATE_A_minmax = [min(UNRATE_A_tmp(1:UNRATE_A_tmp_90idx));max(UNRATE_A_tmp(1:UNRATE_A_tmp_90idx))];
UNRATE_A_minmax = [2; 16];


%RGDP_A_growth_minmax = [-12.517;19.74921]; % ALFRED
%PGDP_A_growth_minmax = [-11.69375;12.94057]; % ALFRED
%UNRATE_A_minmax = [1.2;24.9]; % Table 1 in Stanley Lebergott (1957) Annual Estimates of Unemployment in the United States,
% 1900-1954, In: The Measurement and Behavior of Unemployment, National Bureau of Economic Research. URL: https://www.nber.org/system/files/chapters/c2644/c2644.pdf

data_historical_minmax = [RGDP_A_growth_minmax,PGDP_A_growth_minmax,UNRATE_A_minmax];

startdates = {'PRGDP', 'PRPGDP', 'PRUNEMP'; datenum('1992Q1','yyyyQQ'), datenum('1992Q1','yyyyQQ'), datenum('2009Q2','yyyyQQ')};

%% loop over variables
for d=1:length(DATALABELS)

    datalabel = DATALABELS{d};
    % setup log file
        diaryfilename=fullfile(localtemp, sprintf('kensingtonFitDistributionSPF-%s.log', datalabel));
        if exist(diaryfilename, 'file')
            delete(diaryfilename);
        end
        diary(diaryfilename);

    data_act = SPFimport.data.(datalabel);

    % check dates
    Y = data_act(:,1);
    if ~isequal(Y, year(dates))
        error('inconsistent date vector in the SPF file (years)');
    end
    Q = data_act(:,2);
    if ~isequal(Q, quarter(dates))
        error('inconsistent date vector in the SPF file (quarters)');
    end

    data_probs_SPF = data_act(:,3:end);
    data_probs_CDF = NaN(size(data_probs_SPF));
    fitted_N_params = NaN(size(data_probs_SPF,1),4,2);
    fitted_B_params = NaN(size(data_probs_SPF,1),4,4);
    resnorm_N = NaN(size(data_probs_SPF,1),4);
    resnorm_B = NaN(size(data_probs_SPF,1),4);
    allgrids=NaN(size(dates,1),60); % first column: number of horizons,
                                    % after that repeating the bin edges (as many times as the number of horizons - adjust this if new horizons are added after 2022Q2)

    
    this_tstart = find(dates==startdates{2,strcmpi(startdates(1,:),datalabel)});

    for t=this_tstart:size(dates,1)

        if strcmp(datalabel,'PRGDP')==1

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

        elseif strcmp(datalabel,'PRPGDP')==1

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

        elseif strcmp(datalabel,'PRUNEMP')==1

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
            thisbinlength = length(thisbin) + 1;

            
            for hz = 1 : numhz
                % flip probabilities, as in SPF file the bins are arranged from highest to lowest, then take cumulative sum
                data_probs_PDF_tmp = fliplr(data_probs_SPF(t,(hz - 1) *thisbinlength + 1 : hz * thisbinlength));
                data_probs_CDF_tmp = cumsum(data_probs_PDF_tmp)/100;
                data_probs_CDF(t,(hz - 1) *thisbinlength + 1 : hz * thisbinlength) = data_probs_CDF_tmp;

                % fit normal distribution through CDF (Giordani-Söderlind, 2003, EER)
                [fitted_N_params(t,hz,1:2),resnorm_N(t,hz),residual_N] = fitdist_N_CDF(data_probs_CDF_tmp(1:end-1)',thisbin',lsqnl_optimoptions);
                
                %maxabsres_N = max(abs(residual_N));
                %if maxabsres_N > tol_N
                %    warning("NORMAL %s max. abs. residual is %1.2f, larger than tol_N %1.2f in %s at horizon %d",DATALABELS{d},maxabsres_N,tol_N,datestr(dates(t),'yyyyQQ'),hz)
                %end

                if resnorm_N(t,hz) > tol_SSR_N
                    warning('NORMAL %s SSR is %1.2f, larger than tol_SSR_N %1.2f in %s at horizon %d',DATALABELS{d},resnorm_N(t,hz),tol_SSR_N,datestr(dates(t),'yyyyQQ'),hz);
                end

                % fit generalized beta distribution through CDF (Engelberg et al., 2009, JBES)

                opt_minmax = data_historical_minmax(1:2,d);
                fix_l = false;
                fix_r = false;

                % find first and last bins with nonzero probability
                first_nonzero = find(data_probs_PDF_tmp>0,1,'first');
                last_nonzero = find(data_probs_PDF_tmp>0,1,'last');

                if first_nonzero > 1
                    opt_minmax(1,1) = thisbin(first_nonzero-1); % this is the lb of the *support*
                    fix_l = true; % l parameter is fixed in optimization
                end

                if last_nonzero < thisbinlength
                    opt_minmax(2,1) = thisbin(last_nonzero); % this is the ub of the *support*
                    fix_r = true; % r parameter is fixed in optimization
                end

                [fitted_B_params(t,hz,1:4),resnorm_B(t,hz),residual_B] = fitdist_Beta_CDF(data_probs_CDF_tmp(1:end-1)',thisbin',fix_l,fix_r,opt_minmax,lsqnl_optimoptions);

                %maxabsres_B = max(abs(residual_B));
                %if maxabsres_B > tol_B
                %    warning('BETA %s max. abs. residual is %1.2f, larger than tol_B %1.2f in %s at horizon %d',DATALABELS{d},maxabsres_B,tol_B,datestr(dates(t),'yyyyQQ'),hz);
                %end

                if resnorm_B(t,hz) > tol_SSR_B
                    warning('BETA %s SSR is %1.2f, larger than tol_SSR_B %1.2f in %s at horizon %d',DATALABELS{d},resnorm_B(t,hz),tol_SSR_B,datestr(dates(t),'yyyyQQ'),hz);
                end


            end

        end

    end

    tmp = sum(any(~isnan(allgrids),1)); % keep only up to the last non-NaN column 
    allgrids=allgrids(:,1:tmp);

    if strcmp(datalabel,'PRPGDP')==1
        PGDP.bins=allgrids;
        PGDP.cdf=data_probs_CDF;
        PGDP.dates=dates;
        PGDP.params_N = fitted_N_params;
        PGDP.params_B = fitted_B_params;
        PGDP.resnorm_N = resnorm_N;
        PGDP.resnorm_B = resnorm_B;
        savename_string = 'PGDP';
    elseif strcmp(datalabel,'PRGDP')==1
        RGDP.bins=allgrids;
        RGDP.cdf=data_probs_CDF;
        RGDP.dates=dates;
        RGDP.params_N = fitted_N_params;
        RGDP.params_B = fitted_B_params;
        RGDP.resnorm_N = resnorm_N;
        RGDP.resnorm_B = resnorm_B;
        savename_string = 'RGDP';
    elseif strcmp(datalabel,'PRUNEMP')==1
        UNRATE.bins=allgrids;
        UNRATE.cdf=data_probs_CDF;
        UNRATE.dates=dates;
        UNRATE.params_N = fitted_N_params;
        UNRATE.params_B = fitted_B_params;
        UNRATE.resnorm_N = resnorm_N;
        UNRATE.resnorm_B = resnorm_B;
        savename_string = 'UNRATE';
    end

diary off

%% save results separately
save(sprintf('kensington%sFitDistSPF.mat', savename_string), '-struct', savename_string, '-v7.3');

end

%% save joint file
save('kensingtonFitDistSPF.mat', 'PGDP', 'RGDP', 'UNRATE', '-v7.3');
