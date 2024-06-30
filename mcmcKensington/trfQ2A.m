function [YBARdensityDraws,YCURRENTdensityDraws] = trfQ2A(doNIPA,thisT,datesQ,RT_obsdate_idx,YdensityDraws,RT_vec,Ndraws,Nbar)
%% CGM: here we transform the QoQ (400*logdiff) and Q level forecasts into the SPF concepts (growth rate of calendar-year avg. and annual level)

%% construct predictive densities for current year

% first identify in which quarter we are now
Q_act=datesQ(thisT);

% how many lagged growth rates/levels we HAVE - these are DATA, not forecasts
% in Q1, last year's Q1, Q2, Q3, Q4 were known
% in Q2, last year's Q1, Q2, Q3, Q4 and current year's Q1 were known
% in Q3, last year's Q1, Q2, Q3, Q4 and current year's Q1, Q2 were known
% in Q4, last year's Q1, Q2, Q3, Q4 and current year's Q1, Q2, Q3 were known

% extract last year's Q1-Q4 growth rates/levels
% we actually only use data starting with last year's Q2 when computing
% growth rates but for notational reasons, it's easier to extract Q1 as well
% obsdate is last available observation, row index of RTdata
last_Q4=RT_obsdate_idx-Q_act+1;

% pad values to RTvec if vector (could be matrix if draws for missing values were added)

if doNIPA

    last_year_data=RT_vec(last_Q4-3:last_Q4,:); % only need last year's data when working with growth rates

    C=1;
    for k=2:4
        C=C+exp(sum(last_year_data(2:k,:),1)./400);
    end

end

% extract current year's growth rate/level DATA, if any is available
this_Q1           = last_Q4+1;
this_year_data    = RT_vec(this_Q1:RT_obsdate_idx,:);
if isvector(this_year_data)
    this_year_data = repmat(this_year_data, 1, Ndraws);
end

% extract current year's forecasts
this_year_fc      = YdensityDraws(1:5-Q_act,:);
this_year_data_fc = [this_year_data;this_year_fc];

if doNIPA

    B=ones(1,Ndraws);
    for k=2:4
        B = B + exp(sum(this_year_data_fc(2:k,:),1) ./ 400);
    end

    A = exp((sum(last_year_data(2:4,:),1) + this_year_data_fc(1,:)) ./ 400);

    YCURRENTdensityDraws = (A .* B ./ C - 1) * 100;

else

    YCURRENTdensityDraws = mean(this_year_data_fc,1);
end

%% now go on to next year, and 2 and 3 years... out
future_start=5-Q_act+1; % index in forecast matrix
future_end=future_start+3;  % index in forecast matrix
YBARdensityDraws=NaN(Nbar,Ndraws);

if doNIPA
    A_term1 = this_year_data_fc;
end

ONEs = ones(1,Ndraws);

for hz=1:Nbar

    future_fc = YdensityDraws(future_start:future_end,:);

    if doNIPA

        C = B; % old B becomes C
        B = ONEs;

        for k=2:4
            B = B + exp(sum(future_fc(2:k,:),1)./400);
        end

        A       = exp((sum(A_term1(2:4,:),1) + future_fc(1,:)) ./ 400);
        A_term1 = future_fc; % to be used at next iteration

        YBARdensityDraws(hz,1:Ndraws)=(A.*B./C - 1) * 100;

    else

        YBARdensityDraws(hz,1:Ndraws)=mean(future_fc,1);

    end

    future_start= future_end+1; % to be used at next iteration
    future_end  = future_start+3; % to be used at next iteration


end

end