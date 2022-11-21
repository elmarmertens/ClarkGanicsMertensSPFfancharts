function [GR_stat,GR_location,GR_critval,GR_decision,F_oos_vec] = GRtest(lossdiff,m,bw,sided,alp)
%% Giacomini-Rossi (2010, JAE) Fluctuation test
% inputs:   lossdiff    : (P x 1) times vector of forecast loss differences
%           m           : scalar window size to calculate local average loss differences
%           bw          : bandwidth for Newey-West (1987, ECTA) HAC estimator
%           sided       : 1 or 2, corresponding to one-sided or two-sided test
%           alp         : significance level, must be either 0.05 or 0.10 (these are tabulated)
% outputs:  GR_stat     : Giacomini-Rossi test statistic (either max(abs(F_oos_vec) or max(F_oos_vec))
%           GR_location : location of GR_stat, 
%                         if window size is odd, then center of window, if window size is even, then half window size plus 1
%           GR_critval  : critical value of Giacomini-Rossi test,
%                         linearly interpolating over ratio of window size and OOS sample size
%           GR_decision : 1 if we reject H0, 0 if we fail to reject H0
%           F_oos_vec   : vector of rescaled (standardized by HAC st. dev.) local out-of-sample loss differences

% check input arguments
if any(sided ~= 1 & sided ~= 2)
    error('sided must be 1 or 2!');
end
if any(alp ~= 0.05 & alp ~=0.10)
    error('alp must be 0.05 or 0.10!');
end


P = length(lossdiff);
F_oos_vec = NaN(P,1);
halfwindow = floor(m/2);
if mod(m,2) == 0
    m_even = true;
    t_end = P - halfwindow + 1;
else
    m_even = false;
    t_end = P - halfwindow;
end

HACsig2hat = HAC_est(lossdiff,bw); % GR eq. 2

for t = halfwindow+1 : t_end

    start_idx = t - halfwindow;
    if m_even == true
        end_idx = t + halfwindow - 1;
    else
        end_idx = t + halfwindow;
    end

    F_oos_vec(t,1) = 1/sqrt(m) * sum(lossdiff(start_idx:end_idx,1)) / sqrt(HACsig2hat); % GR eq. 1

end

if sided == 1

    [GR_stat,GR_location] = max(F_oos_vec);
    
    if alp == 0.05
        GR_tabcol = 4;
    elseif alp == 0.10
        GR_tabcol = 5;
    end

elseif sided == 2

    [GR_stat,GR_location] = max(abs(F_oos_vec));

    if alp == 0.05
        GR_tabcol = 2;
    elseif alp == 0.10
        GR_tabcol = 3;
    end

end

% Giacomini-Rossi's Table I below

GR_tab = [0.1, 3.393, 3.170, 3.176, 2.928;
          0.2, 3.179, 2.948, 2.938, 2.676;
          0.3, 3.012, 2.766, 2.770, 2.482;
          0.4, 2.890, 2.626, 2.624, 2.334;
          0.5, 2.779, 2.500, 2.475, 2.168;
          0.6, 2.634, 2.356, 2.352, 2.030;
          0.7, 2.560, 2.252, 2.248, 1.904;
          0.8, 2.433, 2.130, 2.080, 1.740;
          0.9, 2.248, 1.950, 1.975, 1.600];

GR_critval = interp1(GR_tab(:,1),GR_tab(:,GR_tabcol),m/P);

if GR_stat > GR_critval
    GR_decision = 1;
else
    GR_decision = 0;
end

function G0 = HAC_est(y,qn)
% input: y is a T*k matrix and qn is bandwidth
% output: the Newey-West HAC covariance estimator
% formulas are from Hayashi (2000), Section 6.6
T = size(y,1);
ybar = ones(T,1) * ((sum(y))/T);
dy = y - ybar;
G0 = (dy'*dy) / T;
    for j = 1 : qn-1
        gamma = (dy(j+1:T,:)'*dy(1:T-j,:)) / T;
        weight= 1 - j / qn;
        G0 = G0 + (gamma+gamma') * weight;
    end
end

end