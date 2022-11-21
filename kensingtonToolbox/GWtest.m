function [GW_stat,GW_pval]= GWtest(loss1,loss2,bw,cond,h)
%% Giacomini-White (2006, ECTA) test
% input arguments:
% loss1, loss2: (T x 1) vectors of forecast losses
% bw:           scalar bandwidth for HAC estimation
% cond:         true (conditional forecast evaluation) or false (unconditional forecast evaluation)
% h:            (T x q-1) matrix of h_t "test functions", constant must not be included
%% form loss difference
lossdiff=loss1-loss2;
%% get rid of any non-numeric (NaN, Inf) values, valid data are assumed to be contiguous
valid_idx_start=find(isfinite(lossdiff),1,'first');
valid_idx_end=find(isfinite(lossdiff),1,'last');
lossdiff=lossdiff(valid_idx_start:valid_idx_end,1);
T=length(lossdiff);
%% form h matrix
onesT=ones(T,1);
if cond==true
    h=[onesT,h(valid_idx_start:valid_idx_end,1)];
else
    h=onesT;
end
q=size(h,2);
%% form Z matrix
Z=h.*lossdiff;


Zbar=mean(Z)';
Omega=HAC_est(Z,bw);
GW_stat=T*Zbar'*inv(Omega)*Zbar;
GW_pval=1-chi2cdf(GW_stat,q);


function G0 = HAC_est(y,qn)
% input: y is a T*k matrix and qn is bandwidth
% output: the Newey-West HAC covariance estimator
% formulas are from Hayashi (2000), Section 6.6
T=size(y,1);
ybar=ones(T,1)*((sum(y))/T);
dy=y-ybar;
G0=(dy'*dy)/T;
    for j=1:qn-1
        gamma=(dy(j+1:T,:)'*dy(1:T-j,:))/T;
        weight=(1-j/qn);
        G0=G0+(gamma+gamma')*weight;
    end
end



end














