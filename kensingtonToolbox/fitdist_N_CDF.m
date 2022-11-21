function [paramvec_N,resnorm_N,residual_N] = fitdist_N_CDF(cdfvals,cdfpoints,lsqnl_optimoptions)
%% this function fits a normal distribution to histograms by minimizing the SSR between the histogram's CDF and the normal CDF
% cdfvals:    N x 1 vector of cumulative probabilities from histogram
% cdfpoints:   N x 1 vector of histogram bounds (where CDF is evaluated)

cdf_diff = @(paramvec)normcdf(cdfpoints,paramvec(1),paramvec(2))-cdfvals;

% get reasonable starting values
cdf_midpoints = (cdfpoints(2:end,1)+cdfpoints(1:end-1,1))/2;
pdf_vals = cdfvals(2:end,1)-cdfvals(1:end-1,1);
mu0 = cdf_midpoints'*pdf_vals;
sig0 = sqrt(((cdf_midpoints-mu0).^2)'*pdf_vals);

[paramvec_N,resnorm_N,residual_N,~,~] = lsqnonlin(cdf_diff,[mu0;sig0],[-Inf;eps],[Inf,Inf],lsqnl_optimoptions);

end