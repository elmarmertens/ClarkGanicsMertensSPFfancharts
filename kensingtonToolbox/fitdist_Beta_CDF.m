function [paramvec_B,resnorm_B,residual_B] = fitdist_Beta_CDF(cdfvals,cdfpoints,fix_l,fix_r,opt_minmax,lsqnl_optimoptions)
%% this function fits a generalized beta distribution to histograms by minimizing the SSR between the histogram's CDF and the generalized beta CDF
% cdfvals:      N x 1 vector of cumulative probabilities from histogram
% cdfpoints:    N x 1 vector of histogram bounds (where CDF is evaluated)
% fix_l, fix_r: booleans, 1 if left or right endpoint of support is fixed
% opt_minmax:   2 x 1 vector of lb and ub of support

if fix_l == true && fix_r == false % lb of *support* fixed but ub of *support* not fixed

    cdf_diff = @(paramvec)gen_beta_cdf([opt_minmax(1,1);paramvec],cdfpoints)-cdfvals;
    opt_lb = [opt_minmax(1,1)+eps;1+eps;1+eps];
    opt_ub = [opt_minmax(2,1);Inf;Inf];
    paramvec0 = [mean(opt_minmax);1.5;1.5];
    param_case = 1;

elseif fix_l == false && fix_r == true % lb of *support* is not fixed but ub of *support* is fixed

    cdf_diff = @(paramvec)gen_beta_cdf([paramvec(1,1);opt_minmax(2,1);paramvec(2:3,1)],cdfpoints)-cdfvals;
    opt_lb = [opt_minmax(1,1);1+eps;1+eps];
    opt_ub = [opt_minmax(2,1)-eps;Inf;Inf];
    paramvec0 = [mean(opt_minmax);1.5;1.5];
    param_case = 2;

elseif fix_l == true && fix_r == true % both lb and ub of *support* are fixed

    cdf_diff = @(paramvec)gen_beta_cdf([opt_minmax(1:2,1);paramvec(1:2,1)],cdfpoints)-cdfvals;
    opt_lb = [1+eps;1+eps];
    opt_ub = [Inf;Inf];
    paramvec0 = [1.5;1.5];
    param_case = 3;

else % neither the lb nor the ub of the *support* are fixed

    cdf_diff = @(paramvec)gen_beta_cdf(paramvec,cdfpoints)-cdfvals;
    opt_lb = [opt_minmax(1,1);opt_minmax(1,1)+eps;1+eps;1+eps];
    opt_ub = [opt_minmax(2,1)-eps;opt_minmax(2,1);Inf;Inf];
    paramvec0 = [0.8*opt_minmax(1,1) + 0.2*opt_minmax(2,1);0.2*opt_minmax(1,1) + 0.8*opt_minmax(2,1);1.5;1.5];
    param_case = 4;

end

[paramvec_B,resnorm_B,residual_B,~,~] = lsqnonlin(cdf_diff,paramvec0,opt_lb,opt_ub,lsqnl_optimoptions);

switch param_case
    case 1
        paramvec_B = [opt_minmax(1,1);paramvec_B];
    case 2
        paramvec_B = [paramvec_B(1,1);opt_minmax(2,1);paramvec_B(2:3,1)];
    case 3
        paramvec_B = [opt_minmax(1:2,1);paramvec_B(1:2,1)];
end

end