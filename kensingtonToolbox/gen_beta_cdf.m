function [beta_cdf_val] = gen_beta_cdf(paramvec,t_arg_vec)
%% evaluating the generalized beta distribution's CDF parametrized by paramvec at each point in t_arg_vec
% Engelberg et al. (2009, JBES), p. 38, eq. 1

l = paramvec(1); % lower bound of support
r = paramvec(2); % upper bound of support
a = paramvec(3); % shape alpha
b = paramvec(4); % shape beta

num_points = length(t_arg_vec);
beta_cdf_val = NaN(num_points,1);

if l >= r
    %warning('Evaluating generalized beta CDF: lower bound of support must be less than upper bound!');
    return
end

for idx = 1 : num_points

    if t_arg_vec(idx,1) <= l
    
        beta_cdf_val(idx,1) = 0;
    
    elseif t_arg_vec(idx,1) > r
    
        beta_cdf_val(idx,1) = 1;
    
    else
    
        fun_tmp = @(x)((x-l).^(a-1) .* (r-x).^(b-1)) ./ (r-l).^(a+b-1);
        beta_integral = integral(fun_tmp,l,t_arg_vec(idx,1));
        beta_cdf_val(idx,1) = 1/beta(a,b) * beta_integral;
    
    end

end

end