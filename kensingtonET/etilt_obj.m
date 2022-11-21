function [f_val,f_grad,f_hess] = etilt_obj(l,w,momdiff)
%% Entropic tilting objective function, gradient, Hessian
tmp=exp(l'*momdiff);
f_val=tmp*w;  
f_grad=tmp.*momdiff*w;
f_hess=(tmp.*w').*momdiff*momdiff';
end