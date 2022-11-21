function [f_val,f_grad,f_hess] = etilt_obj_pen(l,momdiff,w,penfact)
%% Entropic tilting objective function WITH PENALIZATION, its gradient and Hessian
tmp=exp(l'*momdiff); 
f_val=tmp*w+penfact*sum(l.^2);
f_grad=tmp.*momdiff*w+2*penfact*l;
f_hess=(tmp.*w').*momdiff*momdiff'+2*penfact*eye(length(l));
 end      