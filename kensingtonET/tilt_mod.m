function output = tilt_mod(g_data,w,gamma_init,gbar,prob_tol,point_tol,numprobs,opt_options)
%% Function for entropic, nonparametric tilting
%% INPUTS  
% g_data     : data matrix AFTER transformed by g, size is (m x D) where m is the number of variables, D is the number of draws
%                Note: moment condition is E[g(data)]=gbar
% w          : original weights, size is (D x 1), usually ones(D,1)/D
% gbar       : desired moments, size is (m x 1)
% tol        : tolerance for achieving the desired moments, the program displays a warning if the maximum deviation is larger than tol
% opt_options: optimization options for fminunc
%% OUTPUTS
% output is a structure containing:
% output.gamma      : optimal tilting parameter, size is (m x 1), 'gammastar' in Kruger, Clark and Ravazzolo (2017,JBES)
% output.w          : new weights, size is (1 x D)
% output.klic       : KLIC difference between original and tilted densities
% output.shrink     : amount of L2 shrinkage on the tilting parameter vector in case a penalization was necessary to obtain a solution
% output.convergence: EXITFLAG of minimization, if convergence >= 1 then nu should be zero
% output.w_orig     : original weight vector, usually ones(D,1)/D
% output.data       : the g(data) that was used

%% Get dimensions
[m,D] = size(g_data);  % dim of random vector x sample size
momdiff=g_data-repmat(gbar,1,D); % g(data)-gbar
obj_fun=@(l)etilt_obj(l,w,momdiff);
numshrink=20;
  
%gs_problem=createOptimProblem('fmincon','x0',zeros(m,1),'objective',obj_fun,'options',opt_options);
%gs = GlobalSearch('Display','iter');
%tpoints = CustomStartPointSet(randn(50,m));
%% Minimize objective function
try 
    [gammastar,fval,exitflag] = fminunc(obj_fun,gamma_init,opt_options);%run(gs,gs_problem);
    pen=fval; % base for penalization below
    %gammastar
catch
    [gammastar,fval,exitflag] = fminunc(obj_fun,0.001+gamma_init,opt_options);
    pen=fval; % base for penalization below
    sprintf("CATCH ACTIVE")
end
%% Compute new weights
%gammastar
%pause;
ws = w'.*exp(gammastar'*g_data);
ws = ws/sum(ws);
finalshrink = 0;

% check if moment condition satisfied
empirical_mom = g_data*ws';
absmomdiff=abs(empirical_mom-gbar);
maxabsmomdiff_probs=max(absmomdiff(1:numprobs));
maxabsmomdiff_points=max(absmomdiff(numprobs+1:end));

%% Intervene if some weights are NaN or optimization failed
if any(isnan(ws)) || (exitflag < 1) || any(maxabsmomdiff_probs > prob_tol) || any(maxabsmomdiff_points > point_tol)

        fprintf("Original problem failed - some shrinkage was used to regularize the problem. \n");
        fprintf("Number of NaNs: %d. Exitflag: %d. Maximum absolute moment difference in PROBS: %1.2f in POINTS: %1.2f \n",sum(isnan(ws)),exitflag,maxabsmomdiff_probs,maxabsmomdiff_points);
        shrink = linspace(1e-10,10,numshrink);
        gammastar_mat=NaN(m,numshrink+1);
        gammastar_mat(:,1)=gammastar;
        maxabsmomdiff_vec=NaN(1,numshrink+1);
        if isempty(maxabsmomdiff_points)
            maxabsmomdiff_vec(1,1) = maxabsmomdiff_probs;
        elseif isempty(maxabsmomdiff_probs)
            maxabsmomdiff_vec(1,1) = maxabsmomdiff_points;
        else
            maxabsmomdiff_vec(1,1) = maxabsmomdiff_probs + maxabsmomdiff_points;
        end
        exitflag_vec=NaN(1,numshrink+1);
        exitflag_vec(1,1)=exitflag;

        ok = 0;
        ct = 1;

        % Loop over various shrinkage factors (do until legitimate weights are obtained)
        while ok == 0

                % Minimize objective function
                shrinkct=shrink(ct);
                penfact=pen*shrinkct;
                obj_fun_pen=@(l)etilt_obj_pen(l,momdiff,w,penfact);
                [gammastar,~,exitflag] = fminunc(obj_fun_pen,gamma_init,opt_options);
                exitflag_vec(1,ct+1)=exitflag;
                % Compute new weights
                ws = w'.*exp(gammastar'*g_data);
                ws = ws/sum(ws);
                
                gammastar_mat(:,ct+1)=gammastar;
                % Check if new weights are okay 
                ok_nan = ~any(isnan(ws)); % ~any(isnan(ws)) returns logical 0 is there is any NaN in the weights and logical 1 otherwise
                % check if moment condition satisfied
                empirical_mom = g_data*ws';
                absmomdiff=abs(empirical_mom-gbar);
                maxabsmomdiff_probs=max(absmomdiff(1:numprobs));
                maxabsmomdiff_points=max(absmomdiff(numprobs+1:end));
                if isempty(maxabsmomdiff_points)
                    maxabsmomdiff_vec(1,ct+1) = maxabsmomdiff_probs;
                elseif isempty(maxabsmomdiff_probs)
                    maxabsmomdiff_vec(1,ct+1) = maxabsmomdiff_points;
                else
                    maxabsmomdiff_vec(1,ct+1) = maxabsmomdiff_probs + maxabsmomdiff_points;
                end

                if any(maxabsmomdiff_probs > prob_tol) || any(maxabsmomdiff_points > point_tol)
                    ok_momdiff = 0;
                else
                    ok_momdiff = 1;
                end

                ok_exitflag=exitflag>=1;
                ok=ok_nan*ok_momdiff*ok_exitflag;
                
                %fprintf("Number of NaNs: %d. Exitflag: %d. Maximum absolute moment difference: %1.2f \n",sum(isnan(ws)),exitflag,maxabsmomdiff);

                if ok==0 
                    %fprintf('applying more shrinkage \n');
                    ct = ct + 1;  % apply more shrinkage
                    if ct>20
                        break;
                    end

                end

        end
        
        % based on all runs, find the one leading to lowest maxabsmomdiff
        [~,min_idx]=min(maxabsmomdiff_vec);
        gammastar=gammastar_mat(:,min_idx);
        exitflag=exitflag_vec(1,min_idx);
        ws = w'.*exp(gammastar'*g_data);
        ws = ws/sum(ws);

        if min_idx>1
            finalshrink = shrink(min_idx-1);
        else
            finalshrink=0;
        end
        fprintf("Amount of shrinkage used: %1.2f \n", finalshrink);

end

% Check solution and issue warning if too far from goal
empirical_mom = g_data*ws';
absmomdiff=abs(empirical_mom-gbar);
maxabsmomdiff_probs=max(absmomdiff(1:numprobs));
maxabsmomdiff_points=max(absmomdiff(numprobs+1:end));
if any(maxabsmomdiff_probs > prob_tol) 
    warning('Check solution! PROB max. abs. diff %1.2f while tolerance is %1.2f \n',maxabsmomdiff_probs,prob_tol);
end
if any(maxabsmomdiff_points > point_tol) 
    warning('Check solution! POINT max. abs. diff %1.2f while tolerance is %1.2f \n',maxabsmomdiff_points,point_tol);
end


% Output
output.gamma = gammastar;
output.w = ws;
output.klic =sum(ws.*log(ws./w'));
output.shrink = finalshrink;  
output.convergence = exitflag;

end
 




