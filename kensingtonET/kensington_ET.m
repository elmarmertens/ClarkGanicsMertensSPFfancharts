%% CGM, this script loads the raw CGM results and applies ET on them, then saves the weights
clear; close all; clc;

%% load toolboxes
path(pathdef)
addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/

%#ok<*UNRCH> 

%% settings, CGM results, fitted distributions

% load fitted values
datadir     = fullfile('..', 'kensingtonDataMatfiles');
FitCDFSPF_all = load(fullfile(datadir,'kensingtonFitDistSPF.mat'));

doTargetCurrentYear = false;

models = {'STATEtrendgapSV', 'STATEconst'};

DATALABELS = {'UNRATE', 'RGDP', 'PGDP'};

ET_types = {'binsOnly', 'binsAndMeans',...
            'NMeansOnly', 'NMeansAndVars',...   % normal means only, normal means and variances
            'GBMeansOnly', 'GBMeansAndVars',... % generalized beta means only, generalized beta means and variances
            'GBMeansVarsAndSkew',...            % generalized beta means, variances and skewness
            'GBVarsSkewAndSPFAnnPoints',...     % generalized beta variance and skewness while keeping the SPF annual point forecast as means             
            'SVBins',...                        % tilting CONST to SV probabilities assigned to SPF bins
            'binsOnlyNextYearOnly',...          % same as binsOnly, but only uses next-year SPF bins
            'GBMeansVarsAndSkewNextYearOnly'};  % same as GBMeansVarsAndSkew, but only uses next-year fitted GB moments

MCMCdraws = 3e3;

SPFprob = load(fullfile('..','kensingtonDataMatfiles','kensingtonPROB.mat')); % load all SPF density fcasts

opt_options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off','FiniteDifferenceType','central',...
    'OptimalityTolerance',1e-6,'FunctionTolerance',1e-6,'MaxIterations',1500,'HessianFcn','objective');

relESS_target = 0.4; % warning message will be issued if relESS falls below this threshold

startdates = {'RGDP', 'PGDP', 'UNRATE';
    datenum('1992Q1','yyyyQQ'), datenum('1992Q1','yyyyQQ'), datenum('2009Q2','yyyyQQ')};

%% start loop over models
for m = 1  : length(models)

    modellabel = models{m};

    %% start loop over variables
    for thisD = 1  : length(DATALABELS)

        datalabel = DATALABELS{thisD};
        SPFprob_act = SPFprob.(datalabel);

        % load raw CGM results
        CGMmcmc = matfile(fullfile(localresultsMCMC, sprintf('CGMmcmc-%s%s-Ndraws%d.mat', datalabel, modellabel, MCMCdraws)));
        Ndraws  = CGMmcmc.MCMCdraws * CGMmcmc.Nfedraws;
        init_w  = ones(Ndraws,1)/Ndraws; % initial weight vector, weights are uniform, so all 1/Ndraws
        Nz      = CGMmcmc.Nz;
        Nbar    = CGMmcmc.Nbar;
        dates   = CGMmcmc.dates;
        Zdata   = CGMmcmc.Zdata;
        doNIPA  = CGMmcmc.doNIPA;

        if doTargetCurrentYear
            max_num_fixed_event_fc = Nbar + 1;
        else
            max_num_fixed_event_fc = Nbar;
        end

        % load fitted CDFs
        thisCDFDATA = FitCDFSPF_all.(datalabel);

        % check if dates are the same in fitting results and MCMC results file
        if ~isequal(thisCDFDATA.dates,CGMmcmc.dates)
            warning('dates do NOT match between fitting results and MCMC results!');
        end

        % check if dates are the same in SPF probabilities and MCMC results file
        if ~isequal(SPFprob_act.dates,CGMmcmc.dates)
            warning('dates do NOT match between SPF probabilities and MCMC results!');
        end

        % find subsample where ET is applied
        firstT = find(dates==startdates{2,strcmpi(startdates(1,:),datalabel)});

        %firstT = max(find(~isnan(SPFprob_act.probs(:,1)),1,'first'),CGMmcmc.Tstart);
        lastT  = min(find(~isnan(SPFprob_act.probs(:,1)),1,'last'),CGMmcmc.T);
        numT   = lastT-firstT+1;

        first_last_Date      = dates([firstT;lastT]);
        first_last_DateLabel = datestr(first_last_Date, 'yyyyqq');

        YdensityDraws        = permute(CGMmcmc.YdensityDraws, [2 3 1]);
        YBARdensityDraws     = permute(CGMmcmc.YBARdensityDraws, [2 3 1]);
        YCURRENTdensityDraws = permute(CGMmcmc.YCURRENTdensityDraws, [2 3 1]);

        % save memory by loading only draws from firstT : lastT
        YdensityDraws        = YdensityDraws(:,:,firstT:lastT);
        YBARdensityDraws     = YBARdensityDraws(:,:,firstT:lastT);
        YCURRENTdensityDraws = YCURRENTdensityDraws(:,:,firstT:lastT);

        clear CGMmcmc

        %% start loop over ET types
        if strcmpi(datalabel,'RGDP')
            ET_types_idx_vec = 1:length(ET_types);
        else
            ET_types_idx_vec = 1:length(ET_types)-2;
        end

        for thisET = ET_types_idx_vec

                ETlabel = ET_types{thisET};
                
                if strcmpi(ETlabel,'NMeansOnly')
                    SPFfitdistlabel = 'N';
                    SPFETmom = 1;
                elseif strcmpi(ETlabel,'NMeansAndVars')
                    SPFfitdistlabel = 'N';
                    SPFETmom = 2;
                elseif strcmpi(ETlabel,'GBMeansOnly')
                    SPFfitdistlabel = 'B';
                    SPFETmom = 1;
                elseif strcmpi(ETlabel,'GBMeansAndVars')
                    SPFfitdistlabel = 'B';
                    SPFETmom = 2;
                elseif strcmpi(ETlabel,'GBMeansVarsAndSkew')
                    SPFfitdistlabel = 'B';
                    SPFETmom = 3;
                elseif strcmpi(ETlabel,'GBVarsSkewAndSPFAnnPoints')
                    SPFfitdistlabel = 'B';
                    SPFETmom = 3;
                elseif strcmpi(ETlabel,'GBMeansVarsAndSkewNextYearOnly')
                    SPFfitdistlabel = 'B';
                    SPFETmom = 3;
                end

                if strcmpi(ETlabel,'SVBins') && strcmpi(modellabel,'STATEconst')

                    % connect to SV model results, will need them when tilting CONST to SV
                    CGMmcmcSV = matfile(fullfile(localresultsMCMC, sprintf('CGMmcmc-%s%s-Ndraws%d.mat', datalabel, 'STATEtrendgapSV', MCMCdraws)));
                    
                    YBARdensityDrawsSV     = permute(CGMmcmcSV.YBARdensityDraws, [2 3 1]);
                    YCURRENTdensityDrawsSV = permute(CGMmcmcSV.YCURRENTdensityDraws, [2 3 1]);

                    % save memory by loading only draws from firstT : lastT
                    YBARdensityDrawsSV     = YBARdensityDrawsSV(:,:,firstT:lastT);
                    YCURRENTdensityDrawsSV = YCURRENTdensityDrawsSV(:,:,firstT:lastT);

                elseif strcmpi(ETlabel,'SVBins') && strcmpi(modellabel,'STATEtrendgapSV')
                    fprintf('skipping %s for %s and %s \n',datalabel,modellabel,ETlabel);
                    continue
                end

                if strcmpi(ETlabel,'binsOnly') || strcmpi(ETlabel,'binsAndMeans') || strcmpi(ETlabel,'SVBins') || strcmpi(ETlabel,'binsOnlyNextYearOnly')
                    ET_point_tol  = 0.5; % tolerance for discrepancy between SPF POINT forecasts and ET POINT forecasts
                    ET_prob_tol   = 1;   % tolerance for discrepancy between SPF PROBABILITY and ET PROBABILITY
                else
                    ET_point_tol = 0.01; % tolerance for discrepancy between SPF moments based on fitted CDF and ET moments
                    numprobs     = 0;
                    ET_prob_tol  = [];
                end
    
                % setup log file
                diaryfilename=fullfile(localstore, sprintf('%s-%s-ET%s-Ndraws%d.log', modellabel, datalabel, ETlabel, MCMCdraws));
                if exist(diaryfilename, 'file')
                    delete(diaryfilename);
                end
                diary(diaryfilename);
        
                if strcmpi(ETlabel,'binsOnly') || strcmpi(ETlabel,'binsAndMeans') || strcmpi(ETlabel,'SVBins') || strcmpi(ETlabel,'binsOnlyNextYearOnly')
        
                    [SPF_probs_vec_all,ET_annualprob_hz_all,keep_idx_all,ET_maxabsprobdiff_all,ET_absfixedhzdiff_all,ET_absfixedeventdiff_all,...
                        fixedhz_pointfc_all,fixedevent_pointfc_all,ET_fixedhz_pointfc_hz_all,ET_fixedevent_pointfc_hz_all] ...
                        = deal(cell(lastT,1));
        
                else
        
                    [ET_absmomdiff_all,ET_momcond_all] = deal(cell(lastT,1));
        
                end
        
                ESS_all  = NaN(lastT,1);
                ET_w_all = NaN(lastT,Ndraws);
                
                %% start loop over jump-offs
                for thisT=firstT:lastT
        
                    thisToffset = thisT - firstT + 1;
        
                    thisDate=dates(thisT,1);
                    thisDateLabel = datestr(thisDate, 'yyyyqq');
                    fprintf('********************************** \n');
                    fprintf('now doing %s in %s for %s and %s \n',datalabel,thisDateLabel,modellabel,ETlabel);
        
                    date_idx     = find(SPFprob_act.dates==thisDate); % find SPF round
                    SPF_numhz    = SPFprob_act.bins(date_idx,1); % SPF contains forecasts for this number of horizons
        
                    gridlength   = find(~isnan(SPFprob_act.bins(date_idx,2:end)),1,'last')/SPF_numhz; % number of bin edges
                    SPF_bins_all = SPFprob_act.bins(date_idx,2:2+gridlength*SPF_numhz-1);
        
                    % for how many CALENDAR years out we want to get forecasts beyond the current year, e.g. 3 corresponds to current year, next year, 2 and 3 years out
                    %fc_maxhz = Nbar;
        
                    % annual forecasts (calendar year growth rates for PGDP and RGDP, calendar year average level for UNRATE)
                    all_annual_fc_CGM      = cat(1, YCURRENTdensityDraws(:,:,thisToffset), YBARdensityDraws(:,:,thisToffset));
                    all_annual_fc_CGM_mean = mean(all_annual_fc_CGM,2);

                    % collecting SPF probabilities, NOT all of them HAVE to be used for ET
                    maxSPFhz      = SPFprob_act.bins(date_idx,1); % the SPF in this round contained density forecasts for this many horizons
                    SPF_probs_vec = NaN(maxSPFhz*(gridlength+1),1);
                    SPF_bins      = SPF_bins_all(1:gridlength);

                    if strcmpi(ETlabel,'SVBins') % extract probabilities implied by SV model for the same bins as SPF

                        % annual forecasts (calendar year growth rates for PGDP and RGDP, calendar year average level for UNRATE)
                        all_annual_fc_CGM_SV = cat(1, YCURRENTdensityDrawsSV(:,:,thisToffset), YBARdensityDrawsSV(:,:,thisToffset));
                        SV_probs_mat         = NaN(maxSPFhz,gridlength+1);

                        for hz = 1 : maxSPFhz
                            all_annual_fc_CGM_SV_tmp = all_annual_fc_CGM_SV(hz,1:Ndraws);
                            for bin_idx = 1 : gridlength + 1
                                if bin_idx == 1 % leftmost bin
                                    SV_probs_mat(hz,bin_idx) = 100 * sum(all_annual_fc_CGM_SV_tmp < SPF_bins(bin_idx)) / Ndraws;
                                elseif bin_idx == gridlength + 1 % rightmost bin
                                    SV_probs_mat(hz,bin_idx) = 100 * sum(all_annual_fc_CGM_SV_tmp > SPF_bins(bin_idx-1)) / Ndraws;
                                else % inner bins
                                    SV_probs_mat(hz,bin_idx) = 100 * sum( (all_annual_fc_CGM_SV_tmp > SPF_bins(bin_idx-1)) .* (all_annual_fc_CGM_SV_tmp < SPF_bins(bin_idx))) / Ndraws;
                                end
                            end
                        end

                        SPF_probs_vec = vec(SV_probs_mat'); % this is not SPF but SV

                        clear all_annual_fc_CGM_SV

                    else 
                         for hz = 1 : maxSPFhz
                            SPF_probs_vec((hz-1)*(gridlength+1)+1:hz*(gridlength+1),1)=fliplr(SPFprob_act.probs(date_idx,(hz-1)*(gridlength+1)+1:hz*(gridlength+1)))';
                         end
                    end
        
                    SPF_probs_vec_all{thisT}=SPF_probs_vec;
        
                    if doTargetCurrentYear
                        ET_annual_hz_vec            = 1:SPF_numhz; % 1 is current year, 2 is next, etc.
                    elseif contains(ETlabel,'NextYearOnly')
                        ET_annual_hz_vec            = 2;
                    else
                        ET_annual_hz_vec            = 2:SPF_numhz; % 1 is current year, 2 is next, etc.
                    end
                    ET_annualprob_hz_all{thisT} = ET_annual_hz_vec;
        
                    idx_all_mat=reshape(1:SPF_numhz*(gridlength+1),gridlength+1,SPF_numhz);
                    idx_all_vec=vec(idx_all_mat);
                    drop_idx=vec(idx_all_mat(:,setdiff(1:SPF_numhz,ET_annual_hz_vec))); % indices to be dropped due to horizons not selected for ET
        
                    if strcmpi(ETlabel,'binsOnly') || strcmpi(ETlabel,'binsAndMeans') || strcmpi(ETlabel,'SVBins') || strcmpi(ETlabel,'binsOnlyNextYearOnly')
        
                        % use g_function, which imposes 0-1 on the specific draws
                        g_data=NaN(SPF_numhz*(gridlength+1),Ndraws);
                        for EThz=ET_annual_hz_vec
                            g_data((EThz-1)*(gridlength+1)+1:EThz*(gridlength+1),:) = g_function(all_annual_fc_CGM(EThz,:),SPF_bins,true,[],[],[]);
                        end
        
                        % probabilities sum to 100%, so we are free to drop ONE bin at each horizon WLOG, but we DON'T have to
                        leftbins=0;
                        rightbins=1;
                        drop_left_idx=1:(gridlength+1):SPF_numhz*(gridlength+1);
                        drop_right_idx=(gridlength+1):(gridlength+1):SPF_numhz*(gridlength+1);
                        if leftbins==1 % discard highest bins
                            keep_idx=setdiff(idx_all_vec,[drop_right_idx(ET_annual_hz_vec),drop_idx']);
                        elseif rightbins==1 % discard lowest bins
                            keep_idx=setdiff(idx_all_vec,[drop_left_idx(ET_annual_hz_vec),drop_idx']);
                        else % don't discard any bins
                            keep_idx=setdiff(idx_all_vec,drop_idx');
                        end
            
                        keep_idx_all{thisT}=keep_idx;
                        numprobs = length(keep_idx);
            
                        % we CAN also add the information contained in the fixed-HORIZON point forecasts
                        fixedhz_pointfc = Zdata(thisT, 2:6)'; % expressed in %, so magnitude is similar to probabilities, which could help the numerical optimization
                        if doNIPA % convert from log-difference to simple growth rates
                             fixedhz_pointfc = (exp(fixedhz_pointfc / 100) - 1) * 100;
                        end
                        if strcmpi(ETlabel,'binsAndMeans')
                            fixedhz_pointfc_usehz = 1:5;
                        else
                            fixedhz_pointfc_usehz = [];
                        end
                        fixedhz_pointfc=fixedhz_pointfc(fixedhz_pointfc_usehz); % we can switch off tilting to the fixed-horizon point forecasts
                        fixedhz_pointfc_all{thisT}=fixedhz_pointfc;
                        ET_fixedhz_pointfc_hz_all{thisT}=fixedhz_pointfc_usehz;
            
                        % we CAN also add the information contained in the fixed-EVENT point forecasts starting with the NEXT calendar year
                        % TODO: Zdata contains NaN for next year in Q4 (which is OK if using fixedhz as well)
                        fixedevent_pointfc       = Zdata(thisT, 7:end)';
                        if doNIPA % convert from log-difference to simple growth rates
                             fixedevent_pointfc = (exp(fixedevent_pointfc / 100) - 1) * 100;
                        end
                        if strcmpi(ETlabel,'binsAndMeans')
                            fixedevent_pointfc_usehz = find(~isnan(fixedevent_pointfc));
                        else
                            fixedevent_pointfc_usehz = [];
                        end
                        fixedevent_pointfc       = fixedevent_pointfc(fixedevent_pointfc_usehz);
                        fixedevent_pointfc_all{thisT}       = fixedevent_pointfc;
                        ET_fixedevent_pointfc_hz_all{thisT} = fixedevent_pointfc_usehz;
            
                        % this is used for ET
                        if isempty(fixedhz_pointfc_usehz)
                            fixedhzYdensityDraws = [];
                        else
                            fixedhzYdensityDraws = YdensityDraws(fixedhz_pointfc_usehz,:,thisToffset);
                        end
                        if doNIPA % convert from log-difference to simple growth rates
                            fixedhzYdensityDraws = (exp(fixedhzYdensityDraws / 100) - 1) * 100;
                        end
            
                        g_data_ET = cat(1, g_data(keep_idx,:), fixedhzYdensityDraws, ...
                            all_annual_fc_CGM(1+fixedevent_pointfc_usehz,:));
        
                        clear g_data all_annual_fc_CGM fixedhzYdensityDraws
            
                        % moment conditions to be imposed, gbar
                        ET_momcond = cat(1, SPF_probs_vec(keep_idx), fixedhz_pointfc, ...
                            fixedevent_pointfc);
        
                    else
        
                        % compute mean (and variance) to be imposed in ET
                        if strcmpi(SPFfitdistlabel,'N')
        
                            hz_idx = 1;
                            ETmean = NaN(length(ET_annual_hz_vec),1);
                            ETvar  = NaN(length(ET_annual_hz_vec),1);
                            for EThz = ET_annual_hz_vec
                                thisNormalparams    = vec(thisCDFDATA.params_N(thisT,EThz,:));
                                ETmean(hz_idx,1)    = thisNormalparams(1,1);
                                ETvar(hz_idx,1)     = thisNormalparams(2,1)^2;
                                hz_idx = hz_idx + 1;
                            end
        
                        elseif strcmpi(SPFfitdistlabel,'B')
        
                            hz_idx = 1;
                            ETmean = NaN(length(ET_annual_hz_vec),1);
                            ETvar  = NaN(length(ET_annual_hz_vec),1);
                            ETskew = NaN(length(ET_annual_hz_vec),1);
                            for EThz = ET_annual_hz_vec
                                thisBetaparams      = vec(thisCDFDATA.params_B(thisT,EThz,:));
                                ETmean(hz_idx,1)    = thisBetaparams(1) + (thisBetaparams(2) - thisBetaparams(1)) * (thisBetaparams(3) / (thisBetaparams(3) + thisBetaparams(4)));
                                ETvar(hz_idx,1)     = thisBetaparams(3) * thisBetaparams(4) * (thisBetaparams(2) - thisBetaparams(1))^2 / ((thisBetaparams(3) + thisBetaparams(4))^2 * (thisBetaparams(3) + thisBetaparams(4) + 1));
                                ETskew(hz_idx,1)    = 2 * (thisBetaparams(4) - thisBetaparams(3)) * sqrt(thisBetaparams(3) + thisBetaparams(4) + 1) / ((thisBetaparams(3) + thisBetaparams(4) + 2) * sqrt(thisBetaparams(3) * thisBetaparams(4)));
                                hz_idx = hz_idx + 1;
                            end
                            
                        end

                        if contains(ETlabel,'SPFAnnPoints') % use model-implied means as SPF annual point forecasts in ET
                               
                               hz_idx = 1;
                               for EThz = ET_annual_hz_vec
                                   ETmean(hz_idx,1) = all_annual_fc_CGM_mean(EThz,1);
                                   hz_idx = hz_idx + 1;
                               end

                        end
                        
                        % use g_function which imposes the "standard" ET moment conditions, mean or mean and variance or mean and variance and skewness
                        g_data_ET=NaN((hz_idx-1)*SPFETmom,Ndraws);
                        hz_idx = 1;
                        for EThz=ET_annual_hz_vec
                            g_data_ET((hz_idx-1)*SPFETmom+1:hz_idx*SPFETmom,:) = g_function(all_annual_fc_CGM(EThz,:),[],false,SPFETmom,ETmean(hz_idx),ETvar(hz_idx));
                            hz_idx = hz_idx + 1;
                        end
        
                        clear all_annual_fc_CGM
        
                        if SPFETmom == 1
                            ET_momcond = ETmean;
                        elseif SPFETmom == 2
                            ET_momcond = vec([ETmean, ETvar]');
                        elseif SPFETmom == 3
                            ET_momcond = vec([ETmean, ETvar, ETskew]');
                        end

                        ET_momcond_all{thisT} = ET_momcond;
        
                    end
        
                    % do ET, add moment conditions one by one, using previous gammastar as starting value
                    % potentially more robust than doing the optimization over a large vector
                    % now this is not active, but we can easily turn this back on if we
                    % encounter a numerical issue
                    gamma_init=zeros(size(ET_momcond,1),1);
                    for num_mom_tmp=size(ET_momcond,1)
                        tilt_output = tilt_mod(g_data_ET(1:num_mom_tmp,1:Ndraws),init_w,gamma_init,ET_momcond(1:num_mom_tmp),ET_prob_tol,ET_point_tol,numprobs,opt_options);
                        gamma_init=[tilt_output.gamma;0];
                    end
        
                    tiltmom_ET=g_data_ET*tilt_output.w'; % this calculates the moments based on ET
                    %tiltmom=g_data*tilt_output.w'; % this calculates the HISTOGRAM probabilities based on ET
                    momdiff=tiltmom_ET-ET_momcond;
                    abs_momdiff=abs(momdiff);
        
                    if strcmpi(ETlabel,'binsOnly') || strcmpi(ETlabel,'binsAndMeans') || strcmpi(ETlabel,'SVBins') || strcmpi(ETlabel,'binsOnlyNextYearOnly')
        
                        absmomdiff_prob=abs_momdiff(1:length(keep_idx)); % absolute moment difference for probabilities
                        absmomdiff_fixedhz=abs_momdiff(length(keep_idx)+1:length(keep_idx)+length(fixedhz_pointfc)); % absolute moment difference for fixed-HORIZON forecasts
                        absmomdiff_fixedevent=abs_momdiff(length(keep_idx)+length(fixedhz_pointfc)+1:end); % absolute moment difference for fixed-HORIZON forecasts
                        ETnumbins_annual=length(absmomdiff_prob)/length(ET_annual_hz_vec); % number of bins used in each year
                        max_absmomdiff_prob_vec_tmp=NaN(length(ET_annual_hz_vec),1);
                        for EThz=1:length(ET_annual_hz_vec)
                            max_absmomdiff_prob_vec_tmp(EThz,1)=max(absmomdiff_prob(1+(EThz-1)*ETnumbins_annual:EThz*ETnumbins_annual,1));
                        end
                        max_absmomdiff_prob_vec=NaN(SPF_numhz,1);
                        max_absmomdiff_prob_vec(ET_annual_hz_vec)=max_absmomdiff_prob_vec_tmp; % at each calendar year hz, save mean absolute moment difference for probabilities
                        ET_maxabsprobdiff_all{thisT}=max_absmomdiff_prob_vec;
                        absmomdiff_fixedhz_vec=NaN(5,1);
                        absmomdiff_fixedhz_vec(fixedhz_pointfc_usehz)=absmomdiff_fixedhz;
                        ET_absfixedhzdiff_all{thisT}=absmomdiff_fixedhz_vec;
                        absmomdiff_fixedevent_vec=NaN(max_num_fixed_event_fc,1);
                        absmomdiff_fixedevent_vec(fixedevent_pointfc_usehz)=absmomdiff_fixedevent;
                        ET_absfixedeventdiff_all{thisT}=absmomdiff_fixedevent_vec;

                    else

                        ET_absmomdiff_all{thisT} = abs_momdiff;

                    end

                    ESS=1/sum(tilt_output.w.^2);
                    if ESS/Ndraws<relESS_target
                        warning('relESS is %1.2f at %d for %s ! \n',ESS/Ndraws,thisT,datalabel);
                    end
                    ESS_all(thisT)=ESS;
        
                    ET_w_all(thisT,1:Ndraws)=tilt_output.w;
        
                    fprintf('DONE WITH %s in %s for %s and %s \n',datalabel,thisDateLabel,modellabel,ETlabel);
                    fprintf('********************************** \n');
                end % jump offs
    
            %% here we save the weights, indices and ESS
            if strcmpi(ETlabel,'binsOnly') || strcmpi(ETlabel,'binsAndMeans') || strcmpi(ETlabel,'SVBins') || strcmpi(ETlabel,'binsOnlyNextYearOnly')
    
            save(fullfile(localstore,sprintf('kensington-ET%s-Ndraws%d-%s%s.mat', ETlabel, MCMCdraws, datalabel, modellabel)),...
                'ET_w_all','ESS_all','ET_maxabsprobdiff_all','ET_absfixedhzdiff_all','ET_absfixedeventdiff_all',...
                'fixedhz_pointfc_all','fixedevent_pointfc_all','SPF_probs_vec_all','keep_idx_all',...
                'ET_annualprob_hz_all','ET_fixedhz_pointfc_hz_all','ET_fixedevent_pointfc_hz_all',...
                'first_last_Date','first_last_DateLabel', ...
                'firstT', 'lastT', ...         
                '-v7.3');
    
            else
    
                save(fullfile(localstore,sprintf('kensington-ET%s-Ndraws%d-%s%s.mat', ETlabel, MCMCdraws, datalabel, modellabel)),...
                'ET_w_all','ESS_all','ET_absmomdiff_all','ET_momcond_all',...
                'first_last_Date','first_last_DateLabel', ...
                'firstT', 'lastT', ...         
                '-v7.3');
    
            end

            diary off

        end % ET types
    end % DATALABEL
end % models