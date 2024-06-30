function YAVGdraws = constructYAVGdraws(T, Kavg, Tstart, dates, datalabel, RTdata, YdensityDraws)
    ylags = NaN(T,1,Kavg-1);
    parfor thisT = Tstart : T
        thisDate      = dates(thisT);
        thisDateLabel = datestr(thisDate, 'yyyyqq');
        switch upper(datalabel)
            case {'UNRATE', 'CPI'}
                % use latest vintage, not real-time data
                RT_vintage_idx = length(RTdata.vindates); % pick last vintage
            case {'RGDP', 'PGDP'}
                RT_vintage_idx = find(datenum(RTdata.vindates) == thisDate,1,'first'); % vintage, time of SPF round
            otherwise
                RT_vintage_idx = []; %#ok<NASGU> % to avoid parfor warning
                error('datalabel <<%s>> not supported', datalabel)
        end
        % last observation available at time of SPF round
        RT_obsdate_idx = find(datenum(RTdata.obsdates)==datenum(dateshift(datetime(thisDateLabel,'InputFormat','yyyyQQQ'),'start','quarter',-1)),1,'first');
        % vector of vintage data
        RT_vec = RTdata.data(1:RT_obsdate_idx,RT_vintage_idx);
        
        % collect three lags, note that RT_obsdate_idx is one lag relative to thisT
        ylags(thisT,1,:) = RT_vec(end-(Kavg-1)+1:end); % store lag three first
    end
    
    % collect draws
    ydraws    = permute(YdensityDraws, [1 3 2]);
    YAVGdraws = NaN(size(ydraws));
    % h < Kavg
    for h = 1 : (Kavg - 1)
        YAVGdraws(:,:,h) = (sum(ydraws(:,:,1:h), 3) + sum(ylags(:,1,end-(Kavg-1)+h:end), 3)) / Kavg;
    end
    % h >= Kavg
    for h = Kavg : size(YAVGdraws, 3)
        YAVGdraws(:,:,h) = sum(ydraws(:,:,h+(-(Kavg-1):0)), 3) / Kavg;
    end
    YAVGdraws = permute(YAVGdraws, [1 3 2]);
end