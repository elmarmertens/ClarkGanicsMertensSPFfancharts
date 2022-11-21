function [g_data] = g_function(rawdata,SPF_bins,doBins,ETmom,ETmeans,ETvars)
% this function takes the data
% and applies the g function (Krueger, Clark and Ravazzolo, 2018, JBES)
D=size(rawdata,2);
if doBins % moment conditions are probabilities from histograms

    num_momcond=size(SPF_bins,2)+1;
    tmp=rawdata<SPF_bins(1);
    g_data=NaN(num_momcond,D);
    g_data(1,1:D) = 100 * tmp;
    
    for b=2:num_momcond-1
        g_data(b,1:D) = 100 * ((rawdata>=SPF_bins(b-1)).*(rawdata<SPF_bins(b))); % returns 1 at positions where draw is inside the bounds, 0 otherwise
    end
    g_data(num_momcond,1:D) = 100 * (rawdata>=SPF_bins(end));

else % "standard", either on mean or mean and variance

    switch ETmom

        case 1 % only mean restrictions
    
            g_data = rawdata;
    
        case 2 % mean and variance restrictions
    
            means_mat = repmat(ETmeans,1,D);
            g_data    = [rawdata;(rawdata-means_mat).^2];

        case 3 % mean, variance and skewness restrictions

            means_mat = repmat(ETmeans,1,D);
            std_mat   = repmat(sqrt(ETvars),1,D);
            g_data    = [rawdata;(rawdata-means_mat).^2;((rawdata-means_mat)./std_mat).^3];
            
    end
    
end

end

