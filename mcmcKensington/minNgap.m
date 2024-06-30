function Ngap = minNgap(datalabel, thisDate)
% minNgap: find the minimal Ngap setting for SPF variable datalabel
%
% USAGE: Ngap = minNgap(datalabel, thisDate)
%


%#ok<*DATNM>
if nargin < 2
    thisDate = [];
end

if isempty(thisDate) % full sample choice
    switch upper(datalabel)
        case 'RGDP'
            Ngap = 8;
        case 'PGDP'
            Ngap = 6;
        case {'UNRATE', 'TBILL'}
            Ngap = 10;
        case 'CPI'
            Ngap = 7;
        otherwise
            error('datalabel <<%s>> not recognized', datalabel);
    end
else % group into two buckets (ignoring pre 1981 simplifications)
    switch upper(datalabel)
        case 'RGDP'
            if thisDate < datenum(2009,4,1)
                Ngap = 6;
            else
                Ngap = 8;
            end
        case 'PGDP'
            Ngap = 6;
        case {'UNRATE', 'TBILL'}
            if thisDate < datenum(2009,4,1) 
                Ngap = 6;
            else
                Ngap = 10;
            end
        case 'CPI'
            if thisDate < datenum(2005,10,1)
                Ngap = 6;
            else
                Ngap = 7;
            end
        otherwise
            error('datalabel <<%s>> not recognized', datalabel);
    end
end
