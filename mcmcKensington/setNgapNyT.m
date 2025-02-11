function Ngap = setNgapNyT(datalabel, thisDate)
% setNgapNyT: set Ngap to minimal NyT for a given sample
%
% USAGE: Ngap = setNgapNyT(datalabel, thisDate)
%


%#ok<*DATNM>
if nargin < 2
    thisDate = [];
end

if isempty(thisDate) % full sample choice
    switch upper(datalabel)
        case {'RGDP','UNRATE'}
            Ngap = 17;
        case 'PGDP'
            Ngap = 9;
        case 'CPI'
            Ngap = 13;
        otherwise
            error('datalabel <<%s>> not recognized', datalabel);
    end
else % group into two buckets (ignoring pre 1981 simplifications)
    switch upper(datalabel)
        case {'RGDP','UNRATE'}
            if thisDate < datenum(2009,4,1)
                Ngap = 9;
            else
                Ngap = setNgapNyT(datalabel);
            end
        case 'PGDP'
            Ngap = setNgapNyT(datalabel);
        case 'CPI'
            if thisDate < datenum(2005,10,1)
                Ngap = 9;
            else
                Ngap = setNgapNyT(datalabel);
            end
        otherwise
            error('datalabel <<%s>> not recognized', datalabel);
    end
end
