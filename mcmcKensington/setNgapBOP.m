function Ngap = setNgapBOP(datalabel, thisDate)
% setNgapNyT: set Ngap to minimal beginning of period for longest forecast horizon
%
% USAGE: Ngap = setNgapNyT(datalabel, thisDate)
%


%#ok<*DATNM>

if isempty(thisDate)
    thisDate = datenum(2023,10,1);
end

switch upper(datalabel)
    case {'RGDP','UNRATE'}
        if thisDate < datenum(2009,4,1)
            Ngap = 7;
        else
            Ngap = 14;
        end
    case {'PGDP', 'UNRATEY1'}
        Ngap = 7;
    case 'CPI'
        if thisDate < datenum(2005,10,1)
            Ngap = 7;
        else
            Ngap = 10;
        end
    otherwise
        error('datalabel <<%s>> not recognized', datalabel);
end
