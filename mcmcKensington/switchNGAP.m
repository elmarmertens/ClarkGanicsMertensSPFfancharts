function [Ngap, NGAP] = switchNGAP(NGAP, Ny, datalabel, thisDate)
% SWITCHNGAP ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 11-Jan-2024 18:06:37 $
% $Revision : 1.00 $
% DEVELOPED : 23.2.0.2459199 (R2023b) Update 5
% FILENAME  : switchNGAP.m

narginchk(1, 4);

if nargin < 2
    Ny = 17;
end
if nargin < 3
    datalabel = 'RGDP';
end
if nargin < 4
    thisDate = [];
end

switch NGAP
    case 'Ny'
        Ngap = Ny;
    case 'NyT'
        Ngap = setNgapNyT(datalabel,thisDate);
    case 'NyLess1'
        Ngap = Ny - 1;
    case 'NyLess3'
        Ngap = Ny - 3;
    case {'MinState'}
        Ngap = minNgap(datalabel); 
    case {'MinStateT'}
        Ngap = minNgap(datalabel,thisDate); 
    case {'BOP'}
        Ngap = setNgapBOP(datalabel,thisDate); 
    otherwise % assumes NGAP is numeric
        if isnumeric(NGAP)
            Ngap = NGAP;
            NGAP = 'Ngap';
        else
            error('Ngap choice of <<%s>> not recognized',NGAP)
        end
end