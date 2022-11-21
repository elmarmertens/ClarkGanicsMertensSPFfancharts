function y = wprctile1(X,p,w,hasBeenSorted,type)
%WPRCTILE  Returns weighted percentiles of a sample with six algorithms.
% The idea is to give more emphasis in some examples of data as compared to
% others by giving more weight. For example, we could give lower weights to
% the outliers. The motivation to write this function is to compute percentiles
% for Monte Carlo simulations where some simulations are very bad (in terms of
% goodness of fit between simulated and actual value) than the others and to
% give the lower weights based on some goodness of fit criteria.
%
% USAGE:
%         y = WPRCTILE(X,p)
%         y = WPRCTILE(X,p,w)
%         y = WPRCTILE(X,p,w,hasBeenSorted,type)
%
% INPUT:
%    X -  vector or matrix of the sample data
%    p -  scalar  or a vector of percent values between 0 and 100
%
%    w -  positive weight vector (or matrix) for the sample data. Length of w must be
%         equal to either number of rows or columns of X. 
%         If X is matrix, operates on DIM=1; and, if w is vector 
%         saem weight vector w is used for all columns
%         (DIM=2). If the weights are equal, then WPRCTILE is same as PRCTILE.
%
%  type - an integer between 4 and 9 selecting one of the 6 quantile algorithms.
%         Type 4: p(k) = k/n. That is, linear interpolation of the empirical cdf.
%         Type 5: p(k) = (k-0.5)/n. That is a piecewise linear function where
%                 the knots are the values midway through the steps of the
%                 empirical cdf. This is popular amongst hydrologists. (default)
%                 PRCTILE also uses this formula.
%         Type 6: p(k) = k/(n+1). Thus p(k) = E[F(x[k])].
%                 This is used by Minitab and by SPSS.
%         Type 7: p(k) = (k-1)/(n-1). In this case, p(k) = mode[F(x[k])].
%                 This is used by S.
%         Type 8: p(k) = (k-1/3)/(n+1/3). Then p(k) =~ median[F(x[k])].
%                 The resulting quantile estimates are approximately
%                 median-unbiased regardless of the distribution of x.
%         Type 9: p(k) = (k-3/8)/(n+1/4). The resulting quantile estimates are
%                 approximately unbiased for the expected order statistics
%                 if x is normally distributed.
%
%         Interpolating between the points pk and X(k) gives the sample
%         quantile. Here pk is plotting position and X(k) is order statistics of
%         x such that x(k)< x(k+1) < x(k+2)...
%
% OUTPUT:
%    y - percentiles of the values in X
%        y is Np times Nx 
%
%
% Based on wprctile (matlab file exchange
% (https://www.mathworks.com/matlabcentral/fileexchange/16920-returns-weighted-percentiles-of-a-sample)
%
% Please note that this version of WPRCTILE will not work with NaNs values and
% planned to update in near future to handle NaNs values as missing values.
%
% References: Rob J. Hyndman and Yanan Fan, 1996, Sample Quantiles in Statistical
%             Package, The American Statistician, 50, 4.
%
% HISTORY:
% version 1.0.0, Release 2007/10/16: Initial release
% version 1.1.0, Release 2008/04/02: Implementation of other 5 algorithms and
%                                    other minor improvements of code
%
%
% I appreciate the bug reports and suggestions.
% See also: PRCTILE (Statistical Toolbox)

% Author: Durga Lal Shrestha
% UNESCO-IHE Institute for Water Education, Delft, The Netherlands
% eMail: durgals@hotmail.com
% Website: http://www.hi.ihe.nl/durgalal/index.htm
% Copyright 2004-2007 Durga Lal Shrestha.
% $First created: 16-Oct-2007
% $Revision: 1.1.0 $ $Date: 02-Apr-2008 13:40:29 $
%
% further changes made by Elmar Mertens (formatting, dim=2, hasBeenSorted)
%

% ***********************************************************************

%% check argument list

narginchk(2,5)

if isvector(X)
    X = X(:);
end

[ndraws, nvar] = size(X);

if nargin < 3 || isempty(w)
    % Default weight vector
    w = ones(ndraws,1);
end

% TODO: handle case where w is empty, but hasBeenSorted is true

if nargin < 4 || isempty(hasBeenSorted)
    hasBeenSorted = issorted(X,1);
end

if nargin < 5
    type = 5;
end

%% check inputs
if ~isvector(p) || numel(p) == 0
    error('wprctile:BadPercents', ...
        'P must be a scalar or a non-empty vector.');
elseif any(p < 0 | p > 100) || ~isreal(p)
    error('wprctile:BadPercents', ...
        'P must take real values between 0 and 100');
end
if ~ismatrix(X)
    error('wprctile:InvalidNumberofDimensions','X Must be 2D.')
end

if any(w<0)
    error('wprctile:InvalidWeight', ...
        'w values should be greater than 0');
end

% Check if there are NaN in any of the input
if  any(isnan(p))|| any(isnan(w), 'all')
    error('wprctile:NaNsInput',['This version of WPRCTILE does not work with ' ...
        'NaNs values in p or w.']);
end
%% Figure out which dimension WPRCTILE will work along using weight vector w

if isvector(w)
    if length(w) ~= ndraws
        error('wprctile:InvalidDimension', ...
            'length of w must be equal to draws of X');
    end
else
    if size(w) ~= size(X)
        error('wprctile:InvalidDimension', ...
            'if w is matrix, must be same size as X');
    end
end

%% prepare output matrix
np = length(p);
y  = zeros(np, nvar); 


%% sort data (if needed)
if hasBeenSorted
    sortedXmat = X;
    if isvector(w)
        sortedWmat = repmat(w, 1, nvar);
    else
        sortedWmat = w;
    end
else
    [sortedXmat, indSort] = sort(X,1); 
    % Change w to column vector
    w = w(:);
    sortedWmat = w(indSort); % note: creates matrix in dim==2 case, since w is vector
end

% normalise weight vector such that sum of the weight vector equals to n
sortedWmat = sortedWmat .* ndraws ./ sum(sortedWmat, 1);


%% Work on each column separately because of weight vector
for i=1:nvar
    sortedX = sortedXmat(:,i);
    sortedW = sortedWmat(:,i);              % rearrange the weight according to ind
    k = cumsum(sortedW);           % cumulative weight
    switch type                    % different algorithm to compute percentile
        case 4
            pk = k/ndraws;
        case 5
            pk = (k-sortedW/2)/ndraws;
        case 6
            pk = k/(ndraws+1);
        case 7
            pk = (k-sortedW)/(ndraws-1);
        case 8
            pk = (k-sortedW/3)/(ndraws+1/3);
        case 9
            pk = (k-sortedW*3/8)/(ndraws+1/4);
        otherwise
            error('wprctile:InvalidType', ...
                'Integer to select one of the six quantile algorithm should be between 4 to 9.')
    end

    % to avoid NaN for outside the range, the minimum or maximum values in X are
    % assigned to percentiles for percent values outside that range.
    q  = [0;pk;1];
    xx = [sortedX(1); sortedX; sortedX(end)];

    % Interpolation between q and xx for given value of p
    y(:,i) = interp1q(q,xx,p(:)./100);
end

