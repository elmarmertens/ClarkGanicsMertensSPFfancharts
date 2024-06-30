function Zmvlogscore = calculateLogscore(Ytp1var, Ctp1, Ztp1meanerror, Ztp1varnoise)

    if nargin < 4
        Ztp1varnoise = [];
    end

Nfcstdraws       = size(Ytp1var, 3) * size(Ytp1var, 4);
logscoredraws    = NaN(Nfcstdraws, 1);

[Nztp1, Nstates] = size(Ctp1);
if Nztp1 > Nstates 
    Nztp1         = Nstates;
    Ctp1          = Ctp1(1:Nztp1, :);
    Ztp1meanerror = Ztp1meanerror(1:Nztp1,:);
    if ~isempty(Ztp1varnoise)
        Ztp1varnoise = Ztp1varnoise(1:Nztp1, :, :);
    end
end
Nztp1xlogtwopi = Nztp1 * log(2 * pi);

ii = 0;
for mm = 1 : size(Ytp1var, 4)
    for jj = 1 : size(Ytp1var, 3)
        ii = ii + 1;
        Zvar      = Ctp1 * Ytp1var(:, :, jj, mm) * Ctp1';
        if ~isempty(Ztp1varnoise)
            Zvar = Zvar + diag(Ztp1varnoise(:, jj, mm));
        end
        sqrtZvar  = chol(Zvar, 'lower');
        % [sqrtZvar, flag]  = chol(Zvar, 'lower');
        % if flag ~= 0
        %     sqrtYtp1var = chol(Ytp1var(:, :, jj, mm), "lower");
        %     sqrtZvar    = cholqr(Ctp1 * sqrtYtp1var)
        % end
        logdetvar         = 2 * sum(log(diag(sqrtZvar)));
        Ztilde            = sqrtZvar \ Ztp1meanerror(:, mm);
        logscoredraws(ii) = - .5 * (Nztp1xlogtwopi + logdetvar + sum(Ztilde.^2));
    end
end
maxlogscoredraw = max(logscoredraws);
Zmvlogscore     = log(mean(exp(logscoredraws - maxlogscoredraw))) + maxlogscoredraw;