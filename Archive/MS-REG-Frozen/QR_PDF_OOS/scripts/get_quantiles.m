function [YQ, YQunc, Yh, BQ, BQunc] = get_quantiles(y,x_mat,hh,quantiles)

% Compute moving average for dependent variable
if hh==0
    Yh = y;   
else
    Yh = filter(ones(1,hh)/hh, 1, y);
    Yh(1:hh-1) = NaN;
end

%% Conditional Quantiles
% RHS
Z = [ones(size(x_mat,1),1), x_mat];
% Cut LHS and RHS to forecast horizon
yh = Yh(hh+1:end); % First value for hh=12: Jan-1987, i.e. mean(Jan-1986 to Jan 1987)
Zh = Z(1:end-hh,:); % First value: Jan-1986
% Run Quantile regressions and get predicted quantiles
BQ = NaN(size(Z,2), length(quantiles));
YQ = NaN(size(Zh,1), length(quantiles));
for jq = 1:length(quantiles)
    BQ(:, jq) = rq(Zh, yh, quantiles(jq));
    YQ(:, jq) = Zh * BQ(:, jq);
end

%% Unconditional Quantiles
% RHS
Z = ones(size(x_mat,1),1);
% Cut LHS and RHS to forecast horizon
yh = Yh(hh+1:end); % First value for hh=12: Jan-1987, i.e. mean(Jan-1986 to Jan 1987)
Zh = Z(1:end-hh,:); % First value: Jan-1986
% Run Quantile regressions and get predicted quantiles
BQunc = NaN(size(Z,2), length(quantiles));
YQunc = NaN(size(Zh,1), length(quantiles));
for jq = 1:length(quantiles)
    BQunc(:, jq) = rq(Zh, yh, quantiles(jq));
    YQunc(:, jq) = Zh * BQunc(:, jq);
end
    
end

