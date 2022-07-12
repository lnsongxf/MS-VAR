function [xmean, se] = bw_trend_se(x,bw_bw,nma)


xmean = missmat(rows(x),1);
z = seqa(1,1,rows(x));
tmp = packr([x z]);
x = tmp(:,1);
z = tmp(:,2);
nobs = rows(x);
mtmp = zeros(nobs,1);
trend = seqa(1,1,nobs);

for t=1:nobs
    dt = (trend-t)/bw_bw;
    bw_weight = (15/16)*((1-dt.^2).^2);   % Bi-Weight
    bw_weight = bw_weight.*(abs(dt)< 1);
    bw_weight = bw_weight/sum(bw_weight);
    mtmp(t)=bw_weight'*x;
end


% for ii = 1:rows(z)
%    xmean(z(ii)) = mtmp(z(ii)); 
% end
xmean(z(1):z(rows(z)))= mtmp; 
% Compute SEs
x = x - mtmp;

% Compute ACVs
acv = missmat(nma+1,1);

i = 1;
while i <= nma+1
    acv(i) = x(i:nobs)'*x(1:nobs+1-i)/(nobs+1-i);
    i = i+1;
end

var_tmp = missmat(nobs,1);
se = missmat(rows(xmean),1);

for t=1:nobs
    dt = (trend-t)/bw_bw;
    bw_weight = (15/16)*((1-dt.^2).^2);   % Bi-Weight
    bw_weight = bw_weight.*(abs(dt)< 1);
    bw_weight = bw_weight/sum(bw_weight);
    var_tmp(t) = acv(1)*(bw_weight'*bw_weight);
    i = 2;
    while i <= nma+1
        var_tmp(t) = var_tmp(t) + 2*acv(i)*(bw_weight(i:rows(bw_weight))'*bw_weight(1:rows(bw_weight)+1-i));
        i = i+1;
    end
end

se(z(1):z(rows(z)))=sqrt(var_tmp);








