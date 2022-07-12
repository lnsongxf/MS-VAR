% Construct Local Mean in GDP Growth Rate
clear; close all;

load dln_gdp.mat

%@ Bi-Weight Parameter for local demeaning @
bw_bw = 100; 

%@ Parameters for ses @
nma = 4; 

% Compute Local Means
x = (dln_gdp);
[tmp,se] = bw_trend_se(x,bw_bw,nma);

% Plot
figure(1);
plot([x,tmp])
