%--------------------------------------------------------------------------
% % Housekeeping
%--------------------------------------------------------------------------
clear all;
clc;
close all;

%--------------------------------------------------------------------------
% % Load useful libraries
%--------------------------------------------------------------------------

addpath('td'); % load Enrique M. Quilis Temporal Disaggregation Library


%--------------------------------------------------------------------------
% % Interpolate real GDP
%--------------------------------------------------------------------------

% read real GDP

raw_data_GDPC1 = readtable('csvfiles/GDPC1.csv');
% GDPC1          = table2array(raw_data_GDPC1(find(strcmp(raw_data_GDPC1.DATE,'1965-01-01')==1):find(strcmp(raw_data_GDPC1.DATE,'2007-10-01')==1),2));
GDPC1          = table2array(raw_data_GDPC1(find(raw_data_GDPC1.DATE=='1965-01-01'):find(raw_data_GDPC1.DATE=='2007-10-01'),2));


% read industrial production
raw_data_IP = readtable('csvfiles/INDPRO.csv');
% IP          = table2array(raw_data_IP(find(strcmp(raw_data_IP.DATE,'1965-01-01')==1):find(strcmp(raw_data_IP.DATE,'2007-12-01')==1),2));
IP          = table2array(raw_data_IP(find(raw_data_IP.DATE=='1965-01-01'):find(raw_data_IP.DATE=='2007-12-01'),2));


Y  = GDPC1; % Y: Nx1 ---> vector of low frequency data
x  = IP;    % x: nxp ---> matrix of high frequency indicators (without intercept)
ta = 2;     % type of disaggregation ---> average (index)
sc = 3;     % quarterly to monthly
type = 0;   % estimation method:  (0) weighted least squares  (1) maximum likelihood
opC  = 1;   % no intercept in hf model

% res = chowlin(Y,x,ta,sc,type,opC,[]);
res = chowlin_td(Y,x,ta,sc,type,opC,[]);

monthly_GDP = res.y;


%--------------------------------------------------------------------------
% % Interpolate real GDP deflator
%--------------------------------------------------------------------------

% read GDP Deflator
raw_data_GDPDEF = readtable('csvfiles/GDPDEF.csv');
% GDPDEF          = table2array(raw_data_GDPDEF(find(strcmp(raw_data_GDPDEF.DATE,'1965-01-01')==1):find(strcmp(raw_data_GDPDEF.DATE,'2007-10-01')==1),2));
GDPDEF          = table2array(raw_data_GDPDEF(find(raw_data_GDPDEF.DATE=='1965-01-01'):find(raw_data_GDPDEF.DATE=='2007-10-01'),2));


% read CPIAUCSL
raw_data_CPIAUCSL = readtable('csvfiles/CPIAUCSL.csv');
% CPIAUCSL          = table2array(raw_data_CPIAUCSL(find(strcmp(raw_data_CPIAUCSL.DATE,'1965-01-01')==1):find(strcmp(raw_data_CPIAUCSL.DATE,'2007-12-01')==1),2));
CPIAUCSL          = table2array(raw_data_CPIAUCSL(find(raw_data_CPIAUCSL.DATE=='1965-01-01'):find(raw_data_CPIAUCSL.DATE=='2007-12-01'),2));



% read PPIFGS
raw_data_PPIFGS = readtable('csvfiles/PPIFGS.csv');
% PPIFGS          = table2array(raw_data_PPIFGS(find(strcmp(raw_data_PPIFGS.DATE,'1965-01-01')==1):find(strcmp(raw_data_PPIFGS.DATE,'2007-12-01')==1),2));
PPIFGS          = table2array(raw_data_PPIFGS(find(raw_data_PPIFGS.DATE=='1965-01-01'):find(raw_data_PPIFGS.DATE=='2007-12-01'),2));



Y  = GDPDEF;               % Y: Nx1 ---> vector of low frequency data
x  = [CPIAUCSL,PPIFGS];    % x: nxp ---> matrix of high frequency indicators (without intercept)
ta = 2;     % type of disaggregation ---> average (index)
sc = 3;     % quarterly to monthly
type = 0;   % estimation method:  (0) weighted least squares  (1) maximum likelihood
opC  = 1;   % no intercept in hf model

% res = chowlin(Y,x,ta,sc,type,opC,[]);
res = chowlin_td(Y,x,ta,sc,type,opC,[]);

monthly_GDPDEF = res.y;


%--------------------------------------------------------------------------
% % Commodity Price Index: Obtained from Global Financial Data through the
%   Board of Governors of the Federal Reserve System
%--------------------------------------------------------------------------
raw_data_CPRINDEX = readtable('csvfiles/DJSD_20170720.csv');

% CPRINDEX          = table2array(raw_data_CPRINDEX(find(strcmp(raw_data_CPRINDEX.Date,'01/31/1965')==1):find(strcmp(raw_data_CPRINDEX.Date,'12/31/2007')==1),9));
CPRINDEX          = table2array(raw_data_CPRINDEX(find(raw_data_CPRINDEX.Date=='01/31/1965'):find(raw_data_CPRINDEX.Date=='12/31/2007'),9));


%--------------------------------------------------------------------------
% % Total reserves: Adjusted for Changes in Reserve Requirements
%--------------------------------------------------------------------------
raw_data_TRARR = readtable('csvfiles/TRARR.csv');
% TRARR          = table2array(raw_data_TRARR(find(strcmp(raw_data_TRARR.DATE,'1965-01-01')==1):find(strcmp(raw_data_TRARR.DATE,'2007-12-01')==1),2));
TRARR          = table2array(raw_data_TRARR(find(raw_data_TRARR.DATE=='1965-01-01'):find(raw_data_TRARR.DATE=='2007-12-01'),2));

%--------------------------------------------------------------------------
% % Non-borrowed reserves
%--------------------------------------------------------------------------
raw_data_BOGNONBR = readtable('csvfiles/BOGNONBR.csv');
% BOGNONBR          = table2array(raw_data_BOGNONBR(find(strcmp(raw_data_BOGNONBR.DATE,'1965-01-01')==1):find(strcmp(raw_data_BOGNONBR.DATE,'2007-12-01')==1),2));
BOGNONBR          = table2array(raw_data_BOGNONBR(find(raw_data_BOGNONBR.DATE=='1965-01-01'):find(raw_data_BOGNONBR.DATE=='2007-12-01'),2));

%--------------------------------------------------------------------------
% % Federal Funds Rate
%--------------------------------------------------------------------------
raw_data_FEDFUNDS = readtable('csvfiles/FEDFUNDS.csv');
% FEDFUNDS          = table2array(raw_data_FEDFUNDS(find(strcmp(raw_data_FEDFUNDS.DATE,'1965-01-01')==1):find(strcmp(raw_data_FEDFUNDS.DATE,'2007-12-01')==1),2));
FEDFUNDS          = table2array(raw_data_FEDFUNDS(find(raw_data_FEDFUNDS.DATE=='1965-01-01'):find(raw_data_FEDFUNDS.DATE=='2007-12-01'),2));




%--------------------------------------------------------------------------
% % Export data
%--------------------------------------------------------------------------

% dates = table2array(raw_data_FEDFUNDS(find(strcmp(raw_data_FEDFUNDS.DATE,'1965-01-01')==1):find(strcmp(raw_data_FEDFUNDS.DATE,'2007-12-01')==1),1));
dates = table2array(raw_data_FEDFUNDS(find(raw_data_FEDFUNDS.DATE=='1965-01-01'):find(raw_data_FEDFUNDS.DATE=='2007-12-01'),1));

dataset = table(dates,monthly_GDP,monthly_GDPDEF,CPRINDEX,TRARR,BOGNONBR,FEDFUNDS,IP,CPIAUCSL);

writetable(dataset,'csvfiles/dataset.csv')





