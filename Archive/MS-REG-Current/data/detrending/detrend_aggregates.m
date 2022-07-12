%% Header
% This script detrends monthly aggregates and removes outliers.
% IMPORTANT: No NaNs can exist in middle of series.
%
% Writte by: Pablo Cuba-Borda, Jay Faris, and Dawson Miller
% Last modified: April 27, 2020. Washington, DC
%==========================================================================
% clear; clc;
%
% % Specify dataset structure
% month_export_start = '1970m6';
% end_date = '2020m12';
%
% % Read data
% data = readtable('ppp_us_data_for_ads_20200423_0820.csv');

% Data vintage, sample and country selection
datafilename = '11302020';
sheetuse     = 'US_IPEBPMARKIT';
start_date   = '1973-Jan';
end_date     = '2020-Oct';

% Data bounds in terms of quarter at monthly frequency
month_export_start = '1973M3';
end_date = '2020M12';
bw_bw    = 40;

% Read data
data = readtable(['../' datafilename '.xlsx'],'Sheet',sheetuse);

% Map data
monthly_data  = data.GDPG;

% Compute quarterly average data
data_q    = convert2Quarterly(monthly_data,month_export_start,end_date);


% Censor data before detrending
% Avoids outliers to influece trend
[~,data_q_censored]   = removeOutliers(data_q,4,1);


% Set monthly and quarterly dates
dates_q = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(6:7))/3,'Q'))';
dates_m = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(6:7))/3,'M'))';

%% Get trends and interpolate on monthly data from 1971m4 onwards
[m_trend, m_trend_full]   = fn_local_detrend(data_q_censored,1,'GDP data',bw_bw);

% Detrend variables
detrended.data.raw  = monthly_data  - m_trend;

% Collect trend
trend.q = m_trend;
trend.full = m_trend_full;


