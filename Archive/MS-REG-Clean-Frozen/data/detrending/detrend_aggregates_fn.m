% Header
% This script detrends monthly aggregates and removes outliers.
% IMPORTANT: No NaNs can exist in middle of series.
%
% Writte by: Pablo Cuba-Borda, Jay Faris, and Dawson Miller
% Last modified: April 27, 2020. Washington, DC
%==========================================================================
clear; clc;

% Specify dataset structure
month_export_start = '1970m6';
end_date = '2020m12';

% Read data
data = readtable('tw_fn_data_for_ads_20200423_1138.csv');

% Map data
ip_us = data.data_ip;
rs_us = data.data_rs;
neo_us = data.data_neo;
gdp_us = data.data_gdp;
%ucl_us = data.data_uclaims;

% Compute quarterly average data
ip_us_q    = convert2Quarterly(ip_us,month_export_start,end_date);
rs_us_q    = convert2Quarterly(rs_us,month_export_start,end_date);
neo_us_q   = convert2Quarterly(neo_us,month_export_start,end_date);
gdp_us_q   = nan(length(ip_us_q),1);
%ucl_us_q   = convert2Quarterly(ucl_us,month_export_start,end_date);
gdp_us_tmp = gdp_us(~isnan(gdp_us));
gdp_us_q(1:length(gdp_us_tmp)) = gdp_us_tmp;

% Set monthly and quarterly dates
dates_q = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(6:7))/3,'Q'))';
dates_m = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(6:7))/3,'M'))';

%% Get trends and interpolate on monthly data from 1971m4 onwards
[m_rs_trend, m_rs_trend_full] = fn_local_detrend(rs_us_q,1,'Retail Sales');
[m_ip_trend, m_ip_trend_full] = fn_local_detrend(ip_us_q,1,'Industrial Production');
[m_pmi_trend, m_pmi_trend_full] = fn_local_detrend(neo_us_q,1,'New Export Orders');
[m_gdp_trend, m_gdp_trend_full] = fn_local_detrend(gdp_us_q,0,'GDP');
%[m_ucl_trend, m_ucl_trend_full] = fn_local_detrend(ucl_us_q,1,'Initial Claims / Labor Force');

% Detrend variables
detrended.rs.us  = rs_us  - m_rs_trend(3:end);
detrended.ip.us  = ip_us  - m_ip_trend(3:end);
detrended.neo.us = neo_us - m_pmi_trend(3:end);
detrended.gdp.us = gdp_us - m_gdp_trend(3:end);
%detrended.ucl.us = ucl_us - m_ucl_trend(3:end);

% Collect trend
trend.rs.us  = m_rs_trend;
trend.ip.us  = m_ip_trend;
trend.neo.us = m_pmi_trend;
trend.gdp.us = m_gdp_trend;
trend.gdp.us_full = m_gdp_trend_full;
%trend.ucl.us = m_ucl_trend;

% Remove outliers
detrended.rs.us_no_out   = removeOutliers(detrended.rs.us,4,1);
detrended.ip.us_no_out   = removeOutliers(detrended.ip.us,4,1);
detrended.neo.us_no_out  = removeOutliers(detrended.neo.us,4,1);
detrended.gdp.us_no_out  = removeOutliers(detrended.gdp.us,4,1);
%detrended.ucl.us_no_out  = removeOutliers(detrended.ucl.us,4,1);

%%
%  figure; clf; 
%  plot(detrended.ucl.us_no_out,'-o'); hold on;
%  plot(detrended.ucl.us); 
%%

out=[detrended.rs.us_no_out detrended.ip.us_no_out detrended.neo.us_no_out detrended.gdp.us_no_out];