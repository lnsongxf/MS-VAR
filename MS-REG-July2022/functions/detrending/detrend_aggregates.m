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

% Map data
ip_data  = data.data_ip;
rs_data  = data.data_rs;
neo_data = data.data_neo;
gdp_data = data.data_gdp;
if (ind_us==1 | ind_us==2) & foreign==0
    ucl_data = data.data_uclaims;
end

% Compute quarterly average data
ip_q    = convert2Quarterly(ip_data,month_export_start,end_date);
rs_q    = convert2Quarterly(rs_data,month_export_start,end_date);
neo_q   = convert2Quarterly(neo_data,month_export_start,end_date);
gdp_q   = nan(length(ip_q),1);
gdp_tmp = gdp_data(~isnan(gdp_data));
gdp_q(1:length(gdp_tmp)) = gdp_tmp;
if (ind_us==1 | ind_us==2) & foreign==0
    ucl_q   = convert2Quarterly(ucl_data,month_export_start,end_date);
end


% Censor data before detrending
% Avoids outliers to influece trend
[~,rs_q_censored]   = removeOutliers(rs_q,4,1);
[~,ip_q_censored]   = removeOutliers(ip_q,4,1);
[~,neo_q_censored]  = removeOutliers(neo_q,4,1);
[~,gdp_q_censored]  = removeOutliers(gdp_q,4,1);


% Set monthly and quarterly dates
dates_q = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(6:7))/3,'Q'))';
dates_m = (makeDates(str2num(month_export_start(1:4)),str2num(end_date(1:4)),str2num(month_export_start(end))/3,str2num(end_date(end))/3,'M'))';

%% Get trends and interpolate on monthly data from 1971m4 onwards
[m_rs_trend, m_rs_trend_full]   = fn_local_detrend(rs_q_censored,1,'Retail Sales',bw_bw(1));
[m_ip_trend, m_ip_trend_full]   = fn_local_detrend(ip_q_censored,1,'Industrial Production',bw_bw(2));
[m_pmi_trend, m_pmi_trend_full] = fn_local_detrend(neo_q_censored,1,'New Export Orders',bw_bw(3));
[m_gdp_trend, m_gdp_trend_full] = fn_local_detrend(gdp_q_censored,0,'GDP',bw_bw(4));
if (ind_us==1 | ind_us==2) & foreign==0
    [m_ucl_trend, m_ucl_trend_full] = fn_local_detrend(ucl_q,1,'Initial Claims / Labor Force',bw_bw(5));
end

% Detrend variables
detrended.rs.raw  = rs_data  - m_rs_trend(3:end);
detrended.ip.raw  = ip_data  - m_ip_trend(3:end);
detrended.neo.raw = neo_data - m_pmi_trend(3:end);
detrended.gdp.raw = gdp_data - m_gdp_trend(3:end);
if (ind_us==1 | ind_us==2) & foreign==0
    detrended.ucl.raw = ucl_data - m_ucl_trend(3:end);
end
% Collect trend
trend.rs  = m_rs_trend(3:end);
trend.ip  = m_ip_trend(3:end);
trend.neo = m_pmi_trend(3:end);
trend.gdp.q = m_gdp_trend(3:end);
trend.gdp.full = m_gdp_trend_full(3:end);
if (ind_us==1 | ind_us==2) & foreign==0
    trend.ucl = m_ucl_trend(3:end);
end

% Censore data
[detrended.rs.no_out, detrended.rs.censored]   = removeOutliers(detrended.rs.raw,4,1);
[detrended.ip.no_out, detrended.ip.censored]   = removeOutliers(detrended.ip.raw,4,1);
[detrended.neo.no_out,detrended.neo.censored]  = removeOutliers(detrended.neo.raw,4,1);
[detrended.gdp.no_out,detrended.gdp.censored]  = removeOutliers(detrended.gdp.raw,4,1);

if (ind_us==1 | ind_us==2) & foreign==0
    [detrended.ucl.no_out,detrended.ucl.censored]  = removeOutliers(detrended.ucl.raw,4,1);

    figure; clf;
    plot(detrended.ucl.no_out,'-o'); hold on;
    plot(detrended.ucl.raw);

end
%%
if (ind_us==1 | ind_us==2) & foreign==0
    out=[detrended.rs.censored detrended.ip.censored detrended.neo.censored detrended.gdp.censored detrended.ucl.censored];
end
if (ind_us==0 | ind_us==2) & foreign==1
    out=[detrended.rs.censored detrended.ip.censored detrended.neo.censored detrended.gdp.censored];
end
