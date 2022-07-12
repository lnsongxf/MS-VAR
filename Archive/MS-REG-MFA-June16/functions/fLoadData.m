% Function to transform data to RISE db format
% and construct average growth series:
% delta ybar_{t,t+h}: average growth between t and t+h (in months)
%==================================================================================
function [db, db_full,startdb,enddb, tex]= fLoadData(opt)


% Map from options
hh = opt.hh; 

%% 1. Data

rawdata = readtable(['data/' opt.datafilename '.xlsx'],'Sheet',opt.sheetuse,'ReadVariableNames',true);

% Map Data
Dates    = rawdata.Dates;
PCEPILFE_CCA_IN = rawdata.PCEPILFE_CCA;
FEDFUNDS_IN = rawdata.FEDFUNDS;
FF_IN    = rawdata.FF;
MF_IN    = rawdata.MF;
GDPG_IN  = rawdata.GDPG;
TREND_IN = rawdata.TREND;

% Get data start and end dates
s=find(Dates==datetime(opt.start_date,'InputFormat','yyyy-MMM'));
e=find(Dates==datetime(opt.end_date,'InputFormat','yyyy-MMM'));
e_sim = find(Dates==datetime(opt.end_plot,'InputFormat','yyyy-MMM'));

%% create variables and put them into the RISE time series format
GDP_GROWTH = GDPG_IN(s:e_sim);
GDPG_TREND = TREND_IN(s:e_sim);
FF = FF_IN(s:e_sim);
MF = MF_IN(s:e_sim);
FEDFUNDS = FEDFUNDS_IN(s:e_sim);
PCEPILFE_CCA = PCEPILFE_CCA_IN(s:e_sim);
y  = GDP_GROWTH; 
T  = length(y);

% take care of trends
y_dfm_trend = GDPG_TREND;
Yh_trend = filter(ones(1,hh)/hh, 1, y_dfm_trend);
Yh_trend(1:hh-1) = NaN;
yh_trend = Yh_trend(hh:e_sim); 
y = y-y_dfm_trend;

% Data for estimation
if hh==0
    yh = y;
else
    % Compute future GDP growth
    Yh = filter(ones(1,hh)/hh, 1, y);
    Yh(1:hh-1) = NaN;
    yh = Yh(hh:e_sim);
end

GDP_GROWTH = y;     % Current GDP
GDPGH = yh;         % Average future GDP at horizon H

%Changed to the following:
if opt.dir_or_it ==1
    
    % Estimation sample
    d = [GDPGH(1:e-hh+1,:) PCEPILFE_CCA(1:e-hh+1,:)  FEDFUNDS(1:e-hh+1,:) FF(1:e-hh+1,:) MF(1:e-hh+1,:) GDP_GROWTH(1:e-hh+1,:) y_dfm_trend(1:e-hh+1,:) yh_trend(1:e-hh+1,:)]; % growth rates
    % Full sample
    d_full = [[GDPGH;zeros(hh-1,1)] PCEPILFE_CCA  FEDFUNDS FF MF GDP_GROWTH y_dfm_trend [yh_trend;zeros(hh-1,1)]]; % growth rates
    
elseif opt.dir_or_it==2
    
    % Estimation sample
    d = [[GDPGH(1:e-hh+1,:);zeros(hh-1,1)] PCEPILFE_CCA(1:e,:)  FEDFUNDS(1:e,:)  FF(1:e,:) MF(1:e,:) GDP_GROWTH(1:e,:) y_dfm_trend(1:e,:) [yh_trend(1:e-hh+1,:);zeros(hh-1,1)]]; % growth rates    
    % Full sample
    d_full = [[GDPGH;zeros(hh-1,1)] PCEPILFE_CCA  FEDFUNDS FF MF GDP_GROWTH y_dfm_trend [yh_trend;zeros(hh-1,1)]]; % growth rates

end

% Variable names
vnames={'GDPGH','PCEPILFE_CCA','FEDFUNDS','FF','MF','GDPG','TREND','TRENDH'};

%=========================
%CREATE RISE DATABASE 
%=========================

% Start and end dates for database
% ***Truncate estimation samples hh-months prior to enddate***
% added (hh+1) to adjust for current periot (t)
%
startdb = [num2str(year(datetime(opt.start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(opt.start_date,'InputFormat','yyyy-MMM')))];
if opt.dir_or_it ==1
    enddb   = [num2str(year(datetime(Dates(e-hh+1),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e-hh+1),'InputFormat','yyyy-MMM')))];
elseif opt.dir_or_it ==2
    enddb   = [num2str(year(datetime(Dates(e),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e),'InputFormat','yyyy-MMM')))];
end

% create databases
db=ts(startdb,d,vnames);
db_full=ts(startdb,d_full,vnames);

% transform to structures
db=pages2struct(db);
db_full=pages2struct(db_full);

tex=struct();
tex.GDPGH = 'Average Future Real GDP Growth';
tex.FEDFUNDS    = 'Federal Funds Rate';
tex.PCEPILFE_CCA    = 'Core PCE Inflation (CCA)';
tex.FF    = 'Financial Factor';
tex.MF    = 'Macroeconomic Factor';
tex.GDPG  = 'Real GDP Growth';
tex.TRENDH= 'Average Future Real GDP Growth Trend';
tex.TREND = 'Real GDP Growth Trend';
