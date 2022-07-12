% Function to transform data to RISE db format
% and construct average growth series:
% delta ybar_{t,t+h}: average growth between t and t+h (in months)
%==================================================================================
function [db, db_full,startdb,enddb, tex]= fLoadDataSim(rawdata,SEED,start_date,end_date,hh)

% Preferences
detrend = 0;   % 0 = already detrended (specify trend to be subtracted), 2 = undetrended

%% 1. Data

% Map Data
Dates    = rawdata.Dates;
FF_IN    = rawdata.FF(SEED,:)';
MF_IN    = rawdata.MF(SEED,:)';
GDPG_IN  = rawdata.GDPG(SEED,:)';
TREND_IN = zeros(size(GDPG_IN));

% Get data start and end dates
s=find(Dates==datetime(start_date,'InputFormat','yyyy-MMM'));
e=find(Dates==datetime(end_date,'InputFormat','yyyy-MMM'));


%% create variables and put them into the RISE time series format
GDP_GROWTH = GDPG_IN(s:e);
GDPG_TREND = TREND_IN(s:e);
FF = FF_IN(s:e);
MF = MF_IN(s:e);
y  = GDP_GROWTH; 
T  = length(y);

% take care of trends
if detrend==0
    y_dfm_trend = GDPG_TREND;
    Yh_trend = filter(ones(1,hh)/hh, 1, y_dfm_trend);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    y = y-y_dfm_trend;
elseif detrend==2
    % *** We don't really need this as input file should specify the trend (PC 02/22/21)***
    Yh_trend = zeros(T,1);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    %y = y;
end

% Data for estimation
if hh==0
    yh = y;
else
    % Compute future GDP growth
    Yh = filter(ones(1,hh)/hh, 1, y);
    Yh(1:hh-1) = NaN;
    yh = Yh(hh+1:end);
end

GDP_GROWTH = y;     % Current GDP
GDPGH = yh;         % Average future GDP at horizon H

% Estimation sample
%d = [GDPGH FF(1:end-hh,:) MF(1:end-hh,:) GDP_GROWTH(1:end-hh,:) y_dfm_trend(1:end-hh,:) yh_trend]; % growth rates

% Truncate estimation samples 12-months prior to enddate
d = [GDPGH(1:(e-s+1)-12,:) FF(1:(e-s+1)-12,:) MF(1:(e-s+1)-12,:) GDP_GROWTH(1:(e-s+1)-12,:) y_dfm_trend(1:(e-s+1)-12,:) yh_trend(1:(e-s+1)-12,:)]; % growth rates

% Full sample
d_full = [[GDPGH;zeros(hh,1)] FF MF GDP_GROWTH y_dfm_trend [yh_trend;zeros(hh,1)]]; % growth rates

% Variable names
% *** This was incorrect ***
%vnames={'GDPGH','FF','MF','GDPG','TRENDH','TREND'};
% Correct should be:
vnames={'GDPGH','FF','MF','GDPG','TREND','TRENDH'};
%=========================
%CREATE RISE DATABASE 
%=========================

% Start and end dates for database
% ***Truncate estimation samples 12-months prior to enddate***
startdb = [num2str(year(datetime(start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(start_date,'InputFormat','yyyy-MMM')))];
enddb   = [num2str(year(datetime(Dates(e-hh),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e-hh),'InputFormat','yyyy-MMM')))];

% create databases
db=ts(startdb,d,vnames);
db_full=ts(startdb,d_full,vnames);

% transform to structures
db=pages2struct(db);
db_full=pages2struct(db_full);

tex=struct();
tex.GDPGH = 'Average Future GDP Growth';
tex.FF    = 'Financial Factor';
tex.MF    = 'Macroeconomic Factor';
tex.GDPG  = 'GDP Growth';
tex.TRENDH= 'Average Future GDP Growth Trend';
tex.TREND = 'GDP Growth Trend';


%% plot the data 
% mynames=fieldnames(db);
% 
% figure('name','database');
% for ii=1:numel(mynames)
%     v=mynames{ii};
%     subplot(6,1,ii)
%     plot(db.(v))
%     title(tex.(v))
% end
% xrotate(45)

%% save data for later use
%save('data','db','tex')
%save('data_full','db_full','tex')
