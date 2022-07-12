function [db, db_full,startdb,enddb, tex]= data_exercise(datafilename,sheetuse,start_date,end_date,data_type)

% Preferences
hh      = 12; 
% delta ybar_{t,t+h}: average growth between t and t+h (in months)
detrend = 0;  
% 0 = already detrended (specify trend to be subtracted), 2 = undetrended


%% 1. Data

%check = load('data.mat');

rawdata = readtable(datafilename,'Sheet',sheetuse);

% Map Data
Dates   = rawdata.Dates;

if data_type ==1
    MF_IN   = rawdata.IP;
    GDPG_IN = rawdata.GDPG_US;
    TREND_IN= rawdata.TREND_US;
    FF      = rawdata.EBP;
elseif data_type ==2
    MF_IN   = rawdata.MF_US;
    GDPG_IN = rawdata.GDPG_MARKIT;
    TREND_IN= zeros(size(GDPG_IN));
    FF      = rawdata.FF;
else
    error('Choose between: ''RHS'' or ''LHS'' ')
end

s=find(Dates==datetime(start_date,'InputFormat','yyyy-MMM'));
e=find(Dates==datetime(end_date,'InputFormat','yyyy-MMM'));


%% create variables and put them into the RISE time series format

FF = FF(s:e);
MF = MF_IN(s:e);
GDP_GROWTH = GDPG_IN(s:e);
GDPG_TREND_DFM = TREND_IN(s:e);
y = GDP_GROWTH; 
T = length(y);

% take care of trends
if detrend==0
    y_dfm_trend = GDPG_TREND_DFM;
    Yh_trend = filter(ones(1,hh)/hh, 1, y_dfm_trend);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    y = y-y_dfm_trend;
elseif detrend==2
    Yh_trend = zeros(T,1);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    %y = y;
end

% Data for OLS regression
if hh==0
    yh = y;
else
    Yh = filter(ones(1,hh)/hh, 1, y);
    Yh(1:hh-1) = NaN;
    yh = Yh(hh+1:end);
end

GDP_GROWTH = y;
GDPGH = yh;

d = [GDPGH FF(1:end-hh,:) MF(1:end-hh,:) GDP_GROWTH(1:end-hh,:) y_dfm_trend(1:end-hh,:) yh_trend]; % growth rates
d_full = [[GDPGH;zeros(hh,1)] FF MF GDP_GROWTH y_dfm_trend [yh_trend;zeros(hh,1)]]; % growth rates

vnames={'GDPGH','FF','MF','GDPG','TRENDH','TREND'};

%start='1973M1';
startdb = [num2str(year(datetime(start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(start_date,'InputFormat','yyyy-MMM')))];
enddb = [num2str(year(datetime(Dates(e-12),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e-12),'InputFormat','yyyy-MMM')))];

db=ts(startdb,d,vnames);
db_full=ts(startdb,d_full,vnames);

db=pages2struct(db);
db_full=pages2struct(db_full);

tex=struct();
if data_type ==1
    tex.GDPGH = 'Average Future GDP Growth';
    tex.FF    = 'Excess Bond Premium';
    tex.MF    = 'Industrial Production';
    tex.GDPG  = 'GDP Growth';
    tex.TRENDH= 'Average Future GDP Growth Trend';
    tex.TREND = 'GDP Growth Trend';
elseif data_type ==2
    tex.GDPGH = 'Average Future GDP Growth (Markit)';
    tex.FF    = 'Financial Factor';
    tex.MF    = 'Macroeconomic Factor';
    tex.GDPG  = 'GDP Growth (Markit)';
    tex.TRENDH= 'Average Future GDP Growth Trend';
    tex.TREND = 'GDP Growth Trend';
end



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
