clear; clc; close all

%% 1. Data

load('CHdata.mat');

Dates = datenum(time);
str_sample_init = '1990-01-01';         % Starting date of the sample (include pre-sample)
str_iv_init     = '1994-01-01';         % Starting date of the sample for the proxy
start='1994M1';                         % Declare starting date to RISE
str_sample_end  = '2007-06-01';         % End date of the sample
s=find(time==datetime(str_iv_init,'InputFormat','yyyy-MM-dd'));
e=find(time==datetime(str_sample_end,'InputFormat','yyyy-MM-dd'));


%% create variables and put them into the RISE time series format


d = YYdata(s:e,:);

db=ts(start,d,vnames);
db=pages2struct(db);

tex=struct();
tex.EFFR = 'Federal funds rate';
tex.LIPM = 'Industrial production';
tex.UNRATE = 'Unemployment rate';
tex.LPPI = 'PPI';
tex.BAA10YMOODY = 'BAA Spread';
tex.MHF = 'HFI MP Shock';


%% plot the data 
mynames=fieldnames(db);

figure('name','database');
for ii=1:numel(mynames)
    v=mynames{ii};
    subplot(6,1,ii)
    plot(db.(v))
    title(tex.(v))
end
xrotate(45)

%% save data for later use
save('data','db','tex')

