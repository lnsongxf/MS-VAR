  load Data_CreditDefaults   
  X0 = Data(:,1:4);
  T0 = size(X0,1);
  
  % Convert dates to serial date numbers:

  dates = datenum([dates,ones(T0,2)]);

  % Create plot:
  
  fig = figure;
  ax = axes(fig);
  plot(ax,dates,X0,'LineWidth',2)
  set(ax,'XTick',dates(1:2:end))
  datetick(ax,'x','yyyy','keepticks')
  xlabel(ax,'Year') 
  ylabel(ax,'Level')
  axis(ax,'tight')
 
  % Add recession bands:

  rr=recessionplot;
%% foreign recessions
inputformat = 'yyyy-MMM';
start_date   = '1980-Jan';
end_date     = '2020-Dec';
nlags = 1;
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% 1974m9	1
% 1974m12	1

% 1982m3	1
% 1982m11	1

% 1992m2	1
% 1992m3	1

% 2001m4	1
% 2002m12	1

% 2003m1	1
% 2003m4	1

% 2008m4	1
% 2009m5	1

%%
  fig = figure;
  ax = axes(fig);
  plot(ax,dates_full,randn(491,1),'LineWidth',2)
  set(ax,'XTick',dates_full(1:120:end))
  datetick(ax,'x','yyyy-mmm','keepticks')
  xlabel(ax,'Year') 
  ylabel(ax,'Level')
  axis(ax,'tight')

recdates(1,:) = [datenum(datetime('1974-Sep','InputFormat',inputformat)) ...
    datenum(datetime('1974-Dec','InputFormat',inputformat))];

recdates(2,:) = [datenum(datetime('1982-Mar','InputFormat',inputformat)) ...
    datenum(datetime('1982-Nov','InputFormat',inputformat))];

recdates(3,:) = [datenum(datetime('1992-Feb','InputFormat',inputformat)) ...
    datenum(datetime('1992-Mar','InputFormat',inputformat))];

recdates(4,:) = [datenum(datetime('2001-Apr','InputFormat',inputformat)) ...
    datenum(datetime('2001-Dec','InputFormat',inputformat))];

recdates(5,:) = [datenum(datetime('2002-Dec','InputFormat',inputformat)) ...
    datenum(datetime('2003-Apr','InputFormat',inputformat))];

recdates(6,:) = [datenum(datetime('2008-Apr','InputFormat',inputformat)) ...
    datenum(datetime('2009-May','InputFormat',inputformat))];


rr = recessionplot('recessions',recdates)