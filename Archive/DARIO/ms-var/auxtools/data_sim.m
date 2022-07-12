function [db, db_full, tex]= data_sim(Res,draw,start_date)

% Preferences
hh      = 12; 
% delta ybar_{t+1,t+h+1}: average growth between t and t+h (in months)

%% create variables and put them into the RISE time series format

FF = Res.Res_iterated.FF(draw,:)';
MF = Res.Res_iterated.MF(draw,:)';
GDP_GROWTH = Res.Res_iterated.y_mat(draw,:)';
GDPGH = Res.Res_direct.y_mat(draw,:)';

d = [GDPGH FF(1:end-hh,:) MF(1:end-hh,:) GDP_GROWTH(1:end-hh,:)]; % growth rates
d_full = [[GDPGH;zeros(hh,1)] FF MF GDP_GROWTH]; % growth rates

vnames={'GDPGH','FF','MF','GDPG'};

%start='1973M1';
startdb = [num2str(year(datetime(start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(start_date,'InputFormat','yyyy-MMM')))];

db=ts(startdb,d,vnames);
db_full=ts(startdb,d_full,vnames);

db=pages2struct(db);
db_full=pages2struct(db_full);

tex=struct();
tex.GDPGH = 'Average Future GDP Growth';
tex.FF    = 'Financial Factor';
tex.MF    = 'Macroeconomic Factor';
tex.GDPG  = 'GDP Growth';


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
