% Function to transform data to RISE db format
% and construct average growth series:
% delta ybar_{t,t+h}: average growth between t and t+h (in months)
%==================================================================================
function [db, db_full,dates, startdb,enddb,tex]= fLoadData(opt)

%% LOAD Data
rawdata = readtable(['data/' opt.datafilename '.xlsx'],'Sheet',opt.sheetuse,'ReadVariableNames',true);

% Map dates
inputformat = 'yyyy-MMM';
dates = datetime(rawdata.Dates,'Format',inputformat);

% Get data start and end dates
s=find(dates==opt.start_date);
e=find(dates==opt.end_date);
ep = find(dates==opt.end_plot);

% Map from options
hh = opt.hh;

%% Map variables

% Variable names
varlist = opt.varlist; 

% Map variables and simulation sample
for jj=1:length(varlist)
    data_use.(varlist{jj}) = rawdata.(varlist{jj})(s:ep);
    data_use.([varlist{jj} '_TREND']) = rawdata.([varlist{jj} '_TREND'])(s:ep);
end


for jj=1:length(varlist)
    
    % Collect trends
    y_trend.(varlist{jj}) = data_use.([varlist{jj} '_TREND']);
    
    % Detrend series
    y.(varlist{jj}) = data_use.(varlist{jj})-y_trend.(varlist{jj});
    
    % Compute hh-period ahead average of trends
    Yh_trend = filter(ones(1,hh)/hh, 1, y_trend.(varlist{jj}));
    Yh_trend(1:hh-1) = NaN;
    yh_trend.(varlist{jj}) = [Yh_trend(hh:ep); repmat(Yh_trend(ep),hh-1,1)];
    
    % Compute hh-period ahead average
    if hh==0
        % Detrended series
        yh.(varlist{jj}) = y.(varlist{jj});
        
        % Series with trend
        yh_wt.(varlist{jj}) = data_use.(varlist{jj});
    else
        % Detrended series
        Yh = filter(ones(1,hh)/hh, 1, y.(varlist{jj}));
        Yh(1:hh-1) = NaN;
        yh.(varlist{jj}) = [Yh(hh:ep); zeros(hh-1,1)];
        
        % Series with trend
        Yh_wt = filter(ones(1,hh)/hh, 1, data_use.(varlist{jj}));
        Yh_wt(1:hh-1) = NaN;
        yh_wt.(varlist{jj}) = [Yh_wt(hh:ep); zeros(hh-1,1)];

    end
    
end


%Changed to the following:

if opt.dir_or_it ==1
    
    % Estimation sample
    for jj=1:length(varlist)
        d.(varlist{jj}) = y.(varlist{jj})(1:e-hh+1,:);
        d.([varlist{jj} 'H']) = yh.(varlist{jj})(1:e-hh+1,:);
        d.(['TREND_' varlist{jj}]) = y_trend.(varlist{jj})(1:e-hh+1,:);
        d.(['TRENDH_' varlist{jj} 'H'])= yh_trend.(varlist{jj})(1:e-hh+1,:);
        d.([varlist{jj} '_wt']) = data_use.(varlist{jj})(1:e-hh+1,:);
        d.([varlist{jj} 'H_wt']) = yh_wt.(varlist{jj})(1:e-hh+1,:);

    end
    
    % Full sample
    for jj=1:length(varlist)
        d_full.(varlist{jj}) = y.(varlist{jj});
        d_full.([varlist{jj} 'H']) = yh.(varlist{jj});
        d_full.(['TREND_' varlist{jj} ]) = y_trend.(varlist{jj});
        d_full.(['TRENDH_' varlist{jj} 'H'])= yh_trend.(varlist{jj});
        d_full.([varlist{jj} '_wt']) = data_use.(varlist{jj});
        d_full.([varlist{jj} 'H_wt']) = yh_wt.(varlist{jj});

    end
    
elseif opt.dir_or_it==2
    
    
    % Estimation sample
    for jj=1:length(varlist)
        d.(varlist{jj}) = y.(varlist{jj})(1:e,:);
        d.([varlist{jj} 'H']) = yh.(varlist{jj})(1:e,:);
        d.(['TREND_' varlist{jj} ]) = y_trend.(varlist{jj})(1:e,:);
        d.(['TRENDH_' varlist{jj}  'H'])= yh_trend.(varlist{jj})(1:e,:);
        d.([varlist{jj} '_wt']) = data_use.(varlist{jj})(1:e,:);
        d.([varlist{jj} 'H_wt']) = yh_wt.(varlist{jj})(1:e,:);

    end
    
    
    % Full sample
    for jj=1:length(varlist)
        d_full.(varlist{jj}) = y.(varlist{jj});
        d_full.([varlist{jj} 'H']) = yh.(varlist{jj});
        d_full.(['TREND_' varlist{jj} ]) = y_trend.(varlist{jj});
        d_full.(['TRENDH_' varlist{jj}  'H'])= yh_trend.(varlist{jj});
        d_full.([varlist{jj} '_wt']) = data_use.(varlist{jj});
        d_full.([varlist{jj} 'H_wt']) = yh_wt.(varlist{jj});        
    end
    
end



%=========================
% CREATE RISE DATABASE
%=========================

% Start and end dates for database
% ***Truncate estimation samples hh-months prior to enddate***
% added (hh+1) to adjust for current periot (t)
%


% Variable Names
vnames = fieldnames(d);


startdb = char(datetime(dates(s),'Format','yyyy''M''M'));
if opt.dir_or_it ==1
    enddb   = char(datetime(dates(e-hh+1),'Format','yyyy''M''M'));
elseif opt.dir_or_it ==2
    enddb   = char(datetime(dates(e),'Format','yyyy''M''M'));
end

% create databases
db=ts(startdb,struct2array(d),vnames);
db_full=ts(startdb,struct2array(d_full),vnames);

% transform to structures
db=pages2struct(db);
db_full=pages2struct(db_full);

tex=struct();
for jj=1:length(varlist)
    tex.(varlist{jj}) = opt.vlabels{jj};
end
