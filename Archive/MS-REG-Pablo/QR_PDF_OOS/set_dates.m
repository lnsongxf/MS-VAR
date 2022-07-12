%% Set Dates for Sample
if sample(sam)==1
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-Apr';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2020-Q1';
    end
elseif sample(sam)==2
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-Apr';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q1';
    end
elseif sample(sam)==3
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2007-Dec';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2007-Q4';        
    end
elseif sample(sam)==4
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2007-Dec';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2007-Q4';                
    end
elseif sample(sam)==5
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-Mar';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2020-Q1';                
    end
elseif sample(sam)==6
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-Mar';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q1';                
    end
elseif sample(sam)==7
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-May';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2020-Q1';                
    end
elseif sample(sam)==8
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-May';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q1';                
    end
elseif sample(sam)==9
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-Jun';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2020-Q2';                
    end
elseif sample(sam)==10
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-Jun';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q2';                
    end
elseif sample(sam)==11
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-Jul';
    elseif frequency==4
        startdate = '1973-Q1';
        enddate   = '2020-Q2';                
    end
elseif sample(sam)==12
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-Jul';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q2';                
    end
elseif sample(sam)==13
    if frequency==12
        startdate = '1973-Jan';
        enddate   = '2020-Aug';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q2';                
    end
elseif sample(sam)==14
    if frequency==12
        startdate = '1986-Jan';
        enddate   = '2020-Aug';
    elseif frequency==4
        startdate = '1986-Q1';
        enddate   = '2020-Q2';                
    end
end

% Collect date vectors
dates_dfm   = data_dfm_org.Dates; 
dates_other = data_other_org.Dates;

% Take care of specifications for which data starts later
if specs(ss)==2
    if sample(sam)==1 || sample(sam)==3 || sample(sam)==5 || sample(sam)==7
        if frequency==12
            startdate = '1980-Jan';  
        else
            startdate = '1980-Q1';  
        end
    end
elseif specs(ss)==3
    if frequency==12
        startdate = '1990-Feb'; 
    elseif frequency ==4
        startdate = '1990-Q2'; 
    end
elseif specs(ss)==13
    if sample(sam)<3
        if frequency==12
            enddate = '2020-Mar';
        elseif frequency==4
            enddate = '2020-Q1';
        end
    end
elseif specs(ss)==17
    frequency=4;
    startdate = '1973-Q1';
    enddate   = '2015-Q4';
end

% Set start and end dates
if frequency==12
    inputformat = 'yyyy-MMM';
    dataformat  = 'yyyy-mmm';
    funit = 'Months';
elseif frequency==4
    inputformat = 'yyyy-QQ';
    dataformat  = 'yyyy-qq';
    funit = 'Quarters';
end
start_num       = find(datetime(startdate,'InputFormat',inputformat)==dates_dfm);
end_num         = find(datetime(enddate,'InputFormat',inputformat)==dates_dfm);
start_num_other = find(datetime(startdate,'InputFormat',inputformat)==dates_other);
end_num_other   = find(datetime(enddate,'InputFormat',inputformat)==dates_other);

% Cut sample
data_dfm   = data_dfm_org(start_num:end_num,:);
data_other = data_other_org(start_num_other:end_num_other,:);

% Create date vector
if frequency==12
    dates = datenum((datetime(startdate,'InputFormat',inputformat)+calmonths(p)):calmonths(1):(datetime(enddate,'InputFormat',inputformat)))';
else
    dates = datenum((datetime(startdate,'InputFormat',inputformat)+calquarters(p)):calquarters(1):(datetime(enddate,'InputFormat',inputformat)))';    
end
