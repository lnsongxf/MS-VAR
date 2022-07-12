function [db,varlist]=createData

    %% Create Data
    %==========================================================================
    % variables are in the following order: the description corresponds to the updated dataset
    % 1: Real GDP               (Billions of Chained 2009 Dollars Seasonally Adjusted Annual Rate)
    % 2: GDP Deflator           (Index 2009=100)
    % 3: PCOM                   (Dow Jones Spot Average (Symbol _DJSD) Commodity Price Index, commercially available from Global Financial Data Index 1982-1984=100)
    % 4: Total Reserves         (Billions of Dollars)
    % 5: Non-borrowed Reserves  (Billions of Dollars)
    % 6: Federal Funds Rate     (Annual Rate in Percentages)
    data  = readtable('dataset.csv');
    data   = data(find((data.dates=='1965-01-01')==1):find((data.dates=='2007-06-01')==1),:);
    
    % all variables are in log times 100 except for the federal funds rate that enters the SVAR in annualized percentages
    d  = [100*log(data.monthly_GDP) 100*log(data.monthly_GDPDEF) 100*log(data.CPRINDEX) 100*log(data.TRARR) 100*log(data.BOGNONBR) data.FEDFUNDS];

    % names of the variables
    varlist=struct();
    varlist.GDP='Real GDP';
    varlist.PGDP='GDP Deflator';
    varlist.PCOM='Commodity Prices';
    varlist.TR='Total Reserves';
    varlist.NBR='Non-borrowed Reserves';
    varlist.FFR='Federal Funds Rate';

    % assign data
    vardata.GDP  = d(:,1);
    vardata.PGDP = d(:,2);
    vardata.PCOM = d(:,3);
    vardata.TR   = d(:,4);
    vardata.NBR  = d(:,5);
    vardata.FFR  = d(:,6);

    % create databases
    start_date=sprintf('%0.0dM%0.0d',1965,1);
    db=ts(start_date,d,fieldnames(varlist));

    % transform to structures
    db=pages2struct(db);

    % plot data
    fh=figure('name','Variables in the VAR');
    fields=fieldnames(varlist);
    for iplot=1:length(fields)
        vname=fields{iplot};
        subplot(3,2,iplot)
        plot(db.(vname),'linewidth',2)
        mytitle=strrep([varlist.(vname),'(',vname,')'],'"','');
        title(mytitle)
    end
    
end