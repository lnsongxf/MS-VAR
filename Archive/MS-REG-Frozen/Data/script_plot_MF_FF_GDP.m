clear; clc;

close all;

T = readtable('mGDPMF_US_1973_20200908_uncensored.csv');


DATA = readtable('us_data_for_ads_logs_logs_logs_20200908.csv');


close all;
set(0,'defaulttextinterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

% Adjust Dates
Time4_73_4 = DATA.dates(32:end-7);
Time_73 = datetime(1973,1,1);
Time_20 = datetime(2020,5,1);
ttvec(:,1) = Time_73:Time_20;
Dates_Daily = [year(ttvec), month(ttvec)];
Dates_Monthly = unique(Dates_Daily,'rows');
Time_plot = datetime(Dates_Monthly(:,1),Dates_Monthly(:,2),1, 'Format','MMM-yyyy');

% MAP VARIABLE names
MacroFactor = T.MF;

DataPlot(:,1) = ( DATA.m_ip_detrended(32:end-7) - nanmean(DATA.m_ip_detrended(32:end-7) )) / nanstd( DATA.m_ip_detrended(32:end-7) );
DataPlot(:,2) = ( DATA.m_rs_detrended(32:end-7) - nanmean( DATA.m_rs_detrended(32:end-7) ) ) / nanstd( DATA.m_rs_detrended(32:end-7) );
DataPlot(:,3) = ( DATA.m_neo_detrended(32:end-7) - nanmean(DATA.m_neo_detrended(32:end-7) )) / nanstd( DATA.m_neo_detrended(32:end-7) );
DataPlot(:,4) = -( DATA.m_uclaims_detrended(32:end-7) - nanmean(DATA.m_uclaims_detrended(32:end-7) )) / nanstd( DATA.m_uclaims_detrended(32:end-7) );
DataPlot(:,5) = ( DATA.m_gdp_detrended(32:end-7) - nanmean(DATA.m_gdp_detrended(32:end-7) )) / nanstd( DATA.m_gdp_detrended(32:end-7) );

GDP_norm = (T.dYm-mean(T.dYm(~isnan(T.dYm))))/std(T.dYm(~isnan(T.dYm)));

GDPm_ar = ((((T.dYm./100+1).^(12))-1)*100)+T.trend;

%% Plot common factor and standardized data.u
figure('Units','normalized','Position',[0,0,1,0.33],'Name','MF_data');
h1 = plot(Time_plot,DataPlot(:,1:4),':','LineWidth',1.25); hold on;
h3 = plot(Time_plot,DataPlot(:,5),'o','MarkerSize',3); hold on;
h2 = plot(Time_plot,MacroFactor,'k','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h2 h1(1) h1(2) h1(3) h1(4) h3 ],'Macro factor','Industrial Production','Retail Sales','PMI-NEO','Initial Claims','GDP','Location','southwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

figure('Units','normalized','Position',[0,0,1,0.33],'Name','MF_data_2');
h1 = plot(Time_plot,DataPlot(:,1:4),':','LineWidth',1.25); hold on;
h3 = plot(Time_plot,DataPlot(:,5),'o','MarkerSize',3); hold on;
h2 = plot(Time_plot,MacroFactor,'k','LineWidth',1.5); box on;
h4 = plot(Time_plot,GDP_norm,'b','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h2 h4 h1(1) h1(2) h1(3) h1(4) h3 ],'Macro factor','Montlhy GDP','Industrial Production','Retail Sales','PMI-NEO','Initial Claims','GDP','Location','southwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

figure('Units','normalized','Position',[0,0,1,0.33],'Name','GDP_data');
h3 = plot(Time_plot,(T.gdp_data+T.trend),'o','MarkerSize',3); hold on;
h2 = plot(Time_plot,GDPm_ar,'k','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h2 h3],'Monthly GDP','GDP','Location','southwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

figure('Units','normalized','Position',[0,0,1,0.33],'Name','GDP_data_2');
h2 = plot(Time_plot,T.dYm,'k','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend(h2,'Monthly GDP','Location','southwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

load('FF_data.mat');

figure('Units','normalized','Position',[0,0,1,0.33],'Name','FF_data');
h1 = plot(Time4_73_3,Y_norm,':','LineWidth',1.25); hold on;
h2 = plot(Time4_73_3,Res4_73_3.Z(:,1),'k','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time4_73_3([1 end])); datetick('x','yyyy','keeplimits');
axis tight
recessionplot;lgd = legend([h2 h1(1) h1(2) h1(3) h1(4)],'Financial factor','VXO','EBP','TED spread','CBILL spread','Location','northwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

load('GDPs_data.mat');

figure('Units','normalized','Position',[0,0,1,0.33],'Name','GDP_alt_data');
h1 = plot(Time_plot,GDPs(:,1),'k','LineWidth',1.5); hold on;
h2 = plot(Time_plot,GDPs(:,2),'LineWidth',1.5); hold on;
h3 = plot(Time_plot,GDPs(:,3),'LineWidth',1.5); hold on;
h4 = plot(Time_plot,GDPs(:,4),'LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h1 h2 h3 h4],'Monthly GDP','Weekly Economic Index (Lewis-Mertens-Stock)','IHS Markit MGDP','Stock and Watson Index','Location','southwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

CORRMAT_gdp = NaN(size(GDPs,2),size(GDPs,2));

for row=1:size(GDPs,2)
    for col=1:size(GDPs,2)
        CORRMAT_gdp(row,col) = corr(GDPs(:,row),GDPs(:,col),'row','complete');
    end
end

TABCORR_gdp = array2table(CORRMAT_gdp);
TABCORR_gdp.Properties.VariableNames ={'MonthlyGDP','WEI','Markit','SW'};
TABCORR_gdp.Properties.RowNames ={'MonthlyGDP','WEI','Markit','SW'};


%% Mean and standard deviation
wd = 12;
avg_mf = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));std_mf = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));
avg_ff = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));std_ff = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));
avg_gdp = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));std_gdp = NaN(size(MacroFactor,1)-wd,size(MacroFactor,2));

for rr=wd:size(MacroFactor,1)
    avg_mf(rr-wd+1) = mean(MacroFactor(rr-wd+1:rr,1));
    std_mf(rr-wd+1) = std(MacroFactor(rr-wd+1:rr,1));
    avg_ff(rr-wd+1) = mean(Res4_73_3.Z(rr-wd+1:rr,1));
    std_ff(rr-wd+1) = std(Res4_73_3.Z(rr-wd+1:rr,1));
    avg_gdp(rr-wd+1) = mean(GDPm_ar(rr-wd+1:rr,1));
    std_gdp(rr-wd+1) = std(GDPm_ar(rr-wd+1:rr,1));
end

corr([avg_mf std_mf])
corr([avg_ff std_ff])
corr([avg_gdp std_gdp])
corr([avg_gdp(85:end) std_gdp(85:end)])

figure('Units','normalized','Position',[0,0,1,0.33],'Name','MF_mean_std');
h1 = plot(Time_plot,[NaN(wd-1,1);avg_mf],'k','LineWidth',1.5); hold on;
h2 = plot(Time_plot,[NaN(wd-1,1);std_mf],'b','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h1 h2],'Mean','Std. Deviation','Location','northwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

figure('Units','normalized','Position',[0,0,1,0.33],'Name','GDP_mean_std');
h1 = plot(Time_plot,[NaN(wd-1,1);avg_gdp],'k','LineWidth',1.5); hold on;
h2 = plot(Time_plot,[NaN(wd-1,1);std_gdp],'b','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h1 h2],'Mean','Std. Deviation','Location','northwest');
lgd.FontSize = 16; legend boxoff;
tightfig;

figure('Units','normalized','Position',[0,0,1,0.33],'Name','FF_mean_std');
h1 = plot(Time_plot,[NaN(wd-1,1);avg_ff],'k','LineWidth',1.5); hold on;
h2 = plot(Time_plot,[NaN(wd-1,1);std_ff],'b','LineWidth',1.5); box on;
set(gca,'FontSize',14)
xlim(Time_plot([1 end])); datetick('x','yyyy','keeplimits');
recessionplot;lgd = legend([h1 h2],'Mean','Std. Deviation','Location','northwest');
lgd.FontSize = 16; legend boxoff;
tightfig;





%% Save

h = get(0,'children');
h = sort(h);
for i=1:length(h)
    h = get(0,'children');
    h = sort(h);
    for i=1:length(h)
        set(h(i),'Units','Inches');
        pos = get(h(i),'Position');
        set(h(i),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%         print(h(i),get(h(i),'Name'),'-dpdf','-r0')
        print(h(i),get(h(i),'Name'),'-dpdf','-r0')
    end 
end

%% Correlation Table

DataPlot(:,4) = -DataPlot(:,4);

DataCorr = [MacroFactor DataPlot];

CORRMAT = NaN(size(DataCorr,2),size(DataCorr,2));

for row=1:size(DataCorr,2)
    for col=1:size(DataCorr,2)
        CORRMAT(row,col) = corr(DataCorr(:,row),DataCorr(:,col),'row','complete');
    end
end


TABCORR = array2table(CORRMAT);


TABCORR.Properties.VariableNames ={'MacroFactor','IP','RS','PMI','UCLAIMS','GDP'};
TABCORR.Properties.RowNames ={'MacroFactor','IP','RS','PMI','UCLAIMS','GDP'};