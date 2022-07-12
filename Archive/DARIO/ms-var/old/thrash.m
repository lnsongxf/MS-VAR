close all;
% Index of dates

% Data vintage, sample and country selection
start_date   = '2100-Jan';
end_date     = '2199-Dec';

% Define dates for plots
% Data vintage, sample and country selection
start_plot     = '2100-Feb';
end_plot       = '2198-Dec';
end_plot_full  = '2199-Dec';

% Date formats
inputformat   = 'yyyy-MMM';
dataformat    = 'yyyy-mmm';

sheetuse = 'Sim';

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

% Figure options
FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);
left_color = [0 0 0];
right_color = colors(2,:);

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);


load('results/Params/mode/Simul/Iterated/quantiles_iterated_test.mat');
load('results/Params/mode/Simul/Iterated/quantiles_iterated.mat');


fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(sd:ed), quantiles_iterated.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Alg. 1)');
l2=plot(dates_full(sd:ed), quantiles_iterated.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Alg. 1)');
l3=plot(dates_full(sd:ed), quantiles_iterated_test.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Alg. 2)');
l4=plot(dates_full(sd:ed), quantiles_iterated_test.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Alg. 2)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Quantiles Iterated','Interpreter','Latex','FontSize',16);

axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed)])


load('results/Params/mode/Simul/Direct/quantiles_direct_test.mat');
load('results/Params/mode/Simul/Direct/quantiles_direct.mat');

fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(sd:ed), quantiles_direct.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Alg. 1)');
l2=plot(dates_full(sd:ed), quantiles_direct.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Alg. 1)');
l3=plot(dates_full(sd:ed), quantiles_direct_test.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Alg. 2)');
l4=plot(dates_full(sd:ed), quantiles_direct_test.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Alg. 2)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Quantiles Direct','Interpreter','Latex','FontSize',16);

axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed)])

fig3=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(sd:ed), quantiles_iterated.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Iterated)');
l2=plot(dates_full(sd:ed), quantiles_iterated.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Iterated)');
l3=plot(dates_full(sd:ed), quantiles_direct.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Direct)');
l4=plot(dates_full(sd:ed), quantiles_direct.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Direct)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Quantiles Iterated (Alg.1)','Interpreter','Latex','FontSize',16);

axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed)])

fig4=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(sd:ed), quantiles_iterated_test.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Iterated)');
l2=plot(dates_full(sd:ed), quantiles_iterated_test.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Iterated)');
l3=plot(dates_full(sd:ed), quantiles_direct_test.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Direct)');
l4=plot(dates_full(sd:ed), quantiles_direct_test.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Direct)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Quantiles Iterated (Alg.2)','Interpreter','Latex','FontSize',16);

axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed)])


% 
% 
% 
% 
% 
% if exist('results/Params/mode/Simul/Direct/quantiles_direct_test.mat','file')==2 && exist('results/Params/mode/Simul/Iterated/quantiles_iterated_test.mat','file')==2
%     load('results/Params/mode/Simul/Direct/quantiles_direct_test.mat');
%     load('results/Params/mode/Simul/Iterated/quantiles_iterated_test.mat');
%     
%     
%     fig6=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
%     hold on
%     l1=plot(dates_full(sd:ed), quantiles_direct.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (direct)');
%     l2=plot(dates_full(sd:ed), quantiles_direct.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (direct)');
%     l3=plot(dates_full(sd:ed), quantiles_iterated.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (iterated)');
%     l4=plot(dates_full(sd:ed), quantiles_iterated.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (iterated)');
%     legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
%     set(gca,'children',flipud(get(gca,'children')))
%     hold off
%     ylabel('Percent','interpreter','Latex','fontsize',10)
%     title('Quantiles Direct test and Iterated test','Interpreter','Latex','FontSize',16);
%     
%     axis tight
%     datetick('x','yyyy','keepticks')
%     set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
%     fprintf('\n Correlation 10th quantile: %g',round(corr(quantiles_direct.dYsim_10,quantiles_iterated.dYsim_10),2));
%     fprintf('\n Correlation 90th quantile: %g \n',round(corr(quantiles_direct.dYsim_90,quantiles_iterated.dYsim_90),2));
% end
