% Script to produce figures for slides
% First run msreg_2reg_endoprob.m

close all; 


%%

%==========================================================================
% Quantiles from MS Model: .10, .25, .75, .90
%==========================================================================

% Compute percentiles
dY_sim_25 = prctile(dY_sim',25)';
dY_sim_75 = prctile(dY_sim',75)';
dY_sim_10 = prctile(dY_sim',10)';
dY_sim_90 = prctile(dY_sim',90)';

% Compute means
dY_sim1_mean = mean(dY_sim_1,2);
dY_sim2_mean = mean(dY_sim_2,2);

% Set colors
colors = cbrewer('div', 'RdYlBu', 64);
numticks = 48;
fig=figure; clf;
hold on
% Plot percentiles from MS
l1=plot(dates(sd:ed), dY_sim_10(sd:ed),'Color',colors(5,:),'LineWidth', 3,'DisplayName','10th');
l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'--','Color',colors(15,:),'LineWidth', 2.5,'DisplayName','25th');
l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'--','Color',colors(55,:),'LineWidth', 2.5,'DisplayName','75th');
l4=plot(dates(sd:ed), dY_sim_90(sd:ed),'Color',colors(60,:),'LineWidth', 3,'DisplayName','90th');
rr=recessionplot;
hleg = legend([l1 l2 l3 l4],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('$\bar{\Delta} {y}_{t+1,t+12}$ (\%)','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig,[slides_folder model '_Quantiles_MS'],'-bestfit');
%%
%==========================================================================
% Historical Quantiles QR and MS
%==========================================================================

fig=figure; clf;
%yyaxis left
hold on
% Plot percentiles from MS
l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'--','Color',colors(15,:),'LineWidth', 3,'DisplayName','MS 25th');
l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'--','Color',colors(55,:),'LineWidth', 3,'DisplayName','MS 75th');

% plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
l5=plot(dates(sd:ed),YQ_dist_fut(sd:ed,4),'-','Color',colors(20,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(1)*100) 'th ']);
l6=plot(dates(sd:ed),YQ_dist_fut(sd:ed,14),'-','Color',colors(50,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(3)*100) 'th ']);
ylim([-8 8]);
rr=recessionplot;
hleg = legend([l2 l3 l5 l6],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))

hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig,[slides_folder model '_Quantiles_MS_QR'],'-bestfit');


%%
%==========================================================================
% Density comparison at period1=October-2008 and December-2018
%==========================================================================

md = '2008-Oct';
mydate = find(datenum(md,dataformat)==dates);
od = '2018-Dec';
mydate2 = find(datenum(od,dataformat)==dates);

mcolor = colors(5,:);
ocolor = colors(60,:);

fig=figure;
% At ergodic distribution (Pablo)
[pdf3,xi3]=ksdensity(y_erg_bar(mydate,:));
[pdf4,xi4]=ksdensity(y_erg_bar(mydate2,:));

hold on
plot(xi3,pdf3,'-','Color',mcolor,'LineWidth', 3,'DisplayName',[md ', P(Bad regime) = ' num2str(round(p_reg2_filtered(mydate),2)) ]);
plot(xi4,pdf4,'-.','Color',ocolor,'LineWidth', 3,'DisplayName',[od ', P(Bad regime) = '  num2str(round(p_reg2_filtered(mydate2),2))]);
plot([realized(mydate) realized(mydate)],[0 max([pdf3,pdf4])],'-','Color',mcolor,'LineWidth',1,'DisplayName',[md ', Realized'])
plot([realized(mydate2) realized(mydate2)],[0 max([pdf3,pdf4])],'-','Color',ocolor,'LineWidth',1,'DisplayName',[md ', Realized'])


%plot([0 0],[0 max([pdf,pdf2])],'k-.','LineWidth',2,'HandleVisibility','off')
hold off
legend(gca,'Location','NorthWest','interpreter','Latex')
legend boxoff
ylabel('PDF','fontsize',10,'interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-10 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%%
% Save file
%print('-dpdf',fig,[slides_folder model '_KernelDensity_2008_2018'],'-bestfit');

%==========================================================================
% Risk Assessment GFC
%==========================================================================

period1 = '2007-Aug';
period2 = '2007-Sep';
period3 = '2008-Jan';
period4 = '2008-Mar';
staffforecast = [2.2, 1.4, 1.5, 0.7];
title_use = {'August 2007','September 2007','January 2008','March 2008'};
colors = cbrewer('div', 'RdBu', 8);

% Colors for densities
mcolor = colors(8,:); % Color for QR
ocolor = colors(1,:); % Color for MS

% Extra color palette
colors2 = cbrewer('div', 'RdYlGn', 8);

dcolor = [mcolor;ocolor];

if exist('ResMatch')==0
    load('results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

% Select Time Periods for Which to Plot PDFs
period_pdf = {period1, period2, period3, period4};

PST = ResMatch.PST;

fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -8;
xmax = 8;

for t=1:length(period_pdf)

    tind(t) = find(datenum(char(period_pdf(t)),dataformat)==dates);
    
    % rescale  pdf so that the range is between 0 and 1
    CDF = ResMatch.CST(tind(t),:);
   

    % plot figure for in-sample pdfs
    %if t==1
    %    ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    %end
    subplot(2,2,t)
    % PDF Plot
    hold on;
    l1=plot(ResMatch.YY,PST(tind(t),:)','-','Color',dcolor(1,:),'LineWidth',2,'DisplayName','QR Predictive Density');
    [pdf,xi]=ksdensity(y_erg_bar(tind(t),:));
    l2=plot(xi,pdf,'-','LineWidth', 2,'Color',dcolor(2,:),'DisplayName','MS Predictive Density');
    l3 = plot([realized(tind(t)) realized(tind(t))],[0 yl(2)],'k-','LineWidth',1,'DisplayName',['Realization = ' num2str(round(realized(tind(t)),1)) ' \%']);
    l4 = plot([staffforecast(t) staffforecast(t)],[0 yl(2)],'-.','Color',colors2(8,:),'LineWidth',0.8,'DisplayName',['Staff Forecast = ' num2str(round(staffforecast(t),1)) ' \%']);
    hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','NorthWest');
    set(hleg,'interpreter','Latex'); legend boxoff; 
    ylim([0 0.6]); xlim([-12 6]); ylabel('PDF','Interpreter','Latex');
    title(title_use(t),'Interpreter','Latex');

end
%axis tight
yl = ylim;
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig,[slides_folder model '_RiskAssessment_GFC'],'-bestfit');

%%
% Transition probabilities when coefficients have gamma prior
FF_grid = [-6:0.1:6]'; %[min(FF):0.1:max(FF)]';
p12_mf_fixed = 1./(1+exp(pmode.a12-pmode.b12.*(FF_grid)+pmode.c12.*0*mean(MF)));
p21_mf_fixed = 1./(1+exp(pmode.a21+pmode.b21.*(FF_grid)-pmode.c21.*0*mean(MF)));



% Transition probabilities when coefficients have gamma prior
MF_grid = [-6:0.1:6]'; %[min(MF):0.1:max(MF)]';
p12_ff_fixed = 1./(1+exp(pmode.a12-pmode.b12.*0*mean(FF)+pmode.c12.*(MF_grid)));
p21_ff_fixed = 1./(1+exp(pmode.a21+pmode.b21.*0*mean(FF)-pmode.c21.*(MF_grid)));


%%
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
plot(MF_grid,p12_ff_fixed,'r'); hold on
plot(-FF_grid,p12_mf_fixed,'r--'); hold on
hleg = legend('$p_{12}$: Macro Factor','$p_{12}$: Financial Factor','Interpeter','Latex');
set(hleg,'interpreter','Latex','Location','North','FontSize',14); legend boxoff; 
%xlabel('Macroeconomic Factor','Interpreter','Latex','FontSize',10)
title('Normal-to-Bad','fontsize',12,'interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
axis tight; ylim([0, 0.3]);
xticks([-6 0 6])
xticklabels({'Worse','Average', 'Better'})
yticks([.1 .2 .3])
yticklabels({'0.1','0.1','0.3'})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig,[slides_folder model '_Pr12'],'-bestfit');


fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
plot(MF_grid,p21_ff_fixed,'b'); hold on
plot(-FF_grid,p21_mf_fixed,'b--'); hold on
hleg = legend('$p_{21}$: Macro Factor','$p_{21}$: Financial Factor','Interpeter','Latex');
set(hleg,'interpreter','Latex','Location','North','FontSize',14); legend boxoff; 
%xlabel('Macroeconomic Factor','Interpreter','Latex','FontSize',10)
title('Bad-to-Normal','fontsize',12,'interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
axis tight; ylim([0, 0.3]);
xticks([-6 0 6])
xticklabels({'Worse','Average', 'Better'})
yticks([.1 .2 .3])
yticklabels({'0.1','0.1','0.3'})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig,[slides_folder model '_Pr21'],'-bestfit');
%%


p12_fitted = 1./(1+exp(pmode.a12-pmode.b12.*(FF)+pmode.c12.*(MF)));
p12_fitted_ffonly = 1./(1+exp(pmode.a12-pmode.b12.*(FF)+pmode.c12.*0*(MF)));

p21_fitted = 1./(1+exp(pmode.a21+pmode.b21.*(FF)-pmode.c21.*(MF)));
p21_fitted_ffonly = 1./(1+exp(pmode.a21+pmode.b21.*(FF)-pmode.c21.*0*(MF)));

%% Transition probabilities

fig=figure
set(fig,'defaultAxesColorOrder',[right_color]);
subplot(211)
l1=plot(dates(sd:ed),p12_fitted,'color',colors2(1,:),'DisplayName','$\hat{p}_{12}$: normal-to-bad','LineWidth',3); hold on;
%l2=plot(dates(sd:ed),p12_fitted_ffonly,'k--','DisplayName','$\hat{p}_{12}$ with $f_t$ only');
axis tight; ylim([0, 1]);
rr=recessionplot;
% subplot(212)
% plot(dates(sd:ed),p21_fitted,'r'); hold on;
% plot(dates(sd:ed),p21_fitted_ffonly);
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
%title('Normal-to-Bad','fontsize',12,'interpreter','Latex')
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
yticks([.2 .4 .6 .8])
yticklabels({'0.2','0.4','0.6', '0.8'})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;


%%
% rrdates = [rr(1).Vertices(1) rr(1).Vertices(3) rr(2).Vertices(1) rr(2).Vertices(3)...
%     rr(3).Vertices(1) rr(3).Vertices(3) rr(4).Vertices(1) rr(4).Vertices(3)...
%     rr(5).Vertices(1) rr(5).Vertices(3) rr(6).Vertices(1) rr(6).Vertices(3)  ]'-14;
% 
% 
% [dates(find(dates>=rrdates(1) & dates<=rrdates(2)));
% dates(find(dates>=rrdates(1) & dates<=rrdates(2)));    

%print('-dpdf',fig,[slides_folder model '_Pr12_fitted'],'-bestfit');
%%

fig=figure
set(fig,'defaultAxesColorOrder',[right_color]);
subplot(211)
l1=plot(dates(sd:ed),p21_fitted,'color',colors2(16,:),'DisplayName','$\hat{p}_{21}$: bad-to-normal','LineWidth',3); hold on;
%l2=plot(dates(sd:ed),p21_fitted_ffonly,'k--','DisplayName','$\hat{p}_{21}$ with $f_t$ only');
axis tight; ylim([0, 1]);
rr=recessionplot;
% subplot(212)
% plot(dates(sd:ed),p21_fitted,'r'); hold on;
% plot(dates(sd:ed),p21_fitted_ffonly);
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
%title('Normal-to-Bad','fontsize',12,'interpreter','Latex')
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
yticks([.2 .4 .6 .8])
yticklabels({'0.2','0.4','0.6', '0.8'})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

%print('-dpdf',fig,[slides_folder model '_Pr21_fitted'],'-bestfit');
%%

%*** Compare with smoothed and filtered probabilities ****
fig=figure; clf;
%l1=plot(dates(sd:ed),p_reg2_filtered(sd:ed),'r--','LineWidth',1,'DisplayName','$P(s_t=2)$ filtered'); hold on;
l2=plot(dates(sd:ed),p_reg2_sim,'-','LineWidth',2,'DisplayName','Pseudo Real-Time Simulation of $P(s_t=2)$'); hold on
%l3=plot(dates(sd:ed),p_reg2_simfixed,'r--','LineWidth',2,'DisplayName','Pseudo Real-Time Simulation of $P(s_t=2)$');
ylim([0 1.0]);
rr=recessionplot;
%hleg = legend([l1 l2],'Orientation','Horizontal','Location','N','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight; 
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;

%print('-dpdf',fig,[slides_folder model '_PrReg2'],'-bestfit');
%%

%% Kernel smoothing function October 2008

md = '2007-Aug';
mydate = find(datenum(md,dataformat)==dates);


colors2 = cbrewer('div', 'RdYlGn', 18);

%%
close all;
fig=figure;
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
hold on
l1=plot(xi3,pdf3,'-.','Color',colors2(16,:),'LineWidth', 2,'DisplayName',['Normal regime']);
l2=plot(xi4,pdf4,'--','Color',colors2(1,:),'LineWidth', 2,'DisplayName',['Bad regime']);
%plot(xi,pdf,'k-','LineWidth', 2,'DisplayName',['Full']);
%plot([2.2 2.2],[0 0.6],'b-','LineWidth',0.8);c
hleg = legend([l1 l2],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: August-2007','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-8 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig,[slides_folder model '_Density_Aug2007_1'],'-bestfit');



%%
fig=figure;
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
hold on
l1=plot(xi3,pdf3,'-.','Color',colors2(16,:),'LineWidth', 2,'DisplayName',['Normal regime']);
l2=plot(xi4,pdf4,'--','Color',colors2(1,:),'LineWidth', 2,'DisplayName',['Bad regime']);
l3=plot(xi,pdf,'k-','LineWidth', 3,'DisplayName',['Full distribution']);
%plot([2.2 2.2],[0 0.6],'b-','LineWidth',0.8);c
hleg = legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: August-2007','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-8 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig,[slides_folder model '_Density_Aug2007_2'],'-bestfit');

%% Full distribution + Staff Forecast
fig=figure;
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
hold on
%l1=plot(xi3,pdf3,'-.','Color',colors2(16,:),'LineWidth', 2,'DisplayName',['Normal regime']);
%l2=plot(xi4,pdf4,'--','Color',colors2(1,:),'LineWidth', 2,'DisplayName',['Bad regime']);
l3=plot(xi,pdf,'k-','LineWidth', 3,'DisplayName',['Markov-Switching model']);
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
%annotation('textarrow',[0.75, .75],[0.5,.5],'Color','r','String','Staff');
hleg = legend([l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: August-2007','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-8 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig,[slides_folder model '_Density_Aug2007_3'],'-bestfit');

%%

fig=figure;
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
hold on
l1=plot(xi3,pdf3,'-.','Color',colors2(16,:),'LineWidth', 3,'DisplayName',['Normal regime']);
l2=plot(xi4,pdf4,'--','Color',colors2(1,:),'LineWidth', 3,'DisplayName',['Bad regime']);
l3=plot(xi,pdf,'-','Color',[0.8 0.8 0.8],'LineWidth', 3,'DisplayName',['Markov-Switching model']);
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
%plot([realized(mydate) realized(mydate)],[0 0.7],'-','Color',colors(1,:),'LineWidth',1,'DisplayName',[md ', Realized'])
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
%text(-.65,0.68,'$\mbox{realized} \rightarrow$','Interpreter','Latex','FontSize',12,'Color',colors2(1,:))
%annotation('textarrow',[0.75, .75],[0.5,.5],'Color','r','String','Staff');
hleg = legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: August-2007','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-8 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig,[slides_folder model '_Density_Aug2007_4'],'-bestfit');


%%

tind_use = find(datenum(char(md),dataformat)==dates);
gc = [0.8 0.8 0.8];
fig=figure;
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
hold on
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
%plot([realized(mydate) realized(mydate)],[0 0.7],'-','Color',gc,'LineWidth',0.8,'DisplayName',[md ', Realized'])
%l1=plot(xi3,pdf3,'-.','Color',gc,'LineWidth', 2,'DisplayName',['Normal regime']);
%l2=plot(xi4,pdf4,'--','Color',gc,'LineWidth', 2,'DisplayName',['Bad regime']);
l3=plot(xi,pdf,'k-','LineWidth', 4,'DisplayName',['Markov-Switching model (MS)']);
l4=plot(ResMatch.YY,PST(tind_use,:)','--','Color',dcolor(1,:),'LineWidth',4,'DisplayName','Quantile Regression model (QR)');
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
%text(-0.65,0.68,'$\mbox{realized} \rightarrow$','Interpreter','Latex','FontSize',12,'Color',gc)
%annotation('textarrow',[0.75, .75],[0.5,.5],'Color','r','String','Staff');
hleg = legend([l3 l4],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: August-2007','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-8 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig,[slides_folder model '_Density_Aug2007_5'],'-bestfit');
%% COVID-19 DENSITIES

fig=figure;
[pdf,xi]=ksdensity(y_erg_barcovid_mar13(1,:));
%[pdf3,xi3]=ksdensity(y_erg_barcovid_apr2(1,:));

hold on
l1=plot(xi,pdf,'k-','LineWidth', 4,'DisplayName',['March-13']);
%l2=plot(xi3,pdf3,'LineWidth', 4,'Color',dcolor(1,:),'DisplayName','April-2');
hleg = legend([l1],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: March-2020','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

%print('-dpdf',fig,[slides_folder model '_Density_Mar2020_1'],'-bestfit');
%%
%% COVID-19 DENSITIES

fig=figure;
[pdf,xi]=ksdensity(y_erg_barcovid_mar13(1,:));
[pdf3,xi3]=ksdensity(y_erg_barcovid_apr2(1,:));

hold on
l1=plot(xi,pdf,'k-','LineWidth', 4,'DisplayName',['March-13']);
l2=plot(xi3,pdf3,'LineWidth', 4,'Color',dcolor(1,:),'DisplayName','April-2');
hleg = legend([l1 l2],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: March-2020','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

%print('-dpdf',fig,[slides_folder model '_Density_Mar2020_2'],'-bestfit');


%% COVID-19 DENSITIES with QR
QR_COVID_MAR = load('PDF_vintage___Period=2020-Mar___Scenario=Mar13')
QR_COVID_APR = load('PDF_vintage___Period=2020-Mar___Scenario=Apr2')
fig=figure;
[pdf,xi]=ksdensity(y_erg_barcovid_mar13(1,:));
[pdf3,xi3]=ksdensity(y_erg_barcovid_apr2(1,:));

hold on
l1=plot(xi,pdf,'k-','LineWidth', 3,'DisplayName',['March-13 (MS)']);
l2=plot(xi3,pdf3,'LineWidth', 3,'Color',dcolor(1,:),'DisplayName','April-2 (MS)');
l3=plot(QR_COVID_MAR.vinpdf.X,QR_COVID_MAR.vinpdf.PDF,'k-.','LineWidth', 3,'DisplayName',['March-13 (QR)']);
l4=plot(QR_COVID_APR.vinpdf.X,QR_COVID_APR.vinpdf.PDF,'-.','LineWidth', 3,'Color',dcolor(1,:),'DisplayName','April-2 (QR)');
hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: March-2020','Interpreter','Latex','FontSize',9)
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

%print('-dpdf',fig,[slides_folder model '_Density_Mar2020_MS_QR'],'-bestfit');
