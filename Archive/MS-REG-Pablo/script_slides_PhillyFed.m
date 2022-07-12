 close all

% Figure 1-2
% August 2007 densities
%--------------------------------------------------------------------------

% Options for predictive densities
nDraws      = 10;           % Number of simulations per parameter draw
nParamDraws = options.N;    % Number of parameter draws

% Create data structure
data_use.FF_full       = FF_full;
data_use.MF_full       = MF_full;
data_use.TR_full       = TR_full;
data_use.prob_reg2_new = prob_reg2_new;
data_use.dates         = dates;
data_use.dates_full    = dates_full;
data_use.const         = const;
data_use.normal        = normal;

% Dates for density
md = '2007-Aug';
mydate = find(datenum(md,dataformat)==dates);

% Compute densities for MS model
% results correspond to the posterior parameter draws
outAug2007 = fGetDensityMS(md,data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);

% Predictive distributions
[pdf,xi]=ksdensity(outAug2007.y_erg_bar);
[pdf3,xi3]=ksdensity(outAug2007.y_erg_bar_1);
[pdf4,xi4]=ksdensity(outAug2007.y_erg_bar_2);


%% Define colors
colors = cbrewer('div', 'RdBu', 8);
mcolor = colors(8,:); % Color for QR density
ocolor = colors(1,:); % Color for MS density
dcolor = [mcolor;ocolor];
colors2 = cbrewer('div', 'RdYlGn', 18);

%% QR quantiles
load(['quantiles_QR']);

if exist('ResMatch')==0
    load('results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

PST = ResMatch.PST;
Qplot = [0.25 0.5 0.75];
%%
fig1a=figure;hold on
l3=plot(xi,pdf,'k-','LineWidth', 3,'DisplayName',['Markov-Switching model']);
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
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
set(fig1a,'PaperOrientation','portrait');
set(fig1a, 'PaperSize', figSize);
set(fig1a, 'PaperUnits', 'inches');
set(fig1a, 'Units','inches');
set(fig1a, 'PaperPositionMode', 'auto');
set(fig1a, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig1a,[slides_folder model '_Density_Aug2007_3'],'-bestfit');



fig1b=figure;
hold on
l1=plot(xi3,pdf3,'-.','Color',colors2(16,:),'LineWidth', 3,'DisplayName',['Normal regime']);
l2=plot(xi4,pdf4,'--','Color',colors2(1,:),'LineWidth', 3,'DisplayName',['Bad regime']);
l3=plot(xi,pdf,'-','Color',[0.8 0.8 0.8],'LineWidth', 3,'DisplayName',['Markov-Switching model']);
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
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
set(fig1b,'PaperOrientation','portrait');
set(fig1b, 'PaperSize', figSize);
set(fig1b, 'PaperUnits', 'inches');
set(fig1b, 'Units','inches');
set(fig1b, 'PaperPositionMode', 'auto');
set(fig1b, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig1b,[slides_folder model '_Density_Aug2007_4'],'-bestfit');


% Figure 2
% August 2007 Quantile comparison MS and QR models
%--------------------------------------------------------------------------
tind_use = find(datenum(char(md),dataformat)==dates);
gc = [0.8 0.8 0.8];
fig2=figure;
hold on
plot([2.2 2.2],[0 0.7],'Color',colors2(16,:),'LineWidth',0.8);
l3=plot(xi,pdf,'k-','LineWidth', 4,'DisplayName',['Markov-Switching model (MS)']);
l4=plot(ResMatch.YY,PST(tind_use,:)','--','Color',dcolor(1,:),'LineWidth',4,'DisplayName','Quantile Regression model (QR)');
text(2.3,0.68,'$\leftarrow \mbox{staff''s forecast}$','Interpreter','Latex','FontSize',14,'Color',colors2(16,:))
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
set(fig2,'PaperOrientation','portrait');
set(fig2, 'PaperSize', figSize);
set(fig2, 'PaperUnits', 'inches');
set(fig2, 'Units','inches');
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig2,[slides_folder model '_Density_Aug2007_5'],'-bestfit');

%%
% Figure 3-4
% Endogenous regime transition probabilities
%--------------------------------------------------------------------------

p12_fitted = 1./(1+exp(pmode.a12-pmode.b12.*(FF)+pmode.c12.*(MF)));
p21_fitted = 1./(1+exp(pmode.a21+pmode.b21.*(FF)-pmode.c21.*(MF)));

% Transition probabilities

fig3=figure
set(fig3,'defaultAxesColorOrder',[right_color]);
subplot(211)
l1=plot(dates(sd:ed),p12_fitted,'color',colors2(1,:),'DisplayName','$\hat{p}_{12}$: normal-to-bad','LineWidth',3); hold on;
axis tight; ylim([0, 1]);
rr=recessionplot;
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
yticks([.2 .4 .6 .8])
yticklabels({'0.2','0.4','0.6', '0.8'})
set(fig3,'PaperOrientation','portrait');
set(fig3, 'PaperSize', figSize);
set(fig3, 'PaperUnits', 'inches');
set(fig3, 'Units','inches');
set(fig3, 'PaperPositionMode', 'auto');
set(fig3, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

fig4=figure
set(fig4,'defaultAxesColorOrder',[right_color]);
subplot(211)
l1=plot(dates(sd:ed),p21_fitted,'color',colors2(16,:),'DisplayName','$\hat{p}_{21}$: bad-to-normal','LineWidth',3); hold on;
axis tight; ylim([0, 1]);
rr=recessionplot;
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex'); 
yticks([.2 .4 .6 .8])
yticklabels({'0.2','0.4','0.6', '0.8'})
set(fig4,'PaperOrientation','portrait');
set(fig4, 'PaperSize', figSize);
set(fig4, 'PaperUnits', 'inches');
set(fig4, 'Units','inches');
set(fig4, 'PaperPositionMode', 'auto');
set(fig4, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;

%%
% Figure 5: Pseudo Real-Time Simulation
%
%--------------------------------------------------------------------------
fig5=figure; clf;
l2=plot(dates(sd:ed),p_reg2_sim,'-','LineWidth',2,'DisplayName','Pseudo Real-Time Simulation of $P(s_t=2)$'); hold on
ylim([0 1.0]);
rr=recessionplot;
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
set(fig5,'PaperOrientation','portrait');
set(fig5, 'PaperSize', figSize);
set(fig5, 'PaperUnits', 'inches');
set(fig5, 'Units','inches');
set(fig5, 'PaperPositionMode', 'auto');
set(fig5, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;
%print('-dpdf',fig5,[slides_folder model '_PrReg2'],'-bestfit');
%%
% Figure 6: Evolution of Growth-at-Risk
%
%--------------------------------------------------------------------------

% Compute percentiles
dY_sim_25 = prctile(dY_sim',25)';
dY_sim_75 = prctile(dY_sim',75)';
dY_sim_10 = prctile(dY_sim',10)';
dY_sim_90 = prctile(dY_sim',90)';

% Set colors
colors = cbrewer('div', 'RdYlBu', 64);
numticks = 48;
fig6=figure; clf;
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
set(fig6,'PaperOrientation','portrait');
set(fig6, 'PaperSize', figSize);
set(fig6, 'PaperUnits', 'inches');
set(fig6, 'Units','inches');
set(fig6, 'PaperPositionMode', 'auto');
set(fig6, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;
% Save file
%print('-dpdf',fig6,[slides_folder model '_Quantiles_MS'],'-bestfit');

%%
% Figure 7: MS and QR capture Growth-at-Risk
%
% YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
%--------------------------------------------------------------------------
fig7=figure; clf;
hold on
% Plot percentiles from MS
l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'--','Color',colors(15,:),'LineWidth', 3,'DisplayName','MS 25th');
l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'--','Color',colors(55,:),'LineWidth', 3,'DisplayName','MS 75th');
% plot Quantiles from QR model
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
set(fig7,'PaperOrientation','portrait');
set(fig7, 'PaperSize', figSize);
set(fig7, 'PaperUnits', 'inches');
set(fig7, 'Units','inches');
set(fig7, 'PaperPositionMode', 'auto');
set(fig7, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;

% Save file
%print('-dpdf',fig7,[slides_folder model '_Quantiles_MS_QR'],'-bestfit');

%%
