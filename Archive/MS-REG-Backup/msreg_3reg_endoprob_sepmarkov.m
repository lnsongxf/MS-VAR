%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear
close all
% clc()

saveit=1;
save_mcmc=1;

addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));

%% Load RISE

rise_startup()

%% load the data

load data

%% plot the data

figure;
subplot(3,1,1)
plot(db.GDPGH+db.TRENDH,'linewidth',2)
title('Average Future GDP Growth')
subplot(3,1,2)
plot(db.FF,'linewidth',2)
title('Financial Factor')
subplot(3,1,3)
plot(db.MF,'linewidth',2)
title('Macroeconomic Factor')

varlist={'FF','MF','GDPGH'}';

%% set up Markov chains

switches_coef = {'c(3)','a0(3)'};
switches_vol = {'s(3)'};
model = 'SVAR_StateDependent_3Regimes_SepMarkov';

log_folder = ['Results/Logs/' model '/'];
if exist(log_folder,'dir')==0
    mkdir(log_folder)
end
fig_folder = ['Results/Figures/' model '/'];
if exist(fig_folder,'dir')==0
    mkdir(fig_folder)
end
mcmc_folder = ['Results/MCMC/' model '/'];
if exist(mcmc_folder,'dir')==0
    mkdir(mcmc_folder)
end


markov_chains=struct('name','coef',...
    'number_of_states',3,...
    'controlled_parameters',{switches_coef},...
    'endogenous_probabilities',{{
        'coef_tp_1_2=(1/(1+exp(-b12*(FF)-c12*(MF))))*(1/3)'
        'coef_tp_1_3=(1/(1+exp(-b13*(FF)-c13*(MF))))*(1/3)'
        'coef_tp_2_1=(1/(1+exp(-b21*(FF)-c21*(MF))))*(1/3)'
        'coef_tp_2_3=(1/(1+exp(-b23*(FF)-c23*(MF))))*(1/3)'
        'coef_tp_3_1=(1/(1+exp(-b31*(FF)-c31*(MF))))*(1/3)'
        'coef_tp_3_2=(1/(1+exp(-b32*(FF)-c32*(MF))))*(1/3)'
        }},...
    'probability_parameters',{{'b12','c12',...
                               'b13','c13',...
                               'b21','c21',...
                               'b23','c23',...
                               'b31','c31',...
                               'b32','c32'}});
   
markov_chains(2)=struct('name','vol',...
    'number_of_states',3,...
    'controlled_parameters',{switches_vol},...
    'endogenous_probabilities',{{
        'vol_tp_1_2=(1/(1+exp(-b12_vol*(FF)-c12_vol*(MF))))*(1/3)'
        'vol_tp_1_3=(1/(1+exp(-b13_vol*(FF)-c13_vol*(MF))))*(1/3)'
        'vol_tp_2_1=(1/(1+exp(-b21_vol*(FF)-c21_vol*(MF))))*(1/3)'
        'vol_tp_2_3=(1/(1+exp(-b23_vol*(FF)-c23_vol*(MF))))*(1/3)'
        'vol_tp_3_1=(1/(1+exp(-b31_vol*(FF)-c31_vol*(MF))))*(1/3)'
        'vol_tp_3_2=(1/(1+exp(-b32_vol*(FF)-c32_vol*(MF))))*(1/3)'
        }},...
    'probability_parameters',{{'b12_vol','c12_vol',...
                               'b13_vol','c13_vol',...
                               'b21_vol','c21_vol',...
                               'b23_vol','c23_vol',...
                               'b31_vol','c31_vol',...
                               'b32_vol','c32_vol'}});
   


%% Create the VAR

nlags=1;

exog={};

constant=true;

panel=[];

sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
lin_restr={
    % first equation or "financial factor" equation
    %----------------------------------
    'a0(1,GDPGH)=0'
    'a1(1,GDPGH)=0'
    %'a1(1,MF)=0'
    % second equation or "macroeconomic factor" equation
    %-----------------------------------
    'a0(2,GDPGH)=0'
    'a0(2,FF)=0'
    'a1(2,GDPGH)=0'
    % third equation or "GDP Growth" equation
    %----------------------------------
    'a1(3,GDPGH)=0'
    'a1(3,FF)=0'
    'a1(3,MF)=0'
    };

restrictions=[lin_restr];

%% set priors

% priors for the VAR coefficients
var_prior=svar.prior_template();
%--------------------------------
% Sims and Zha (1998) Normal-Wishart Prior Hyperparameters
% L1 : Overall tightness
% L2 : Cross-variable specific variance parameter
% L3 : Speed at which lags greater than 1 converge to zero
% L4 : tightness on deterministic/exogenous terms
% L5 : covariance dummies(omega)
% L6 : co-persistence (lambda)
% L7 : Own-persistence (mu)
var_prior.type='sz';
var_prior.L1=1;
var_prior.L5=1; 
var_prior.L6=0; % cointegration
var_prior.L7=0; % unit root
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both

% priors for the coef transition probabilities
%------------------------------------------------
switch_prior=struct();
% coefficients
switch_prior.b12={0,0,2,'normal'};
switch_prior.c12={0,0,2,'normal'};
switch_prior.b13={0,0,2,'normal'};
switch_prior.c13={0,0,2,'normal'};
switch_prior.b21={0,0,2,'normal'};
switch_prior.c21={0,0,2,'normal'};
switch_prior.b23={0,0,2,'normal'};
switch_prior.c23={0,0,2,'normal'};
switch_prior.b31={0,0,2,'normal'};
switch_prior.c31={0,0,2,'normal'};
switch_prior.b32={0,0,2,'normal'};
switch_prior.c32={0,0,2,'normal'};
% volatilities
switch_prior.b12_vol={0,0,2,'normal'};
switch_prior.c12_vol={0,0,2,'normal'};
switch_prior.b13_vol={0,0,2,'normal'};
switch_prior.c13_vol={0,0,2,'normal'};
switch_prior.b21_vol={0,0,2,'normal'};
switch_prior.c21_vol={0,0,2,'normal'};
switch_prior.b23_vol={0,0,2,'normal'};
switch_prior.c23_vol={0,0,2,'normal'};
switch_prior.b31_vol={0,0,2,'normal'};
switch_prior.c31_vol={0,0,2,'normal'};
switch_prior.b32_vol={0,0,2,'normal'};
switch_prior.c32_vol={0,0,2,'normal'};



prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

%% Find posterior mode

%options = optimset('TolX',1e-12);
if exist([mcmc_folder 'sv.mat'],'file')==2
    load([mcmc_folder 'sv.mat'])    
else
    %sv=estimate(sv0,db,{'1973M1','2019M5'},prior,restrictions,'fmincon',false,'optimset',options);
    sv=estimate(sv0,db,{'1973M1','2019M5'},prior,restrictions,'fmincon');
    if save_mcmc==1
        save([mcmc_folder 'sv.mat'],'sv')
    end
end

%% estimates

pmode=posterior_mode(sv)

%% Printing estimates

dfile = ([log_folder 'structural.txt']);
if exist(dfile, 'file') 
    delete(dfile); 
end
diary(dfile)
print_structural_form(sv)
diary off

%% Printing solution

dfile = ([log_folder 'reduced.txt']);
if exist(dfile, 'file') 
    delete(dfile); 
end
diary(dfile)
print_solution(sv)
diary off



%% plots probabilities against data

%{
%plot_probabilities(sv)
plot_data_against_probabilities(sv,'regime')
%}

%% Ergodic Fitted Values

[Resids,Fits]=residuals(sv);
fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
fit3 = Fits.GDPGH.data(:,3)+db.TRENDH.data(nlags+1:end);
fit4 = Fits.GDPGH.data(:,4)+db.TRENDH.data(nlags+1:end);
fit5 = Fits.GDPGH.data(:,5)+db.TRENDH.data(nlags+1:end);
fit6 = Fits.GDPGH.data(:,6)+db.TRENDH.data(nlags+1:end);
fit7 = Fits.GDPGH.data(:,7)+db.TRENDH.data(nlags+1:end);
fit8 = Fits.GDPGH.data(:,8)+db.TRENDH.data(nlags+1:end);
fit9 = Fits.GDPGH.data(:,9)+db.TRENDH.data(nlags+1:end);

[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;    
p_reg3 = f.smoothed_regime_probabilities.regime3.data;
p_reg4 = f.smoothed_regime_probabilities.regime4.data;
p_reg5 = f.smoothed_regime_probabilities.regime5.data;    
p_reg6 = f.smoothed_regime_probabilities.regime6.data;    
p_reg7 = f.smoothed_regime_probabilities.regime7.data;
p_reg8 = f.smoothed_regime_probabilities.regime8.data;    
p_reg9 = f.smoothed_regime_probabilities.regime9.data;    


% coeff_fit1 = fit1.*p_reg1 + fit2.*p_reg2 + fit3.*p_reg3; 
% coeff_fit2 = fit4.*p_reg4 + fit5.*p_reg5 + fit6.*p_reg6; 
% coeff_fit3 = fit7.*p_reg7 + fit8.*p_reg8 + fit9.*p_reg9; 

inputformat = 'yyyy-MMM';
dates = datenum((datetime('1973-Jan','InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime('2019-May','InputFormat',inputformat)))';
dataformat  = 'yyyy-mmm';

FontSize = 16;
numticks = 48;  
figSize = [12 6];
linestyle = {'-','--',':'};

start_plot = '1973-Feb';
end_plot   = '2019-May';
sd = find(datenum(start_plot,dataformat)==dates);
ed = find(datenum(end_plot,dataformat)==dates);

y_fit = fit1.*p_reg1 + fit2.*p_reg2 + fit3.*p_reg3 + fit4.*p_reg4 + fit5.*p_reg5 + fit6.*p_reg6 + fit7.*p_reg7 + fit8.*p_reg8 + fit9.*p_reg9;

realized = db.GDPGH.data + db.TRENDH.data;
realized = realized(nlags+1:end);

fig=figure;
hold on
plot(dates(sd:ed),y_fit(sd:ed),'linewidth',2)
plot(dates(sd:ed),realized(sd:ed),'k','linewidth',2)
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
legend('Ergodic Fitted Value','Realized','interpreter','Latex')
legend boxoff
axis tight
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Ergodic_Fitted'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Ergodic_Fitted']));
end

fig=figure;
hold on
plot(dates,(realized(sd:ed))-y_fit(sd:ed),'r','linewidth',2)
plot(dates,zeros(length(dates(sd:ed))),'k-.','linewidth',2)
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
ylim([-6.5 6.5])
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Ergodic_Difference'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Ergodic_Difference']));
end


disp(round(sum((realized-y_fit).^2),2))


%% Regime Probabilities Against Data

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close

% coef 
% bad = regime 1
% normal = regime 3
% good = regime 2

% vol
% low = regime 3
% moderate = regime 1
% high = regime 2

coef_1 = [p_reg3 p_reg1 p_reg2];
coef_2 = [p_reg6 p_reg4 p_reg5];
coef_3 = [p_reg9 p_reg7 p_reg8];
reg_prob = [coef_2 coef_3 coef_1];

figSize2 = [14 12];
numticks = 48*2;
left_color = [0 0 0];

fig=figure;

% Good Regime - Low Variance
right_color = colors(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,1,[0.125,0.1]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg1_LB(sd:ed) p_reg1_M(sd:ed) p_reg1_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,1),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
%legend('Financial','Macroeconomic','Location','Best','interpreter','Latex')
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.1 spPos(1,2)-0.005 0.35 0.3],'String','\textbf{Good Regime \& Low Variance}','interpreter','Latex',titleSettings{:})

% Good Regime - Mild Variance
right_color = colors(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,2,[0.125,0.1]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg1_LB(sd:ed) p_reg1_M(sd:ed) p_reg1_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,2),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
%legend('Financial','Macroeconomic','Location','Best','interpreter','Latex')
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.05 spPos(1,2)-0.005 0.35 0.3],'String','\textbf{Good Regime \& Mild Variance}','interpreter','Latex',titleSettings{:})

% Good Regime - High Variance
right_color = colors(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,3,[0.125,0.1]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg1_LB(sd:ed) p_reg1_M(sd:ed) p_reg1_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,3),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
%legend('Financial','Macroeconomic','Location','Best','interpreter','Latex')
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.025 spPos(1,2)-0.005 0.35 0.3],'String','\textbf{Good Regime \& High Variance}','interpreter','Latex',titleSettings{:})

% Normal Regime - Low Variance
right_color = colors(1,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,4,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,4),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.1 spPos(1,2)-0.06 0.35 0.3],'String','\textbf{Normal Regime \& Low Variance}','interpreter','Latex',titleSettings{:})

% Normal Regime - Mild Variance
right_color = colors(1,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,5,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,5),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.05 spPos(1,2)-0.06 0.35 0.3],'String','\textbf{Normal Regime \& Mild Variance}','interpreter','Latex',titleSettings{:})

% Normal Regime - High Variance
right_color = colors(1,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,6,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,6),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.025 spPos(1,2)-0.06 0.35 0.3],'String','\textbf{Normal Regime \& High Variance}','interpreter','Latex',titleSettings{:})

% Bad Regime - Low Variance
right_color = colors(2,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,7,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg3_LB(sd:ed) p_reg3_M(sd:ed) p_reg3_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,7),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.1 spPos(1,2)-0.12 0.35 0.3],'String','\textbf{Bad Regime \& Low Variance}','interpreter','Latex',titleSettings{:})

% Bad Regime - Mild Variance
right_color = colors(2,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,8,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg3_LB(sd:ed) p_reg3_M(sd:ed) p_reg3_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,8),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.05 spPos(1,2)-0.12 0.35 0.3],'String','\textbf{Bad Regime \& Mild Variance}','interpreter','Latex',titleSettings{:})

% Bad Regime - High Variance
right_color = colors(2,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,3,9,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg3_LB(sd:ed) p_reg3_M(sd:ed) p_reg3_UB(sd:ed)],2)
plot(dates(sd:ed),reg_prob(sd:ed,9),'Linewidth',1.5)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
axis tight
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
hold off
spPos = cat(1,sfig.Position);
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',FontSize,'FontWeight','bold', 'FontName', 'Times New Roman'};
annotation('textbox','Position',[spPos(1,1)-0.025 spPos(1,2)-0.12 0.35 0.3],'String','\textbf{Bad Regime \& High Variance}','interpreter','Latex',titleSettings{:})


set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize2);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize2(1)/2 figSize2(2)/15 figSize2(1) figSize2(2)]);
tightfigadv;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end


%% Time Series of Quantiles - QR Comparison

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(2,:);colr(1,:);colr(5,:)];
close

load(['quantiles_QR']);
Qplot = [0.25 0.5 0.75];
quantiles_dist = 0.1:0.05:0.9; % quantiles to construct distribution
indq = ismember(round(quantiles_dist,2),Qplot);
indq = find(indq==1);

fig=figure;
%yyaxis left
hold on
plot(dates(sd:ed), fit1(sd:ed),'Color',colors(1,:),'LineWidth', 1.5)
plot(dates(sd:ed), fit7(sd:ed),'Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed), fit4(sd:ed),'Color',colors(3,:),'LineWidth', 1.5)
j=0;
for iquantile=indq
    j=j+1;
    plot(dates(sd:ed),YQ_dist_fut(sd:ed,iquantile),'-.','Color',colors(j,:),'LineWidth',1.5);
    hold on
end
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {'Bad Regime','Normal Regime','Good Regime',[num2str(Qplot(1)*100) 'th Quantile'],'Median',[num2str(Qplot(3)*100) 'th Quantile']};
%yyaxis right
%plot(dates(sd:ed),p_reg1(sd:ed),'k-.','LineWidth', 1.25)
%ylabel('Probability','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
legend(leg,'Orientation','Vertical','Location','South','interpreter','Latex')
legend boxoff
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
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Quantiles_MSFitted_QR'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
end

%% Impulse responses

%{

myirfs=irf(sv);

snames=fieldnames(myirfs);

figure('name','Simple impulse response functions')

iter=0;

jter=0;

for ishk=1:numel(snames)
    
    shk=snames{ishk};
    
    for iv=1:numel(varlist)
        
        v=varlist{iv};
        
        iter=iter+1;
        
        subplot(3,3,iter)
        
        plot(myirfs.(shk).(v),'linewidth',2)
        
        if iter<=3,title(v),end
        
        if any(iter==[1,4,7])
            
            jter=jter+1;
            
            ylabel([shk,'(',varlist{jter},' shock)'])
            
        end
        
    end
    
end
%}

%% Generalized impulse responses

%{

shock_names=[];

irf_periods=[];

params=[];

girf_setup=struct();% girf_setup=struct('nsims',300);

mygirfs=irf(sv,shock_names,irf_periods,params,girf_setup);

figure('name','Generalized impulse response functions')

iter=0;

jter=0;

for ishk=1:numel(snames)
    
    shk=snames{ishk};
    
    for iv=1:numel(varlist)
        
        v=varlist{iv};
        
        iter=iter+1;
        
        subplot(3,3,iter)
        
        plot(mygirfs.(shk).(v),'linewidth',2)
        
        if iter<=3,title(v),end
        
        if any(iter==[1,4,7])
            
            jter=jter+1;
            
            ylabel([shk,'(',varlist{jter},' shock)'])
            
        end
        
    end
    
end

%}

%% Posterior sampling

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);
vcov = utils.cov.nearest(vcov); %ensure positive definite  VCV

options=struct();
options.alpha=0.234;
options.thin=10;
options.burnin=10^3;
options.N=2*10^4;
options.nchain=1;

lb(~isfinite(lb))=-500;
ub(~isfinite(ub))=500;

if exist([mcmc_folder 'results.mat'],'file')==2
    load([mcmc_folder 'results.mat'])    
else
    results=mh_sampler(ff,lb,ub,options,x0,vcov);
    if save_mcmc==1
        save([mcmc_folder 'results.mat'],'results')
    end
end


%% Marginal data density

%{
mddobj=mdd(results,ff,lb,ub);

dfile =[log_folder 'mdd.txt'];
if exist(dfile, 'file') 
    delete(dfile); 
end
diary(dfile)
 
fprintf('Importance sampling::%0.2f\n',is(mddobj,[],mdd.global_options))

fprintf('Reciprocal Importance sampling::%0.2f\n',ris(mddobj,[],mdd.global_options))

fprintf('Meng and Wong''s bridge::%0.2f\n',bridge(mddobj,true,mdd.global_options))

fprintf('Ulrich Mueller::%0.2f\n',mueller(mddobj,[],mdd.global_options))

fprintf('Laplace::%0.2f\n',laplace(mddobj))

fprintf('Sims, Waggoner and Zha::%0.2f\n',swz(mddobj,[],mdd.global_options))

fprintf('Laplace MCMC::%0.2f\n',laplace_mcmc(mddobj))

fprintf('Chib and Jeliazkov::%0.2f\n',cj(mddobj,[],mdd.global_options))

diary off
%}

%% irfs distributions

%{
shock_names=[];2

params=[results.pop.x];

myirfs2=irf(sv,shock_names,[],params);

% plotting of fan chargs
for ireg=1:sv.nregs
    
    reg=sprintf('regime_%d',ireg);
    
    figure('name',['Simple impulse response functions(',reg,')'])
    
    iter=0;
    
    jter=0;
    
    for ishk=1:numel(snames)
        
        shk=snames{ishk};
        
        for iv=1:numel(varlist)
            
            v=varlist{iv};
            
            iter=iter+1;
            
            subplot(3,3,iter)
            
            ffr=fanchart(myirfs2.(shk).(v),[30,50,70,90]);
            
            plot_fanchart(ffr.(reg),'r')
            
            if iter<=3,title(v),end
            
            if any(iter==[1,4,7])
                
                jter=jter+1;
                
                ylabel([shk,'(',varlist{jter},' shock)'])
                
            end
            
        end
        
    end
    
end
%}


%% Out-of sample forecasts at the mode

%{

shock_uncertainty=false;

nsteps=12;

mycast=forecast(sv,[],[],[],nsteps,shock_uncertainty);

do_plot_unconditional_forecasts(mycast,sv)
%}

%% Residuals and Fitted Values

%{
[Resids,Fits]=residuals(sv);

% Residuals
res_names=fieldnames(Resids);

titel='Residuals';

figure('name',titel);
    
thisname='shock_1';

plot(Resids.(thisname),'linewidth',2)

title(thisname)
        
[~,h]=sup_label(titel,'t');

set(h,'fontsize',12)

% Fitted Values
fit_names=varlist{end};
titel='Fitted Values';
figure('name',titel);    
thisname=fit_names;
plot(Fits.(thisname),'linewidth',2)
hold on
plot(db.GDPGH,'k','linewidth',2)
hold off
title(thisname)
    
    
[~,h]=sup_label(titel,'t');

set(h,'fontsize',12)

% Ergodic Fitted Values
y_fit = fit1.*p_reg1 + fit2.*p_reg2 + fit3.*p_reg3 + fit4.*p_reg4 + fit5.*p_reg5 + fit6.*p_reg6 + fit7.*p_reg7 + fit8.*p_reg8 + fit9.*p_reg9;
start='1973M1';
yfit=ts(start,y_fit);

titel='Ergodic Fitted Values';
figure('name',titel);
thisname=varlist{end};
plot(yfit,'linewidth',2)
hold on
plot(db.(thisname)+db.TRENDH,'k','linewidth',2)
hold off
legend('Ergodic Fitted Values','Realized')
legend boxoff
xrotate(45)
title(thisname)
%}


%% Posterior sampling

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);

options=struct();
options.alpha=0.234;
options.thin=10;
options.burnin=10^3;
options.N=2*10^4;
options.nchain=1;

lb(~isfinite(lb))=-500;
ub(~isfinite(ub))=500;

if exist([mcmc_folder 'results.mat'],'file')==2
    load([mcmc_folder 'results.mat'])    
else
    results=mh_sampler(ff,lb,ub,options,x0,vcov);
    if save_mcmc==1
        save([mcmc_folder 'results.mat'],'results')
    end
end

%% MCMC diagnostics

pnames=fieldnames(pmode);

a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

if options.nchain>1
    chain_id = 2;
    myresults = results{1,chain_id};
    mcmcobj=mcmc(myresults,pnames,{1:options.N},a2tilde_to_a);
else
    mcmcobj=mcmc(results,pnames,{1:options.N},a2tilde_to_a);    
end

myList={'coef_tp_1_2','coef_tp_2_1'};

chain_id=[];

%{
for ii=1:numel(myList)
    
    v=myList{ii};
    
    figure('name',['Univariate convergence diagnostics for ',v]);
    
    subplot(2,2,1),autocorrplot(mcmcobj,v,chain_id,40); title('autocorrelogram')
    
    subplot(2,2,2),densplot(mcmcobj,v,chain_id,250); title('density')
    
    subplot(2,2,3),meanplot(mcmcobj,v,chain_id); title('recursive mean')
    
    %subplot(3,2,4),psrf_plot(mcmcobj,v); title('pot. scale red. factor')
    
    subplot(2,2,4),traceplot(mcmcobj,v,chain_id,20); title('trace')
    
    pause
    
    close
    
end
%}

%% Marginal data density

mddobj=mdd(results,ff,lb,ub);

dfile =[log_folder 'mdd.txt'];
if exist(dfile, 'file') 
    delete(dfile); 
end
diary(dfile)

fprintf('Importance sampling::%0.2f\n',is(mddobj,[],mdd.global_options))

fprintf('Reciprocal Importance sampling::%0.2f\n',ris(mddobj,[],mdd.global_options))

fprintf('Meng and Wong''s bridge::%0.2f\n',bridge(mddobj,true,mdd.global_options))

fprintf('Ulrich Mueller::%0.2f\n',mueller(mddobj,[],mdd.global_options))

fprintf('Laplace::%0.2f\n',laplace(mddobj))

fprintf('Sims, Waggoner and Zha::%0.2f\n',swz(mddobj,[],mdd.global_options))

fprintf('Laplace MCMC::%0.2f\n',laplace_mcmc(mddobj))

fprintf('Chib and Jeliazkov::%0.2f\n',cj(mddobj,[],mdd.global_options))

diary off

%% Solve

params=[results.pop.x];
myparams=a2tilde_to_a(params);

T = length(db.GDPGH.data)-nlags;
p_reg1 = zeros(T,options.N);
p_reg2 = zeros(T,options.N);
p_reg3 = zeros(T,options.N);
p_reg4 = zeros(T,options.N);
p_reg5 = zeros(T,options.N);
p_reg6 = zeros(T,options.N);
p_reg7 = zeros(T,options.N);
p_reg8 = zeros(T,options.N);
p_reg9 = zeros(T,options.N);
y_erg = zeros(T,options.N);

parfor i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
    
    [Resids,Fits]=residuals(sv,myparams(:,i));
    [~,~,~,f]=filter(sv,myparams(:,i));
    
    fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
    fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
    fit3 = Fits.GDPGH.data(:,3)+db.TRENDH.data(nlags+1:end);
    fit4 = Fits.GDPGH.data(:,4)+db.TRENDH.data(nlags+1:end);
    fit5 = Fits.GDPGH.data(:,5)+db.TRENDH.data(nlags+1:end);
    fit6 = Fits.GDPGH.data(:,6)+db.TRENDH.data(nlags+1:end);
    fit7 = Fits.GDPGH.data(:,7)+db.TRENDH.data(nlags+1:end);
    fit8 = Fits.GDPGH.data(:,8)+db.TRENDH.data(nlags+1:end);
    fit9 = Fits.GDPGH.data(:,9)+db.TRENDH.data(nlags+1:end);
    resid1 = Resids.shock_1.data(:,1);
    resid2 = Resids.shock_1.data(:,2);
    resid3 = Resids.shock_1.data(:,3);
    resid4 = Resids.shock_1.data(:,4);
    resid5 = Resids.shock_1.data(:,5);
    resid6 = Resids.shock_1.data(:,6);
    resid7 = Resids.shock_1.data(:,7);
    resid8 = Resids.shock_1.data(:,8);
    resid9 = Resids.shock_1.data(:,9);
    y1 = fit1 + sqrt(myparams(end-2,i))*randn; 
    y2 = fit2 + sqrt(myparams(end-1,i))*randn; 
    y3 = fit3 + sqrt(myparams(end,i))*randn; 
    y4 = fit4 + sqrt(myparams(end-2,i))*randn; 
    y5 = fit5 + sqrt(myparams(end-1,i))*randn; 
    y6 = fit6 + sqrt(myparams(end,i))*randn; 
    y7 = fit7 + sqrt(myparams(end-2,i))*randn; 
    y8 = fit8 + sqrt(myparams(end-1,i))*randn; 
    y9 = fit9 + sqrt(myparams(end,i))*randn; 
    
    p_reg1(:,i) = f.smoothed_regime_probabilities.regime1.data;
    p_reg2(:,i) = f.smoothed_regime_probabilities.regime2.data;    
    p_reg3(:,i) = f.smoothed_regime_probabilities.regime3.data;    
    p_reg4(:,i) = f.smoothed_regime_probabilities.regime4.data;
    p_reg5(:,i) = f.smoothed_regime_probabilities.regime5.data;    
    p_reg6(:,i) = f.smoothed_regime_probabilities.regime6.data;    
    p_reg7(:,i) = f.smoothed_regime_probabilities.regime7.data;
    p_reg8(:,i) = f.smoothed_regime_probabilities.regime8.data;    
    p_reg9(:,i) = f.smoothed_regime_probabilities.regime9.data;    
    y_erg(:,i) = y1.*p_reg1(:,i) + y2.*p_reg2(:,i) + y3.*p_reg3(:,i)...
        + y4.*p_reg4(:,i) + y5.*p_reg5(:,i) + y6.*p_reg6(:,i) ...
        + y7.*p_reg7(:,i) + y8.*p_reg8(:,i) + y9.*p_reg9(:,i);
    
    
end

inputformat = 'yyyy-MMM';
dates = datenum((datetime('1973-Jan','InputFormat',inputformat)):calmonths(1):(datetime('2019-May','InputFormat',inputformat)))';
dataformat  = 'yyyy-mmm';


q10 = prctile(y_erg,10,2);
q50 = prctile(y_erg,50,2);
q90 = prctile(y_erg,90,2);


%% Time Series of Quantiles

%{
FontSize = 16;
numticks = 48;  
figSize = [12 6];

start_plot = '1973-Jan';
end_plot   = '2019-May';
sd = find(datenum(start_plot,dataformat)==dates);
ed = find(datenum(end_plot,dataformat)==dates);

colors = [[139,0,0]/255;[0,0,0]/255;[0,0,128]/255];

fig=figure;
hold on
plot(dates(sd:ed), q10(sd:ed),'Color',colors(1,:),'LineWidth', 2)
plot(dates(sd:ed), q50(sd:ed),'--','Color',colors(2,:),'LineWidth', 2)
plot(dates(sd:ed), q90(sd:ed),':','Color',colors(3,:),'LineWidth', 2)
hold off
leg = {'10th Quantile','Median','90th Quantile'};
ylabel('Percent','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
rr=recessionplot;
legend(leg,'Orientation','Horizontal','Location','South','interpreter','Latex')
legend boxoff
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
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Quantiles'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles']));
end
%}

%% Time Series of Quantiles - QR Comparison

%{
load(['quantiles_QR']);
Qplot = [0.10 0.5 0.90];
quantiles_dist = 0.1:0.05:0.9; % quantiles to construct distribution
indq = ismember(round(quantiles_dist,2),Qplot);
indq = find(indq==1);

q1 = prctile(y_erg,Qplot(1)*100,2);
q2 = prctile(y_erg,Qplot(2)*100,2);
q3 = prctile(y_erg,Qplot(3)*100,2);


fig=figure;
%yyaxis left
hold on
plot(dates(sd:ed), q1(sd:ed),'-.','Color',colors(1,:),'LineWidth', 2)
plot(dates(sd:ed), q2(sd:ed),'-.','Color',colors(2,:),'LineWidth', 2)
plot(dates(sd:ed), q3(sd:ed),'-.','Color',colors(3,:),'LineWidth', 2)
j=0;
for iquantile=1:3
    j = j+1;
    plot(dates(sd:ed),YQ_dist_fut(sd:ed,indq(iquantile)),'Color',colors(j,:),'LineWidth',2);
    hold on
end
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'g','linewidth',2)
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {[num2str(Qplot(1)*100) 'th - MS'],'Median - MS',[num2str(Qplot(3)*100) 'th - MS'],[num2str(Qplot(1)*100) 'th - QR'],'Median - QR',[num2str(Qplot(3)*100) 'th - QR'],'Realized'};
%yyaxis right
%plot(dates(sd:ed),p_reg1(sd:ed),'k-.','LineWidth', 1.25)
%ylabel('Probability','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
legend(leg,'Orientation','Vertical','Location','South','interpreter','Latex')
legend boxoff
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
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Quantiles_MSFitted'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted']));
end
%}

%% Time Series of Quantiles - Zoomed

%{
numticks = 12;  

start_plot = '2006-Jan';
end_plot   = '2010-Dec';
sdz = find(datenum(start_plot,dataformat)==dates);
edz = find(datenum(end_plot,dataformat)==dates);

fig=figure;
hold on
plot(dates(sdz:edz), q1(sdz:edz),'-.','Color',colors(1,:),'LineWidth', 2)
plot(dates(sdz:edz), q2(sdz:edz),'-.','Color',colors(2,:),'LineWidth', 2)
plot(dates(sdz:edz), q3(sdz:edz),'-.','Color',colors(3,:),'LineWidth', 2)
j=0;
for iquantile=1:3
    j = j+1;
    plot(dates(sd:ed),YQ_dist_fut(sd:ed,indq(iquantile)),'Color',colors(j,:),'LineWidth',2);
    hold on
end
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'g','linewidth',2)
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {[num2str(Qplot(1)*100) 'th - MS'],'Median - MS',[num2str(Qplot(3)*100) 'th - MS'],[num2str(Qplot(1)*100) 'th - QR'],'Median - QR',[num2str(Qplot(3)*100) 'th - QR'],'Realized'};
legend(leg,'Orientation','Vertical','Location','Best','interpreter','Latex')
legend boxoff
ax=gca;
ax.XTick = datenum(dates(sdz:numticks:edz));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sdz), dates(edz)])
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
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Quantiles_Zoom'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_Zoom']));
end
%}

%% Density for Specific Period

%{
colors = [0    0.4470    0.7410;    
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;    
        0.4940    0.1840    0.5560;...
        0.4660 0.6740 0.1888];

%[~,mydate] = min(q10);
%md = datestr(dates(mydate),dataformat);

md = '2008-Mar';
mydate = find(datenum(md,dataformat)==dates);
od = '2018-Dec';
mydate2 = find(datenum(od,dataformat)==dates);

mcolor = colors(1,:);
ocolor = colors(2,:);

% Histogram
fig=figure;
h1=histogram(y_erg(mydate,:),'Normalization','pdf','FaceColor',mcolor);
hold on
h2=histogram(y_erg(mydate2,:),'Normalization','pdf','FaceColor',ocolor);
hold off
axis tight
legend({md,od},'Location','NorthWest','interpreter','Latex')
legend boxoff
ylabel('PDF','fontsize',10,'interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Histogram'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Histogram']));
end

% Kernel smoothing function estimate of density
fig=figure;
[pdf,xi]=ksdensity(y_erg(mydate,:));
plot(xi,pdf,'Color',mcolor,'LineWidth', 2);
hold on
[pdf2,xi2]=ksdensity(y_erg(mydate2,:));
plot(xi2,pdf2,'Color',ocolor,'LineWidth', 2);
hold off
legend({md,od},'Location','NorthWest','interpreter','Latex')
legend boxoff
ylabel('PDF','fontsize',10,'interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'KernelDensity'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'KernelDensity']));
end

%}

%% Comparison with Quantile Regression Densities

close all;
colors = [0    0.4470    0.7410;    
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;    
        0.4940    0.1840    0.5560;...
        0.4660 0.6740 0.1888];
    
dcolor = [colors(5,:);colors(2,:);colors(1,:)];

if exist('ResMatch')==0
    load('/mcr/data_bfi/m1fxl02/3. Miscellanea/Coronavirus VIX/Text_Paper/Figures_Monthly/US/horizon_12/results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

% coef 
% bad = regime 1
% normal = regime 3
% good = regime 2

% vol
% low = regime 3
% moderate = regime 1
% high = regime 2

% Select Time Periods for Which to Plot PDFs
period_pdf = {'2005-Mar','2008-Apr','2011-Dec'}; 
for tt=1:length(period_pdf)
    tind(tt) = find(datenum(char(period_pdf(tt)),dataformat)==dates);
end

clc;
% good, low volatility
p1 = coef_2(:,1);
p1(tind(1))
% bad, high volatility
p3 = coef_1(:,3);
p3(tind(2))
% normal, mild volatility
p2 = coef_3(:,2);
p2(tind(3))



% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST(tind,:);

fig=figure;
xmin = -8; 
xmax = 8; 
for t=1:length(period_pdf)
        % plot figure for in-sample pdfs    
    if t==1
        ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    end
    % PDF Plot
    plot(ax1,ResMatch.YY,PST(tind(t),:)','-','Color',dcolor(t,:),'LineWidth',2,'DisplayName',[char(period_pdf(t)) ' - QR']);
    hold(ax1,'on')
    [pdf,xi]=ksdensity(y_erg(tind(t),:));
    plot(ax1,xi,pdf,'-.','LineWidth', 2,'Color',dcolor(t,:),'DisplayName',[char(period_pdf(t)) ' - MS-SVAR']);
end
set(ax1, 'YLim', [0, 1])
yl = ylim;
for t=1:length(period_pdf)
    plot(ax1,[realized(tind(t)) realized(tind(t))],[0 yl(2)],'--','Color',dcolor(t,:),'LineWidth',2,'DisplayName',[char(period_pdf(t)) ' - Realized'])
end        
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
set(ax1, 'XLim', [xmin, xmax])
set(ax1,'XTick',xmin:2:xmax)
set(ax1,'TickLabelInterpreter','Latex')
set(ax1, 'FontName', 'Times New Roman');
set(ax1, 'FontSize', FontSize);
set(ax1,'Layer','top')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'DensityQR'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end


%% Uncertainty Around Estimates

params=[results.pop.x];

params_M = prctile(params,50,2);
params_UB = prctile(params,95,2);
params_LB = prctile(params,5,2);

print_structural_form(sv)

pmode=posterior_mode(sv);
F = fieldnames(pmode);
C = num2cell(struct2cell(pmode),3);
C = cellfun(@(c)[c{:}],C,'uni',0);
params_pmode = vertcat(C{:});

params_M_Struct=a2tilde_to_a(params_M);
params_UB_Struct=a2tilde_to_a(params_UB);
params_LB_Struct=a2tilde_to_a(params_LB);

dfile =[log_folder 'CI.txt'];
if exist(dfile, 'file') 
    delete(dfile); 
end
diary(dfile)

disp('coef')
disp('bad = regime 1')
disp('normal = regime 3')
disp('good = regime 2')

disp('vol')
disp('low = regime 3')
disp('moderate = regime 1')
disp('high = regime 2')

pmode

params_CI = [params_pmode-(params_M_Struct-params_LB_Struct) params_pmode params_pmode+(params_UB_Struct-params_M_Struct)];

disp('Transition Probabilities')
round(params_CI(end-26:end-3,:),2)
disp('Coefficients')
round(params_CI(1:end-27,:),2)
disp('Volatilities')
round((params_CI(end-2:end,:)),2) % no need to take sqrt in structural form

diary off



