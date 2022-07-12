%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear
close all
% clc()

saveit = 1; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const = 0; % 1= have a constant in transition probability
normal = 0; % 1 = use normal distribution, 0 = gamma distribution

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

switches = {'c(3)','a0(3)','s(3)'};
model = 'SVAR_StateDependent_3Regimes';

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

if const==1
    if normal==1
        prob_fct = {{
                'sync_tp_1_2=(1/(1+exp(a12-b12*(FF)-c12*(MF))))*(1/3)'
                'sync_tp_1_3=(1/(1+exp(a13-b13*(FF)-c13*(MF))))*(1/3)'
                'sync_tp_2_1=(1/(1+exp(a21-b21*(FF)-c21*(MF))))*(1/3)'
                'sync_tp_2_3=(1/(1+exp(a23-b23*(FF)-c23*(MF))))*(1/3)'
                'sync_tp_3_1=(1/(1+exp(a31-b31*(FF)-c31*(MF))))*(1/3)'
                'sync_tp_3_2=(1/(1+exp(a32-b32*(FF)-c32*(MF))))*(1/3)'
                }};        
    elseif normal==0
        prob_fct = {{
                'sync_tp_1_2=(1/(1+exp(a12+b12*(FF)-c12*(MF))))*(1/3)'
                'sync_tp_1_3=(1/(1+exp(a13+b13*(FF)-c13*(MF))))*(1/3)'
                'sync_tp_2_1=(1/(1+exp(a21-b21*(FF)+c21*(MF))))*(1/3)'
                'sync_tp_2_3=(1/(1+exp(a23+b23*(FF)-c23*(MF))))*(1/3)'
                'sync_tp_3_1=(1/(1+exp(a31-b31*(FF)+c31*(MF))))*(1/3)'
                'sync_tp_3_2=(1/(1+exp(a32-b32*(FF)+c32*(MF))))*(1/3)'
                }};
    end
    prob_params = {{'a12','a13','a21','a23','a31','a32','b12','c12','b13','c13','b21','c21','b23','c23','b31','c31','b32','c32'}};
elseif const==0
    if normal==1
        prob_fct = {{
                'sync_tp_1_2=(1/(1+exp(-b12*(FF)-c12*(MF))))*(1/3)'
                'sync_tp_1_3=(1/(1+exp(-b13*(FF)-c13*(MF))))*(1/3)'
                'sync_tp_2_1=(1/(1+exp(-b21*(FF)-c21*(MF))))*(1/3)'
                'sync_tp_2_3=(1/(1+exp(-b23*(FF)-c23*(MF))))*(1/3)'
                'sync_tp_3_1=(1/(1+exp(-b31*(FF)-c31*(MF))))*(1/3)'
                'sync_tp_3_2=(1/(1+exp(-b32*(FF)-c32*(MF))))*(1/3)'
                }};        
    elseif normal==0
        prob_fct = {{
                'sync_tp_1_2=(1/(1+exp(+b12*(FF)-c12*(MF))))*(1/3)'
                'sync_tp_1_3=(1/(1+exp(+b13*(FF)-c13*(MF))))*(1/3)'
                'sync_tp_2_1=(1/(1+exp(-b21*(FF)+c21*(MF))))*(1/3)'
                'sync_tp_2_3=(1/(1+exp(+b23*(FF)-c23*(MF))))*(1/3)'
                'sync_tp_3_1=(1/(1+exp(-b31*(FF)+c31*(MF))))*(1/3)'
                'sync_tp_3_2=(1/(1+exp(-b32*(FF)+c32*(MF))))*(1/3)'
                }};        
    end
    prob_params = {{'b12','c12','b13','c13','b21','c21','b23','c23','b31','c31','b32','c32'}};
end


markov_chains=struct('name','sync',...
    'number_of_states',3,...
    'controlled_parameters',{switches},...
    'endogenous_probabilities',prob_fct,...
    'probability_parameters',prob_params);
   


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
% var_prior.L1=1;
% var_prior.L5=1; 
% var_prior.L6=0; % cointegration
% var_prior.L7=0; % unit root
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both

% priors for the sync transition probabilities
%------------------------------------------------
switch_prior=struct();
if const==1
    if normal==1
        switch_prior.a12={0,0,2,'normal'};
        switch_prior.a13={0,0,2,'normal'};
        switch_prior.a21={0,0,2,'normal'};
        switch_prior.a23={0,0,2,'normal'};
        switch_prior.a31={0,0,2,'normal'};
        switch_prior.a32={0,0,2,'normal'};
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
    elseif normal==0
        switch_prior.a12={0.5,0.5,0.5,'normal'};
        switch_prior.a13={0.5,0.5,0.5,'normal'};
        switch_prior.a21={0.5,0.5,0.5,'normal'};
        switch_prior.a23={0.5,0.5,0.5,'normal'};
        switch_prior.a31={0.5,0.5,0.5,'normal'};
        switch_prior.a32={0.5,0.5,0.5,'normal'};
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b13={0.5,0.5,0.25,'gamma'};
        switch_prior.c13={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
        switch_prior.b23={0.5,0.5,0.25,'gamma'};
        switch_prior.c23={0.5,0.5,0.25,'gamma'};
        switch_prior.b31={0.5,0.5,0.25,'gamma'};
        switch_prior.c31={0.5,0.5,0.25,'gamma'};
        switch_prior.b32={0.5,0.5,0.25,'gamma'};
        switch_prior.c32={0.5,0.5,0.25,'gamma'};
    end
elseif const==0
    if normal==1
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
%         switch_prior.b12={0.25,0.25,0.1,'normal'};
%         switch_prior.c12={-0.25,-0.25,0.1,'normal'};
%         switch_prior.b13={0.25,0.25,0.1,'normal'};
%         switch_prior.c13={-0.25,-0.25,0.1,'normal'};
%         switch_prior.b21={-0.25,0.25,0.1,'normal'};
%         switch_prior.c21={0.25,0.25,0.1,'normal'};
%         switch_prior.b23={0.25,0.25,0.1,'normal'};
%         switch_prior.c23={-0.25,-0.25,0.1,'normal'};
%         switch_prior.b31={-0.25,-0.25,0.1,'normal'};
%         switch_prior.c31={0.25,0.25,0.1,'normal'};
%         switch_prior.b32={-0.25,-0.25,0.1,'normal'};
%         switch_prior.c32={0.25,0.25,0.1,'normal'};
    elseif normal==0
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b13={0.5,0.5,0.25,'gamma'};
        switch_prior.c13={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
        switch_prior.b23={0.5,0.5,0.25,'gamma'};
        switch_prior.c23={0.5,0.5,0.25,'gamma'};
        switch_prior.b31={0.5,0.5,0.25,'gamma'};
        switch_prior.c31={0.5,0.5,0.25,'gamma'};
        switch_prior.b32={0.5,0.5,0.25,'gamma'};
        switch_prior.c32={0.5,0.5,0.25,'gamma'};
    end
end


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
[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;    
p_reg3 = f.smoothed_regime_probabilities.regime3.data;    

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

y_fit = fit1.*p_reg1 + fit2.*p_reg2 + fit3.*p_reg3;

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

figSize2 = [10 12];
numticks = 48;
left_color = [0 0 0];

fig=figure;

% Good Regime
right_color = colors(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,1,[0.125,0.1]);
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
plot(dates(sd:ed),p_reg3(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.005 0.35 0.3],'String','\textbf{Good Regime}','interpreter','Latex',titleSettings{:})

% Normal Regime
right_color = colors(1,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,2,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg2(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.06 0.35 0.3],'String','\textbf{Normal Regime}','interpreter','Latex',titleSettings{:})

% Bad Regime
right_color = colors(2,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,3,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg3_LB(sd:ed) p_reg3_M(sd:ed) p_reg3_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg1(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.12 0.35 0.3],'String','\textbf{Bad Regime}','interpreter','Latex',titleSettings{:})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize2);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize2(1)/2 figSize2(2)/30 figSize2(1) figSize2(2)]);
tightfigadv;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end


%% Regime Probabilities Against Data & RHS

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close


figSize2 = [10 12];
numticks = 48;
left_color = [0 0 0];

fig=figure;

% Good Regime
right_color = colors(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,1,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
plot(dates(sd:ed),db.MF.data(sd:ed),'k:','Color',colors(2,:),'LineWidth', 1.5)
%plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg1_LB(sd:ed) p_reg1_M(sd:ed) p_reg1_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg3(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.005 0.35 0.3],'String','\textbf{Good Regime}','interpreter','Latex',titleSettings{:})

% Normal Regime
right_color = colors(1,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,2,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
plot(dates(sd:ed),db.MF.data(sd:ed),'k:','Color',colors(2,:),'LineWidth', 1.5)
%plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg2(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.06 0.35 0.3],'String','\textbf{Normal Regime}','interpreter','Latex',titleSettings{:})

% Bad Regime
right_color = colors(2,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
sfig=subplot_tight(3,1,3,[0.125,0.1]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
plot(dates(sd:ed),db.MF.data(sd:ed),'k:','Color',colors(2,:),'LineWidth', 1.5)
%plot(dates(sd:ed),db.GDPGH.data(sd:ed)+db.TRENDH.data(sd:ed),'k-.','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
hold off
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg3_LB(sd:ed) p_reg3_M(sd:ed) p_reg3_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg1(sd:ed),'Linewidth',1.5)
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
annotation('textbox','Position',[spPos(1,1)+0.25 spPos(1,2)-0.12 0.35 0.3],'String','\textbf{Bad Regime}','interpreter','Latex',titleSettings{:})
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize2);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize2(1)/2 figSize2(2)/15 figSize2(1) figSize2(2)]);
tightfigadv;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs_X'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs_X']));
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
plot(dates(sd:ed), fit2(sd:ed),'Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed), fit3(sd:ed),'Color',colors(3,:),'LineWidth', 1.5)
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
set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Quantiles_MSFitted_QR'],'-bestfit');
   %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
end


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

myList={'sync_tp_1_2','sync_tp_2_1'};

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


%% Ergodic Distribution

params=[results.pop.x];
myparams=a2tilde_to_a(params);

T = length(db.GDPGH.data)-nlags;
p_reg1 = zeros(T,options.N);
p_reg2 = zeros(T,options.N);
p_reg3 = zeros(T,options.N);
y_erg = zeros(T,options.N);

parfor i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
    
    [Resids,Fits]=residuals(sv,myparams(:,i));
    [~,~,~,f]=filter(sv,myparams(:,i));
    
    fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
    fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
    fit3 = Fits.GDPGH.data(:,3)+db.TRENDH.data(nlags+1:end);
    resid1 = Resids.shock_1.data(:,1);
    resid2 = Resids.shock_1.data(:,2);
    resid3 = Resids.shock_1.data(:,3);
    y1 = fit1 + sqrt(myparams(end-2,i))*randn; 
    y2 = fit2 + sqrt(myparams(end-1,i))*randn; 
    y3 = fit3 + sqrt(myparams(end,i))*randn; 
    
    p_reg1(:,i) = f.smoothed_regime_probabilities.regime1.data;
    p_reg2(:,i) = f.smoothed_regime_probabilities.regime2.data;    
    p_reg3(:,i) = f.smoothed_regime_probabilities.regime3.data;    
    y_erg(:,i) = y1.*p_reg1(:,i) + y2.*p_reg2(:,i) + y3.*p_reg3(:,i);
    
    
end

inputformat = 'yyyy-MMM';
dates = datenum((datetime('1973-Jan','InputFormat',inputformat)):calmonths(1):(datetime('2019-May','InputFormat',inputformat)))';
dataformat  = 'yyyy-mmm';

p_reg1_LB = prctile(p_reg1,5,2);
p_reg1_M = prctile(p_reg1,50,2);
p_reg1_UB = prctile(p_reg1,95,2);

p_reg2_LB = prctile(p_reg2,5,2);
p_reg2_M = prctile(p_reg2,50,2);
p_reg2_UB = prctile(p_reg2,95,2);

p_reg3_LB = prctile(p_reg3,5,2);
p_reg3_M = prctile(p_reg3,50,2);
p_reg3_UB = prctile(p_reg3,95,2);

q10 = prctile(y_erg,10,2);
q50 = prctile(y_erg,50,2);
q90 = prctile(y_erg,90,2);



%% Density for Specific Period

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


%% Comparison with Quantile Regression Densities

colors = [0    0.4470    0.7410;    
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;    
        0.4940    0.1840    0.5560;...
        0.4660 0.6740 0.1888];
    
dcolor = [colors(5,:);colors(2,:);colors(1,:)];

% good = regime 1
% bad = regime 3
% normal = regime 2

if exist('ResMatch')==0
    load('/mcr/data_bfi/m1fxl02/3. Miscellanea/Coronavirus VIX/Text_Paper/Figures_Monthly/US/horizon_12/results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

% Select Time Periods for Which to Plot PDFs
% pick representative for bad, normal and good
%period_pdf = {'1997-Mar','2008-Mar','2018-Mar'}; 
period_pdf = {'2005-Mar','2008-Apr','2011-Dec'}; 
for tt=1:length(period_pdf)
    tind(tt) = find(datenum(char(period_pdf(tt)),dataformat)==dates);
end

clc;
p_reg1(tind(1))
p_reg3(tind(2))
p_reg2(tind(3))

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
disp('bad = regime 3')
disp('normal = regime 2')
disp('good = regime 1')

pmode

params_CI = [params_pmode-(params_M_Struct-params_LB_Struct) params_pmode params_pmode+(params_UB_Struct-params_M_Struct)];

disp('Transition Probabilities')
round(params_CI(20:31,:),2)
disp('Coefficients')
round([params_CI(1:19,:); params_CI(32:end-2,:)],2)
disp('Volatilities')
round((params_CI(end-2:end,:)),2) % no need to take sqrt in structural form

diary off
