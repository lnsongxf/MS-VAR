%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all;
% clc()

saveit    = 1; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

alternative = 1; % =1 if alternative test 

% Model label
if alternative ==0
    model_label = ['_const_' num2str(const) '_normal_20201108_' num2str(normal)];
elseif alternative ==1
    model_label = ['_const_' num2str(const) '_normal_20201108_alternative_' num2str(normal)];
end


% Important paths
% addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));
addpath(genpath('RISE_toolbox'));


%% Load RISE
rise_startup()

%% load the data
if alternative ==0
    load data_20201108
elseif alternative ==1
    load data_20201108_alternative
end


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
model = 'SVAR_StateDependent_2Regimes_SepMarkov';

% Folder for storing results
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

slides_folder = ['Slides/FiguresPhillyFed/'];
if exist(slides_folder,'dir')==0
    mkdir(slides_folder)
end


if const==1
    if normal==1
        prob_fct = {{

        % Format is [name_tp_1_2]
        'sync_tp_1_2=1/(1+exp(a12-b12*(FF)-c12*(MF)))'
        'sync_tp_2_1=1/(1+exp(a21-b21*(FF)-c21*(MF)))'
        }};
    
    prob_fct_vol = {{

        % Format is [name_tp_1_2]
        'vol_tp_1_2=1/(1+exp(a12_vol-b12_vol*(FF)-c12_vol*(MF)))'
        'vol_tp_2_1=1/(1+exp(a21_vol-b21_vol*(FF)-c21_vol*(MF)))'
        }};
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'
            }};
        prob_fct_vol = {{
            'vol_tp_1_2=1/(1+exp(a12_vol-b12_vol*(FF)+c12_vol*(MF)))'
            'vol_tp_2_1=1/(1+exp(a21_vol+b21_vol*(FF)-c21_vol*(MF)))'
            }};
    end
    prob_params = {{'a12','a21','b12','c12','b21','c21'}};
    prob_params_vol = {{'a12_vol','a21_vol','b12_vol','c12_vol','b21_vol','c21_vol'}};
elseif const==0
    if normal==1
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)-c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(-b21*(FF)-c21*(MF)))'
            }};
        prob_fct_vol = {{
            'vol_tp_1_2=1/(1+exp(-b12_vol*(FF)-c12_vol*(MF)))'
            'vol_tp_2_1=1/(1+exp(-b21_vol*(FF)-c21_vol*(MF)))'
            }};
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(+b21*(FF)-c21*(MF)))'
            }};
        prob_fct_vol = {{
            'vol_tp_1_2=1/(1+exp(-b12_vol*(FF)+c12_vol*(MF)))'
            'vol_tp_2_1=1/(1+exp(+b21_vol*(FF)-c21_vol*(MF)))'
            }};        
    end
    prob_params = {{'b12','c12','b21','c21'}};
    prob_params_vol = {{'b12_vol','c12_vol','b21_vol','c21_vol'}};
end

markov_chains=struct('name','sync',...
    'number_of_states',2,...
    'controlled_parameters',{switches_coef},...   % these correspond to the parameters that are switching
    'endogenous_probabilities',prob_fct,...
    'probability_parameters',prob_params);

markov_chains(2)=struct('name','vol',...
    'number_of_states',2,...
    'controlled_parameters',{switches_vol},...   % these correspond to the parameters that are switching
    'endogenous_probabilities',prob_fct_vol,...
    'probability_parameters',prob_params_vol);

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

% priors for the sync transition probabilities
%------------------------------------------------
switch_prior=struct();
if const==1
    if normal==1
        %         coefficients
        switch_prior.a12={0,0,2,'normal'};
        switch_prior.b12={0,0,2,'normal'};
        switch_prior.c12={0,0,2,'normal'};
        switch_prior.a21={0,0,2,'normal'};
        switch_prior.b21={0,0,2,'normal'};
        switch_prior.c21={0,0,2,'normal'};
        %         volatilities
        switch_prior.a12_vol={0,0,2,'normal'};
        switch_prior.b12_vol={0,0,2,'normal'};
        switch_prior.c12_vol={0,0,2,'normal'};
        switch_prior.a21_vol={0,0,2,'normal'};
        switch_prior.b21_vol={0,0,2,'normal'};
        switch_prior.c21_vol={0,0,2,'normal'};
    elseif normal==0
        %         coefficients
        switch_prior.a12={0.5,0.5,0.5,'normal'};
        switch_prior.a21={0.5,0.5,0.5,'normal'};
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
        %         volatilities
        switch_prior.a12_vol={0.5,0.5,0.5,'normal'};
        switch_prior.a21_vol={0.5,0.5,0.5,'normal'};
        switch_prior.b12_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.c12_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.b21_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.c21_vol={0.5,0.5,0.25,'gamma'};
    end
elseif const==0
    if normal==1
        %         coefficients
        switch_prior.b12={-0.25,-0.25,0.1,'normal'};
        switch_prior.c12={0.25,0.25,0.1,'normal'};
        switch_prior.b21={0.25,0.25,0.1,'normal'};
        switch_prior.c21={-0.25,-0.25,0.1,'normal'};
        %         volatilities
        switch_prior.b12_vol={-0.25,-0.25,0.1,'normal'};
        switch_prior.c12_vol={0.25,0.25,0.1,'normal'};
        switch_prior.b21_vol={0.25,0.25,0.1,'normal'};
        switch_prior.c21_vol={-0.25,-0.25,0.1,'normal'};
    elseif normal==0
        %         coefficients
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
        %         volatilities
        switch_prior.b12_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.c12_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.b21_vol={0.5,0.5,0.25,'gamma'};
        switch_prior.c21_vol={0.5,0.5,0.25,'gamma'};
    end
end

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

%% Find posterior mode

%options = optimset('TolX',1e-12);
if exist([mcmc_folder 'sv' model_label '.mat'],'file')==2
    load([mcmc_folder 'sv' model_label '.mat'])
else
    %sv=estimate(sv0,db,{'1973M1','2019M5'},prior,restrictions,'fmincon',false,'optimset',options);
    sv=estimate(sv0,db,{'1973M1','2019M10'},prior,restrictions,'fmincon');
    if save_mcmc==1
        save([mcmc_folder 'sv' model_label '.mat'],'sv')
    end
end

%% Store estimates

pmode=posterior_mode(sv)

%% Printing estimates

dfile = ([log_folder 'structural_' model_label '.txt']);
if exist(dfile, 'file')
    delete(dfile);
end
diary(dfile)
print_structural_form(sv)
diary off


%% Printing solution

dfile = ([log_folder 'reduced' model_label '.txt']);
if exist(dfile, 'file')
    delete(dfile);
end
diary(dfile)
print_solution(sv)
diary off

%% plots probabilities against data

% %{
%plot_probabilities(sv)
%plot_data_against_probabilities(sv,'regime')
%}


%% Compute Ergodic Fitted Values

[Resids,Fits]=residuals(sv);
fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
fit3 = Fits.GDPGH.data(:,3)+db.TRENDH.data(nlags+1:end);
fit4 = Fits.GDPGH.data(:,4)+db.TRENDH.data(nlags+1:end);


% Smoothed and Filtered probabilities
[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;
p_reg3 = f.smoothed_regime_probabilities.regime3.data;
p_reg4 = f.smoothed_regime_probabilities.regime4.data;

p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
p_reg3_filtered = f.filtered_regime_probabilities.regime3.data;
p_reg4_filtered = f.filtered_regime_probabilities.regime4.data;

%----------------------------------------------------
% Figure: Fitted value and realization of dependent variable
%----------------------------------------------------

inputformat = 'yyyy-MMM';
dates = datenum((datetime('1973-Jan','InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime('2019-Oct','InputFormat',inputformat)))';
dataformat  = 'yyyy-mmm';

FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};

start_plot = '1973-Feb';
end_plot   = '2019-Oct';
sd = find(datenum(start_plot,dataformat)==dates);
ed = find(datenum(end_plot,dataformat)==dates);

% y_fit = fit1.*p_reg1 + fit2.*p_reg2;
y_fit = fit1.*p_reg1 + fit2.*p_reg2 + fit3.*p_reg3 + fit4.*p_reg4;

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
title('Fitted and Realized GDP','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Ergodic_Fitted'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Ergodic_Fitted']));
end

%%
%----------------------------------------------------
% Figure: Realized innovations to GDP equation
%----------------------------------------------------

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
title('Realized Innovation to GDP Equation','FontSize',16','Interpreter','Latex');
tightfig;
if alternative ==0
    load data_20201108
elseif alternative ==1
    load data_20201108_alternative
end
if saveit==1
    print('-dpdf',fig,[fig_folder 'Ergodic_Difference'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Ergodic_Difference']));
end

disp(round(sum((realized-y_fit).^2),2))

%%

%----------------------------------------------------
% Figure: Regime Probabilities Against Data and Factors
%----------------------------------------------------

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close

figSize2 = [10 12];
numticks = 48;
left_color = [0 0 0];
right_color = colors(2,:);

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),realized(sd:ed),'k','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg2(sd:ed),'LineWidth', 1.5)
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
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
title('Smoothed Regime 2 Probabilities and Realized GDP','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs2'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close

figSize2 = [10 12];
numticks = 48;
left_color = [0 0 0];
right_color = colors(2,:);

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),realized(sd:ed),'k','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg3(sd:ed),'LineWidth', 1.5)
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
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
title('Smoothed Regime 3 Probabilities and Realized GDP','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs3'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close

figSize2 = [10 12];
numticks = 48;
left_color = [0 0 0];
right_color = colors(2,:);

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
plot(dates(sd:ed),realized(sd:ed),'k','LineWidth', 1.5)
ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
yyaxis right
%plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
plot(dates(sd:ed),p_reg4(sd:ed),'LineWidth', 1.5)
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
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
title('Smoothed Regime 4 Probabilities and Realized GDP','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs4'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end

% %%
% 
% %----------------------------------------------------
% % Figure: Regime Probabilities Against Data and Factors
% %----------------------------------------------------
% figure;
% set(gca,'ColorOrder','factory');
% colr = get(gca,'ColorOrder');
% colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
% close
% 
% numticks = 48;
% left_color = [0 0 0];
% right_color = colors(2,:);
% 
% fig=figure;
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% hold on
% yyaxis left
% hold on
% plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
% plot(dates(sd:ed),db.MF.data(sd:ed),'c:','LineWidth', 1.5)
% % plot(dates(sd:ed),realized(sd:ed),'k','LineWidth', 1.5)
% ylabel('$\bar{\Delta} y_{t+1,t+12}$','interpreter','Latex','fontsize',10)
% yyaxis right
% %plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
% plot(dates(sd:ed),p_reg2(sd:ed),'LineWidth', 1.5)
% ylabel('Probability','interpreter','Latex','fontsize',10)
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% axis tight
% rr=recessionplot;
% legend('$f_t$','$m_t$','interpreter','Latex','Orientation','Horizontal')
% legend boxoff
% axis tight
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
% tightfig;
% if saveit==1
%     print('-dpdf',fig,[fig_folder 'RegimeProbs_X'],'-bestfit');
%     %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs_X']));
% end



% %% Time Series of Quantiles - QR Comparison
% %----------------------------------------------------
% % Figure: QR quantiles vs MS regime fitted values
% %----------------------------------------------------
% 
% figure;
% set(gca,'ColorOrder','factory');
% colr = get(gca,'ColorOrder');
% colors = [colr(2,:);colr(1,:);colr(5,:)];
% close
% 
% load(['quantiles_QR']);
% Qplot = [0.25 0.5 0.75];
% quantiles_dist = 0.1:0.05:0.9; % quantiles to construct distribution
% indq = ismember(round(quantiles_dist,2),Qplot);
% indq = find(indq==1);
% 
% fig=figure;
% %yyaxis left
% hold on
% plot(dates(sd:ed), fit2(sd:ed),'Color',colors(1,:),'LineWidth', 1.5)
% plot(dates(sd:ed), fit1(sd:ed),'Color',colors(2,:),'LineWidth', 1.5)
% j=0;
% for iquantile=indq
%     j=j+1;
%     plot(dates(sd:ed),YQ_dist_fut(sd:ed,iquantile),'-.','Color',colors(j,:),'LineWidth',1.5);
%     hold on
% end
% hold off
% ylabel('Percent','interpreter','Latex','fontsize',10)
% leg = {'Bad Regime','Normal Regime',[num2str(Qplot(1)*100) 'th Quantile'],'Median',[num2str(Qplot(3)*100) 'th Quantile']};
% axis tight
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% legend(leg,'Orientation','Vertical','Location','South','interpreter','Latex')
% legend boxoff
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% tightfig;
% if saveit==1
%     print('-dpdf',fig,[fig_folder 'Quantiles_MSFitted_QR'],'-bestfit');
%     %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
% end


%% Recursive updating of s(t) from simulations of Markov-chain

% Number of draws for simulting predictive density
nDraws = 20000;

% Set a 0.9 threshold to switch from good to bad regime
prob_reg34_new = zeros(1,length(p_reg2));
prob_reg24_new = zeros(1,length(p_reg2));
prob_reg34_new((p_reg3(sd:ed)+p_reg4(sd:ed))>0.9) = 1; % prob of being in low growth
prob_reg24_new((p_reg2(sd:ed)+p_reg4(sd:ed))>0.9) = 1; % prob of being in high variance



% Draw shocks for simulation
rng(123);
shocks_sim = randn(length(prob_reg34_new),nDraws);

% Collect MF, FF and trend
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);
TR = db.TRENDH.data(nlags+1:end);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws

    for tt= 1:length(prob_reg34_new)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

        if tt<12
            st_lag = 1; % initialize in bad regime (Depends on when the sample starts)
            st_lag_vol = 1; % initialize in high vol regime (Depends on when the sample starts)

            % Set st
            st_use = 1;
            st_use_vol = 1;
            

        else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag = prob_reg34_new(tt-11);
            st_lag_vol =  prob_reg24_new(tt-11);

            % Get transition probabilities from t-11:tt
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF(tt-11:tt))-pmode.c12*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF(tt-11:tt))-pmode.c21*(MF(tt-11:tt))));
                p12_vol = 1./(1+exp(pmode.a12_vol-pmode.b12_vol*(FF(tt-11:tt))-pmode.c12_vol*(MF(tt-11:tt))));
                p21_vol = 1./(1+exp(pmode.a21_vol-pmode.b21_vol*(FF(tt-11:tt))-pmode.c21_vol*(MF(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
                p12_vol = 1./(1+exp(pmode.a12_vol-pmode.b12_vol.*(FF(tt-11:tt))+pmode.c12_vol.*(MF(tt-11:tt))));
                p21_vol = 1./(1+exp(pmode.a21_vol+pmode.b21_vol.*(FF(tt-11:tt))-pmode.c21_vol.*(MF(tt-11:tt))));
            end

            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;
            p11_vol = ones(12,1) - p12_vol;
            p22_vol = ones(12,1) - p21_vol;            


            % Draw regimes from t-12:t 
            for tt2 = 1:12

                % Draw the Markov Chain for period t-12:t
                udraw = rand(1);
                udraw_vol = rand(1);

                if st_lag == 0 && st_lag_vol == 0 % started in high growth / low vol

                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol > p11_vol(tt2)
                        st_vol(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag == 0 && st_lag_vol == 1 % started in high growth / low vol
                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol > p22_vol(tt2)
                        st_vol(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol(tt2,1) = 1;  % don't switch and remain in high vol
                    end   
               elseif st_lag == 1 && st_lag_vol == 0 % started in low growth / low vol
                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol > p11_vol(tt2)
                        st_vol(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag == 1 && st_lag_vol == 1 % started in low growth / high vol
                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol > p22_vol(tt2)
                        st_vol(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol(tt2,1) = 1;  % don't switch and remain in high vol
                    end
               end

                st_lag = st(tt2,1);
                st_lag_vol = st_vol(tt2,1);
            end

            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use = st(tt2);
            % st = 0 (low vol), st = 1 (high vol)
            st_use_vol = st_vol(tt2);
        end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        IND_good_vol = 1 - st_use_vol;
        
        IND_hg_lv = IND_good*IND_good_vol;
        IND_hg_hv = IND_good*st_use_vol;
        IND_lg_lv = st_use*IND_good_vol;
        IND_lg_hv = st_use*st_use_vol;


        % Good regime in period t / low vol
        dY_sim_1(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF(tt)...
            - pmode.a0_3_2_sync_1*MF(tt) + pmode.s_3_3_vol_1*shocks_sim(tt,nsim))...
            + TR(tt);

        % Good regime in period t / high vol
        dY_sim_2(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF(tt)...
            - pmode.a0_3_2_sync_1*MF(tt)+ pmode.s_3_3_vol_2*shocks_sim(tt,nsim))...
            + TR(tt);
        
         % Bad regime in period t / low vol       
        dY_sim_3(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF(tt)...
            - pmode.a0_3_2_sync_2*MF(tt) + pmode.s_3_3_vol_1*shocks_sim(tt,nsim))...
            + TR(tt);

        % Bad regime in period t / high vol
        dY_sim_4(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF(tt)...
            - pmode.a0_3_2_sync_2*MF(tt)+ pmode.s_3_3_vol_2*shocks_sim(tt,nsim))...
            + TR(tt);
        
        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_hg_lv + dY_sim_2(tt,nsim)*IND_hg_hv + dY_sim_3(tt,nsim)*IND_lg_lv + dY_sim_4(tt,nsim)*IND_lg_hv;

        St_sim1(tt,nsim) = IND_hg_lv;
        St_sim2(tt,nsim) = IND_hg_hv;
        St_sim3(tt,nsim) = IND_lg_lv;
        St_sim4(tt,nsim) = IND_lg_hv;
        
        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt) + dY_sim_3(tt,nsim)*p_reg3_filtered(tt) + dY_sim_4(tt,nsim)*p_reg4_filtered(tt);

    end

    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end

% Ergodic probability of regime 1;
p_reg1_sim(:,1) =  sum(St_sim1,2)/size(sum(St_sim1),2);
p_reg2_sim(:,1) =  sum(St_sim2,2)/size(sum(St_sim2),2);
p_reg3_sim(:,1) =  sum(St_sim3,2)/size(sum(St_sim3),2);
p_reg4_sim(:,1) =  sum(St_sim4,2)/size(sum(St_sim4),2);


% Compute percentiles
dY_sim_25 = prctile(dY_sim',25)';
dY_sim_75 = prctile(dY_sim',75)';
dY_sim_10 = prctile(dY_sim',10)';
dY_sim_90 = prctile(dY_sim',90)';

% Compute Regime specific means
dY_sim1_mean = mean(dY_sim_1,2);
dY_sim2_mean = mean(dY_sim_2,2);



%% Recursive updating of s(t) from simulations of Markov-chain (FULL SET)

if alternative ==0
    load data_full_20201108
elseif alternative ==1
    load data_full_20201108_alternative
end

dates_full = datenum((datetime('1973-Jan','InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime('2020-Oct','InputFormat',inputformat)))';
end_plot_full   = '2020-Oct';
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Number of draws for simulting predictive density
nDraws = 20000;

% Set a 0.9 threshold to switch from good to bad regime
prob_reg34_new_full = zeros(1,length(dates_full));
prob_reg24_new_full = zeros(1,length(dates_full));
prob_reg34_new_full((p_reg3(sd:ed)+p_reg4(sd:ed))>0.9) = 1; % prob of being in low growth
prob_reg24_new_full((p_reg2(sd:ed)+p_reg4(sd:ed))>0.9) = 1; % prob of being in high variance



% Draw shocks for simulation
rng(123);
shocks_sim_full = randn(length(prob_reg34_new_full),nDraws);

% Collect MF, FF and trend
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);
TR_full = db_full.TRENDH.data(nlags+1:end);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws

    for tt= 1:length(prob_reg34_new_full)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

        if tt<12
            st_lag_full = 1; % initialize in bad regime (Depends on when the sample starts)
            st_lag_vol_full = 1; % initialize in high vol regime (Depends on when the sample starts)

            % Set st
            st_use_full = 1;
            st_use_vol_full = 1;
            
        else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag_full = prob_reg34_new_full(tt-11);
            st_lag_vol_full =  prob_reg24_new_full(tt-11);

            % Get transition probabilities from t-11:tt
            if normal
                % Transition probabilities when coefficients have normal prior
                p12_full = 1./(1+exp(pmode.a12-pmode.b12*(FF_full(tt-11:tt))-pmode.c12*(MF_full(tt-11:tt))));
                p21_full = 1./(1+exp(pmode.a21-pmode.b21*(FF_full(tt-11:tt))-pmode.c21*(MF_full(tt-11:tt))));
                p12_vol_full = 1./(1+exp(pmode.a12_vol-pmode.b12_vol*(FF_full(tt-11:tt))-pmode.c12_vol*(MF_full(tt-11:tt))));
                p21_vol_full = 1./(1+exp(pmode.a21_vol-pmode.b21_vol*(FF_full(tt-11:tt))-pmode.c21_vol*(MF_full(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12_full = 1./(1+exp(pmode.a12-pmode.b12.*(FF_full(tt-11:tt))+pmode.c12.*(MF_full(tt-11:tt))));
                p21_full = 1./(1+exp(pmode.a21+pmode.b21.*(FF_full(tt-11:tt))-pmode.c21.*(MF_full(tt-11:tt))));
                p12_vol_full = 1./(1+exp(pmode.a12_vol-pmode.b12_vol.*(FF_full(tt-11:tt))+pmode.c12_vol.*(MF_full(tt-11:tt))));
                p21_vol_full = 1./(1+exp(pmode.a21_vol+pmode.b21_vol.*(FF_full(tt-11:tt))-pmode.c21_vol.*(MF_full(tt-11:tt))));
            end

            % Compute the probability of remaining in a given regime
            p11_full = ones(12,1) - p12_full;
            p22_full = ones(12,1) - p21_full;
            p11_vol_full = ones(12,1) - p12_vol_full;
            p22_vol_full = ones(12,1) - p21_vol_full;            


            % Draw regimes from t-12:t 
            for tt2 = 1:12

                % Draw the Markov Chain for period t-12:t
                udraw_full = rand(1);
                udraw_vol_full = rand(1);

                if st_lag_full == 0 && st_lag_vol_full == 0 % started in high growth / low vol

                    if udraw_full > p11_full(tt2)
                        st_full(tt2,1) = 1; % switch from good to bad
                    else
                        st_full(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol_full > p11_vol_full(tt2)
                        st_vol_full(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol_full(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag_full == 0 && st_lag_vol_full == 1 % started in high growth / low vol
                    if udraw_full > p11_full(tt2)
                        st_full(tt2,1) = 1; % switch from good to bad
                    else
                        st_full(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol_full > p22_vol_full(tt2)
                        st_vol_full(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol_full(tt2,1) = 1;  % don't switch and remain in high vol
                    end   
               elseif st_lag_full == 1 && st_lag_vol_full == 0 % started in low growth / low vol
                    if udraw_full > p22_full(tt2)
                        st_full(tt2,1) = 0; % switch from bad to good
                    else
                        st_full(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol_full > p11_vol_full(tt2)
                        st_vol_full(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol_full(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag_full == 1 && st_lag_vol_full == 1 % started in low growth / high vol
                    if udraw_full > p22_full(tt2)
                        st_full(tt2,1) = 0; % switch from bad to good
                    else
                        st_full(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol_full > p22_vol_full(tt2)
                        st_vol_full(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol_full(tt2,1) = 1;  % don't switch and remain in high vol
                    end
               end

                st_lag_full = st_full(tt2,1);
                st_lag_vol_full = st_vol_full(tt2,1);
            end

            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use_full = st_full(tt2);
            % st = 0 (low vol), st = 1 (high vol)
            st_use_vol_full = st_vol_full(tt2);
        end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good_full = 1 - st_use_full;
        IND_good_vol_full = 1 - st_use_vol_full;
        
        IND_hg_lv_full = IND_good_full*IND_good_vol_full;
        IND_hg_hv_full = IND_good_full*st_use_vol_full;
        IND_lg_lv_full = st_use_full*IND_good_vol_full;
        IND_lg_hv_full = st_use_full*st_use_vol_full;


        % Good regime in period t / low vol
        dY_sim_1_full(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF_full(tt)...
            - pmode.a0_3_2_sync_1*MF_full(tt) + pmode.s_3_3_vol_1*shocks_sim_full(tt,nsim))...
            + TR_full(tt);

        % Good regime in period t / high vol
        dY_sim_2_full(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF_full(tt)...
            - pmode.a0_3_2_sync_1*MF_full(tt)+ pmode.s_3_3_vol_2*shocks_sim_full(tt,nsim))...
            + TR_full(tt);
        
         % Bad regime in period t / low vol       
        dY_sim_3_full(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF_full(tt)...
            - pmode.a0_3_2_sync_2*MF_full(tt) + pmode.s_3_3_vol_1*shocks_sim_full(tt,nsim))...
            + TR_full(tt);

        % Bad regime in period t / high vol
        dY_sim_4_full(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF_full(tt)...
            - pmode.a0_3_2_sync_2*MF_full(tt)+ pmode.s_3_3_vol_2*shocks_sim_full(tt,nsim))...
            + TR_full(tt);
        
        dY_sim_full(tt,nsim) = dY_sim_1_full(tt,nsim)*IND_hg_lv_full + dY_sim_2_full(tt,nsim)*IND_hg_hv_full + dY_sim_3_full(tt,nsim)*IND_lg_lv_full + dY_sim_4_full(tt,nsim)*IND_lg_hv_full;

        St_sim1_full(tt,nsim) = IND_hg_lv_full;
        St_sim2_full(tt,nsim) = IND_hg_hv_full;
        St_sim3_full(tt,nsim) = IND_lg_lv_full;
        St_sim4_full(tt,nsim) = IND_lg_hv_full;
        
        % Weight regime simulation using filtered probabilities
%         dY_sim_weighted(tt,nsim) = dY_sim_1_full(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2_full(tt,nsim)*p_reg2_filtered(tt) + dY_sim_3_full(tt,nsim)*p_reg3_filtered(tt) + dY_sim_4_full(tt,nsim)*p_reg4_filtered(tt);

    end

    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end

% Ergodic probability of regime 1;
p_reg1_sim_full(:,1) =  sum(St_sim1_full,2)/size(sum(St_sim1_full),2);
p_reg2_sim_full(:,1) =  sum(St_sim2_full,2)/size(sum(St_sim2_full),2);
p_reg3_sim_full(:,1) =  sum(St_sim3_full,2)/size(sum(St_sim3_full),2);
p_reg4_sim_full(:,1) =  sum(St_sim4_full,2)/size(sum(St_sim4_full),2);


% Compute percentiles
dY_sim_25_full = prctile(dY_sim_full',25)';
dY_sim_75_full = prctile(dY_sim_full',75)';
dY_sim_10_full = prctile(dY_sim_full',10)';
dY_sim_90_full = prctile(dY_sim_full',90)';

% Compute Regime specific means
dY_sim1_mean_full = mean(dY_sim_1_full,2);
dY_sim2_mean_full = mean(dY_sim_2_full,2);


% %%
% %----------------------------------------------------
% % Figure: Compare Real-Time Estimate of Prob of Bad Regime 
% %----------------------------------------------------
% 
% fig=figure; clf;
% l1=plot(dates(sd:ed),p_reg2,'r-','DisplayName','p(bad regime) filtered'); hold on;
% l2=plot(dates(sd:ed),p_reg2_sim,'--','DisplayName','p(bad regime) simulated');
% for row=1:length(regime2_shade)
%  patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [0 0, 1 1], [0.8 0.8 0.8],'EdgeColor','none'); hold on;
% end
% hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
% set(gca,'children',flipud(get(gca,'children')))
% axis tight
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% title('Probability of Bad Regime: Filtered vs Real-Time','FontSize',16','Interpreter','Latex');
% tightfig;

% %%
% %----------------------------------------------------
% % Figure: Compare Real-Time Estimate of Prob of Bad Regime (FULL SET)
% %----------------------------------------------------
% 
% 
% fig=figure; clf;
% l2=plot(dates_full(sd:ed_full),p_reg2_sim_full,'-','LineWidth', 1.5,'DisplayName','p(bad regime) simulated'); hold on;
% for row=1:length(regime2_shade_full)
%  patch([regime2_shade_full(row,1) regime2_shade_full(row,2) regime2_shade_full(row,2) regime2_shade_full(row,1)], [0 0, 1 1], [0.8 0.8 0.8],'EdgeColor','none'); hold on;
% end
% % hleg = legend(l2,'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
% set(gca,'children',flipud(get(gca,'children')))
% axis tight
% ax=gca;
% ax.XTick = datenum(dates_full(sd:numticks:ed_full));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% % title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
% tightfig;
% if saveit==1
%     print('-dpdf',fig,[fig_folder 'Prob_bad_regime'],'-bestfit');
%     %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
% end

% %% ========================================================================
% % PLOT FIGURE 10. 25. 75. 90. percentiles
% % =========================================================================
% 
% % Set colors
% %colors = cbrewer('qual', 'Set2', 8);
% colors = cbrewer('div', 'RdYlBu', 64);
% 
% fig=figure; clf;
% %yyaxis left
% hold on
% % Plot percentiles from MS
% l1=plot(dates(sd:ed), dY_sim_10(sd:ed),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
% l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
% l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
% l4=plot(dates(sd:ed), dY_sim_90(sd:ed),'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
% 
% % plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% % is the 75th Quantile. Other quantiles are in quantiles_dist.
% l5=plot(dates(sd:ed),YQ_dist_fut(sd:ed,4),'--','Color',colors(20,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(1)*100) 'th Quantile']);
% l6=plot(dates(sd:ed),YQ_dist_fut(sd:ed,14),'--','Color',colors(50,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(3)*100) 'th Quantile']);
% % Plot shades of bad regime
% % for row=1:length(regime2_shade)
% % patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [-8 -8, 8 8], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
% % end
% hleg = legend([l1 l2 l3 l4 l5 l6],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
% set(gca,'children',flipud(get(gca,'children')))
% hold off
% ylabel('Percent','interpreter','Latex','fontsize',10)
% axis tight
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% tightfig;

%% Just MS FULL



%==========================================================================
% Quantiles from MS Model: .10, .25, .75, .90
%==========================================================================

% Set colors
colors = cbrewer('div', 'RdYlBu', 64);
numticks = 48;
fig=figure; clf;
hold on
% Plot percentiles from MS
l1=plot(dates_full(sd:ed_full), dY_sim_10_full(sd:ed_full),'Color',colors(5,:),'LineWidth', 3,'DisplayName','10th');
l2=plot(dates_full(sd:ed_full), dY_sim_25_full(sd:ed_full),'--','Color',colors(15,:),'LineWidth', 2.5,'DisplayName','25th');
l3=plot(dates_full(sd:ed_full), dY_sim_75_full(sd:ed_full),'--','Color',colors(55,:),'LineWidth', 2.5,'DisplayName','75th');
l4=plot(dates_full(sd:ed_full), dY_sim_90_full(sd:ed_full),'Color',colors(60,:),'LineWidth', 3,'DisplayName','90th');
rr=recessionplot;
hleg = legend([l1 l2 l3 l4],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('$\bar{\Delta} {y}_{t+1,t+12}$ (\%)','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates_full(ed_full)])
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
print('-dpdf',fig,[slides_folder model '_Quantiles_MS_full'],'-bestfit');

%%

start_plot2 = '2007-Jan';
sd2 = find(datenum(start_plot2,dataformat)==dates);

fig=figure; clf;
%yyaxis left
hold on
% Plot percentiles from MS
l1=plot(dates_full(sd2:ed_full), dY_sim_10_full(sd2:ed_full),'Color',colors(5,:),'LineWidth', 2.5,'DisplayName','10th');
l2=plot(dates_full(sd2:ed_full), dY_sim_25_full(sd2:ed_full),'--','Color',colors(15,:),'LineWidth', 2.5,'DisplayName','25th');
l3=plot(dates_full(sd2:ed_full), dY_sim_75_full(sd2:ed_full),'--','Color',colors(45,:),'LineWidth', 2.5,'DisplayName','75th');
l4=plot(dates_full(sd2:ed_full), dY_sim_90_full(sd2:ed_full),'Color',colors(55,:),'LineWidth', 2.5,'DisplayName','90th');
rr=recessionplot;
hleg = legend([l1 l2 l3 l4],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates_full(sd2:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd2), dates_full(ed_full)])
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
print('-dpdf',fig,[slides_folder model '_Quantiles_MS_recent'],'-bestfit');


% %% ========================================================================
% % PLOT FIGURE 25. 75. percentiles vs Expected value in each regime
% % =========================================================================
% 
% % Set colors
% colors = cbrewer('div', 'RdYlBu', 64);
% colors2 = cbrewer('qual', 'Set1', 8);
% 
% fig=figure; clf;
% %yyaxis left
% hold on
% l1=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
% l2=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
% l3=plot(dates(sd:ed), dY_sim1_mean(sd:ed),'-.','Color',colors2(2,:),'LineWidth', 1.5,'DisplayName','$E(Good Regime)$');
% l4=plot(dates(sd:ed), dY_sim2_mean(sd:ed),'--','Color',colors2(4,:),'LineWidth', 1.5,'DisplayName','$E(Bad Regime)$');
% 
% % Plot shades of bad regime
% for row=1:length(regime2_shade)
% patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [-8 -8, 8 8], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
% end
% 
% % Format plot
% hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
% set(gca,'children',flipud(get(gca,'children')))
% ylabel('Percent','interpreter','Latex','fontsize',10)
% axis tight
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% tightfig;

% %%
% %========================================================================
% % Percentile of the expected value in regime 1 and 2
% %=========================================================================
% for tt=1:length(prob_reg2_new)
%     %[,Ydens] = ksdensity(dY_sim(tt,:));
% 
%     % Fit the period-t predictive distribution
%     pdY = fitdist(dY_sim(tt,:)','Kernel','Kernel','epanechnikov');
% 
%     percentile_reg1(tt,1) = cdf(pdY,dY_sim1_mean(tt,1));
%     percentile_reg2(tt,1) = cdf(pdY,dY_sim2_mean(tt,1));
% end
% 
% 
% 
% colors = cbrewer('div', 'RdYlBu', 64);
% 
% fig=figure;
% %yyaxis left
% hold on
% l1=plot(dates(sd:ed), percentile_reg1(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','prctl Mean of Good regime');
% l2=plot(dates(sd:ed), percentile_reg2(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','prctl Mean of Bad regime');
% l3=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg1(sd:ed)),'--','Color',colors(45,:),'LineWidth', 2);
% l4=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg2(sd:ed)),'--','Color',colors(15,:),'LineWidth', 2);
% % 
% % % Plot shades of bad regime
% % for row=1:length(regime2_shade)
% % patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [0 0, 1 1], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
% % end
% 
% % Format plot
% hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
% set(gca,'children',flipud(get(gca,'children')))
% ylabel('Percent','interpreter','Latex','fontsize',10)
% axis tight
% ax=gca;
% ax.XTick = datenum(dates(sd:numticks:ed));
% datetick('x','yyyy','keepticks')
% set(gca, 'XLim', [dates(sd), dates(ed)])
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% hold off
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
% tightfig;

%%
%==========================================================================
% Posterior sampling
%==========================================================================

pnames=fieldnames(pmode);

a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);

options=struct();
options.alpha=0.234;
options.thin=10;
options.burnin=10^3; %^3
options.N=2*10^4; %^4
options.nchain=1;

lb(~isfinite(lb))=-500;
ub(~isfinite(ub))=500;

if exist([mcmc_folder 'results' model_label '.mat'],'file')==2
    load([mcmc_folder 'results' model_label '.mat'])
else
    if save_mcmc==1
        results=mh_sampler(ff,lb,ub,options,x0,vcov);
        save([mcmc_folder 'results' model_label '.mat'],'results')
    end
end

%==========================================================================
% COMPUTE OBJECTS FROM POSTERIOR DISTRIBUTION
%==========================================================================

% Collect parameters
params=[results.pop.x];

myparams=a2tilde_to_a(params);

p_reg1_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
p_reg2_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
p_reg3_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
p_reg4_mat = zeros(length(db.GDPGH.data)-nlags,options.N);

y_erg = zeros(length(db.GDPGH.data)-nlags,options.N);


nDraws2 = 1;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_2_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_3_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_4_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_weighted_mat = NaN(length(MF),nDraws2,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Simulating model at posterior draws...");

% COMPUTE PREDICTIVE DENSITY
tic
parfor i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_

    paramas_use = myparams(:,i);

    % Compute regime specific densities
    [Resids,Fits]=residuals(sv,paramas_use);
    [~,~,~,f]=filter(sv,paramas_use);

    fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
    fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
    fit3 = Fits.GDPGH.data(:,3)+db.TRENDH.data(nlags+1:end);
    fit4 = Fits.GDPGH.data(:,4)+db.TRENDH.data(nlags+1:end);
    resid1 = Resids.shock_3.data(:,1);
    resid2 = Resids.shock_3.data(:,2);
    resid3 = Resids.shock_3.data(:,3);
    resid4 = Resids.shock_3.data(:,4);

    y1 = fit1 + (paramas_use(end-1))*randn;
    y2 = fit2 + (paramas_use(end))*randn;
    y3 = fit1 + (paramas_use(end-1))*randn;
    y4 = fit2 + (paramas_use(end))*randn;
    

    p_reg1_mat(:,i) = f.smoothed_regime_probabilities.regime1.data;
    p_reg2_mat(:,i) = f.smoothed_regime_probabilities.regime2.data;
    p_reg3_mat(:,i) = f.smoothed_regime_probabilities.regime3.data;
    p_reg4_mat(:,i) = f.smoothed_regime_probabilities.regime4.data;
    y_erg(:,i) = y1.*p_reg1_mat(:,i) + y2.*p_reg2_mat(:,i) + y3.*p_reg3_mat(:,i) + y4.*p_reg4_mat(:,i);

    %fprintf('\n *** Processing parameter draw = %i of %i ****', i, nParamDraws);


    % Compute predictive density for each period (t, nSims,nParamDraws)
    [dY_erg_1_mat(:,:,i),dY_erg_2_mat(:,:,i),dY_erg_3_mat(:,:,i),dY_erg_4_mat(:,:,i), dY_erg_mat(:,:,i),dY_erg_weighted_mat(:,:,i)] = fPredictiveDist_2reg_endoprob_2sepMarkov(MF,FF,TR,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws2,normal,1,1);

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);



% Percentiles of Predictive Density
p_reg1_LB = prctile(p_reg1_mat,5,2);
p_reg1_M = prctile(p_reg1_mat,50,2);
p_reg1_UB = prctile(p_reg1_mat,95,2);

p_reg2_LB = prctile(p_reg2_mat,5,2);
p_reg2_M = prctile(p_reg2_mat,50,2);
p_reg2_UB = prctile(p_reg2_mat,95,2);

q10 = prctile(y_erg,10,2);
q50 = prctile(y_erg,50,2);
q90 = prctile(y_erg,90,2);

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_bar = (reshape(dY_erg_mat,length(MF),nDraws2*nParamDraws));
y_erg_1_bar = (reshape(dY_erg_1_mat,length(MF),nDraws2*nParamDraws));
y_erg_2_bar = (reshape(dY_erg_2_mat,length(MF),nDraws2*nParamDraws));


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

pmode

params_CI = [params_pmode-(params_M_Struct-params_LB_Struct) params_pmode params_pmode+(params_UB_Struct-params_M_Struct)];

disp('Transition Probabilities')
round(params_CI(20:25,:),2)
disp('Coefficients')
round([params_CI(1:19,:);params_CI(26:end-2,:)],2)
disp('Volatilities')
round((params_CI(end-1:end,:)),2) % no need to take sqrt in structural form

diary off

%%

%==========================================================================
% Prepare for density plots
%==========================================================================

colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];


md = '2008-Oct';
mydate = find(datenum(md,dataformat)==dates);
od = '2018-Dec';
mydate2 = find(datenum(od,dataformat)==dates);

mcolor = colors(2,:);
ocolor = colors(1,:);

% Kernel smoothing function October 2008

md = '2007-Aug';
mydate = find(datenum(md,dataformat)==dates);

fig=figure;
% Ergodic
[pdf,xi]=ksdensity(y_erg_bar(mydate,:));
% Simulated at posterior mode
[pdf11,xi11]=ksdensity(dY_sim_1(mydate,:));
[pdf22,xi22]=ksdensity(dY_sim_2(mydate,:));
[pdf33,xi33]=ksdensity(dY_sim_3(mydate,:));
[pdf44,xi44]=ksdensity(dY_sim_4(mydate,:));

[pdf5,xi5]=ksdensity(dY_sim(mydate,:));
hold on
plot(xi11,pdf11,'--','Color',[1 0.9 0.1],'LineWidth', 1.5,'DisplayName',['High growth, low vol']);
plot(xi22,pdf22,'--','Color',[0.5 0.9 0.9],'LineWidth', 1.5,'DisplayName',['High growth, high vol']);
plot(xi33,pdf33,'--','Color',[0 0.9 0.1],'LineWidth', 1.5,'DisplayName',['Low growth, low vol']);
plot(xi44,pdf44,'--','Color',[0.5 0.1 0.1],'LineWidth', 1.5,'DisplayName',['Low growth, high vol']);
plot(xi5,pdf5,'k-','LineWidth', 2,'DisplayName',['Full']);
plot([2.2 2.2],[0 0.6],'b-','LineWidth',0.8);
title('August-2007 Predictive Density 1-Year-Ahead GDP Growth','Interpreter','Latex','FontSize',12)
hold off
legend(gca,'Location','NorthWest','interpreter','Latex')
legend boxoff
ylabel('PDF','fontsize',10,'interpreter','Latex')
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
 if saveit==1
     print('-dpdf',fig,[fig_folder 'KernelDensity_October2008'],'-bestfit');

     print('-dpdf',fig,[slides_folder model '_KernelDensity_October2008'],'-bestfit');

     %saveas(fig,sprintf('%s.png',[fig_folder 'KernelDensity']));
 end


%% Density comparison at period1=October-2008 and December-2018
fig=figure;
% At fitted distribution (Francesca)
[pdf,xi]=ksdensity(y_erg(mydate,:));
[pdf2,xi2]=ksdensity(y_erg(mydate2,:));

% At ergodic distribution (Pablo)
[pdf3,xi3]=ksdensity(y_erg_bar(mydate,:));
[pdf4,xi4]=ksdensity(y_erg_bar(mydate2,:));

% At posterior mode (Pablo)
[pdf5,xi5]=ksdensity(dY_sim(mydate,:));
[pdf6,xi6]=ksdensity(dY_sim(mydate2,:));

hold on
plot(xi,pdf,'Color',mcolor,'LineWidth', 2,'DisplayName',[md '-MS Ergodic (Fra)']);
%plot([realized(mydate) realized(mydate)],[0 max([pdf,pdf2])],'--','Color',mcolor,'LineWidth',2,'DisplayName',[md ' - Realized'])
plot(xi2,pdf2,'Color',ocolor,'LineWidth', 2,'DisplayName',[od '-MS Ergodic (Fra)']);

plot(xi3,pdf3,'--','Color',mcolor,'LineWidth', 2,'DisplayName',[md '-MS Ergodic (Pablo)']);
plot(xi4,pdf4,'--','Color',ocolor,'LineWidth', 2,'DisplayName',[od '-MS Ergodic (Pablo)']);


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
 if saveit==1
     print('-dpdf',fig,[fig_folder 'KernelDensity_2008_2018'],'-bestfit');
     %saveas(fig,sprintf('%s.png',[fig_folder 'KernelDensity']));
 end


%-------------------------------------------------------------------
% Density comparison at period1=October-2008 and December-2018
%-------------------------------------------------------------------

fig=figure;
% At fitted distribution (Francesca)
%[pdf,xi]=ksdensity(y_erg(mydate,:));
%[pdf2,xi2]=ksdensity(y_erg(mydate2,:));

% At ergodic distribution (Pablo)
[pdf3,xi3]=ksdensity(y_erg_bar(mydate,:));
[pdf4,xi4]=ksdensity(y_erg_bar(mydate2,:));

% At posterior mode (Pablo)
[pdf5,xi5]=ksdensity(dY_sim(mydate,:));
[pdf6,xi6]=ksdensity(dY_sim(mydate2,:));

hold on
plot(xi5,pdf5,'Color',mcolor,'LineWidth', 2,'DisplayName',[md '-MS Posterior Mode']);
%plot([realized(mydate) realized(mydate)],[0 max([pdf,pdf2])],'--','Color',mcolor,'LineWidth',2,'DisplayName',[md ' - Realized'])
plot(xi6,pdf6,'Color',ocolor,'LineWidth', 2,'DisplayName',[od '-MS Posterior Mode']);

plot(xi3,pdf3,'--','Color',mcolor,'LineWidth', 2,'DisplayName',[md '-MS Ergodic']);
plot(xi4,pdf4,'--','Color',ocolor,'LineWidth', 2,'DisplayName',[od '-MS Ergodic']);


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
xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
 if saveit==1
     print('-dpdf',fig,[fig_folder 'KernelDensity_ergodic_vs_pmode'],'-bestfit');
     %saveas(fig,sprintf('%s.png',[fig_folder 'KernelDensity']));
 end

%%
%------------------------------------------------------------------- 
% Expected Shortfall and Longrise
%-------------------------------------------------------------------

for tt=1:size(y_erg_bar,1)
    pdY = fitdist(y_erg_bar(tt,:)','Kernel','Kernel','epanechnikov');
    shorfall_cutoff = icdf(pdY,0.05);
    longrise_cutoff = icdf(pdY,0.95);
    shortfall_vec   = y_erg_bar(tt,y_erg_bar(tt,:)<shorfall_cutoff);
    longrise_vec    = y_erg_bar(tt,y_erg_bar(tt,:)>longrise_cutoff);

    EST(tt,1) = mean(shortfall_vec);
    ELR(tt,1) = mean(longrise_vec);
    if mod(tt,100)==0
    fprintf('\n Computing EST and ERT for obs %i ',tt);
    end
end

%% ========================================================================
% Figure: ESF AND ELR
% =========================================================================
% Set colors
colors = cbrewer('div', 'RdYlBu', 64);

fig=figure; clf;
hold on
% Plot percentiles from MS
l1=plot(dates(sd:ed), EST(sd:ed),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Expected Shortfall');
l2=plot(dates(sd:ed), ELR(sd:ed),'--','Color',colors(15,:),'LineWidth', 2,'DisplayName','Expected Longrise');
rr=recessionplot;
hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
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


%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

md = '2007-Sep';
mydate = find(datenum(md,dataformat)==dates);
od = '2007-Sep';
mydate2 = find(datenum(od,dataformat)==dates);

mcolor = colors(1,:);
ocolor = colors(2,:);

if exist('ResMatch')==0
    load('results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};

for tt=1:length(period_pdf)
    tind(tt) = find(datenum(char(period_pdf(tt)),dataformat)==dates);
end

% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST(tind,:);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -8;
xmax = 8;
for t=1:length(period_pdf)
    % plot figure for in-sample pdfs
    if t==1
        ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    end
    % PDF Plot
    %yyaxis left
    hold on;
    plot(ax1,ResMatch.YY,PST(tind(t),:)','-','Color',dcolor(1,:),'LineWidth',2,'DisplayName',[char(period_pdf(t)) ' - QR']);
    %hold(ax1,'on')
    %yyaxis right
    [pdf,xi]=ksdensity(y_erg_bar(tind(t),:));
    plot(ax1,xi,pdf,'-.','LineWidth', 2,'Color',dcolor(2,:),'DisplayName',[char(period_pdf(t)) ' - MS-SVAR']);
end
axis tight
yl = ylim;
%for t=1:length(period_pdf)
%    plot(ax1,[realized(tind(t)) realized(tind(t))],[0 yl(2)],'--','Color',dcolor(t,:),'LineWidth',2,'DisplayName',[char(period_pdf(t)) ' - Realized'])
%end
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


%-------------------------------------------------------------------
% COVID DISTRIBUTIONS
%-------------------------------------------------------------------
% Collect parameters
params=[results.pop.x];

myparams=a2tilde_to_a(params);

p_reg1_matcovid = zeros(length(db.GDPGH.data)-nlags,options.N);
p_reg2_matcovid = zeros(length(db.GDPGH.data)-nlags,options.N);
y_erg_covid = zeros(length(db.GDPGH.data)-nlags,options.N);


% Index for March 2019
t_march_2019 = find(datenum('2019-Mar',dataformat)==dates);



 
% Covid_y = zeros(length(db_full.GDPGH.data),options.N);
% 
inputformat = 'yyyy-MMM';
covid_dates = datenum((datetime('1973-Jan','InputFormat',inputformat)):calmonths(1):(datetime('2020-May','InputFormat',inputformat)))';


%yd = '2020-Apr';
%mydate3 = find(datenum(yd,dataformat)==covid_dates);

%-------------------------------------------------------------------
% Distributions for MAR-19 VINTAGE
%-------------------------------------------------------------------

% Observables from Apr-2019 to March-2020
FF_COVID = db_full.FF.data(find(datenum('2019-Apr',dataformat)==covid_dates):find(datenum('2020-Mar',dataformat)==covid_dates));
MF_COVID = db_full.MF.data(find(datenum('2019-Apr',dataformat)==covid_dates):find(datenum('2020-Mar',dataformat)==covid_dates));
TR_COVID = db_full.TRENDH.data(find(datenum('2019-Apr',dataformat)==covid_dates):find(datenum('2020-Mar',dataformat)==covid_dates));

% Read FF AND MF FROM VINTAGE FILES
DATA_MAR13 = readtable('Data/MacroRisk_March13_US_only.xlsx','Sheet','DFM_73_Monthly');
FF_COVID = DATA_MAR13.FF(556:567);
MF_COVID = DATA_MAR13.MF_US(556:567);
TR_COVID = DATA_MAR13.TREND_US(556:567);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_3_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_4_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for Mar13...");

tic
for i=1:options.N
[dY_erg_1_matcovid(:,:,i),dY_erg_2_matcovid(:,:,i),dY_erg_3_matcovid(:,:,i),dY_erg_4_matcovid(:,:,i), dY_erg_matcovid(:,:,i)] = ...
    fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF_COVID,FF_COVID,TR_COVID,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,prob_reg34_new(t_march_2019),prob_reg24_new(t_march_2019),i);
updateWaitbar(); %#ok<PFBNS>
end


time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);


% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_mar13 = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));


%-------------------------------------------------------------------
% Distributions for APR2-2020 VINTAGE
%-------------------------------------------------------------------

% Read FF AND MF FROM VINTAGE FILES
DATA_APR2 = readtable('Data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
FF_COVID = DATA_APR2.FF(556:567);
MF_COVID = DATA_APR2.MF_US(556:567);
TR_COVID = DATA_APR2.TREND_US(556:567);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_3_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_4_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);


% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for Apr2...");

tic
for i=1:options.N
[dY_erg_1_matcovid(:,:,i),dY_erg_2_matcovid(:,:,i),dY_erg_3_matcovid(:,:,i),dY_erg_4_matcovid(:,:,i), dY_erg_matcovid(:,:,i)] = ...
    fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF_COVID,FF_COVID,TR_COVID,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,prob_reg34_new(t_march_2019),prob_reg24_new(t_march_2019),i);
updateWaitbar(); %#ok<PFBNS>
end

time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_apr2 = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));

%%
%-----------------------------------------------------------------
% Figure: COVID-19 Builup of risk
%-----------------------------------------------------------------
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
% xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;



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


%%
%-------------------------------------------------------------------
% Distributions for NOV8-2020 VINTAGE
%-------------------------------------------------------------------

% Index for September 2019
t_sep_2019 = find(datenum('2019-Sep',dataformat)==dates);


% Read FF AND MF FROM VINTAGE FILES
% DATA_APR2 = readtable('Data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
% FF_COVID = DATA_APR2.FF(556:567);
% MF_COVID = DATA_APR2.MF_US(556:567);
% TR_COVID = DATA_APR2.TREND_US(556:567);
FF_COVID = db_full.FF.data(end-11-1:end-1);
MF_COVID = db_full.MF.data(end-11-1:end-1);
TR_COVID = db_full.TRENDH.data(end-11-1:end-1);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for September (Nov 8)...");

tic
for i=1:options.N
[dY_erg_1_matcovid(:,:,i),dY_erg_2_matcovid(:,:,i),dY_erg_3_matcovid(:,:,i),dY_erg_4_matcovid(:,:,i), dY_erg_matcovid(:,:,i)] = ...
    fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF_COVID,FF_COVID,TR_COVID,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,prob_reg34_new(t_sep_2019),prob_reg24_new(t_sep_2019),i);
updateWaitbar(); %#ok<PFBNS>
end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_nov8_sep = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));

%% Index for October 2019
t_oct_2019 = find(datenum('2019-Oct',dataformat)==dates);


% Read FF AND MF FROM VINTAGE FILES
% DATA_APR2 = readtable('Data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
% FF_COVID = DATA_APR2.FF(556:567);
% MF_COVID = DATA_APR2.MF_US(556:567);
% TR_COVID = DATA_APR2.TREND_US(556:567);
FF_COVID = db_full.FF.data(end-11:end);
MF_COVID = db_full.MF.data(end-11:end);
TR_COVID = db_full.TRENDH.data(end-11:end);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for October (Nov 8)...");

tic
for i=1:options.N
[dY_erg_1_matcovid(:,:,i),dY_erg_2_matcovid(:,:,i),dY_erg_3_matcovid(:,:,i),dY_erg_4_matcovid(:,:,i), dY_erg_matcovid(:,:,i)] = ...
    fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF_COVID,FF_COVID,TR_COVID,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,prob_reg34_new(t_sep_2019),prob_reg24_new(t_sep_2019),i);
updateWaitbar(); %#ok<PFBNS>
end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_nov8_oct = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));

%% Index for March 2019

% Read FF AND MF FROM VINTAGE FILES
% DATA_APR2 = readtable('Data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
% FF_COVID = DATA_APR2.FF(556:567);
% MF_COVID = DATA_APR2.MF_US(556:567);
% TR_COVID = DATA_APR2.TREND_US(556:567);
FF_COVID = db_full.FF.data(end-11-7:end-7);
MF_COVID = db_full.MF.data(end-11-7:end-7);
TR_COVID = db_full.TRENDH.data(end-11-7:end-7);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for March (Nov 8)...");

tic
for i=1:options.N
[dY_erg_1_matcovid(:,:,i),dY_erg_2_matcovid(:,:,i),dY_erg_3_matcovid(:,:,i),dY_erg_4_matcovid(:,:,i), dY_erg_matcovid(:,:,i)] = ...
    fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF_COVID,FF_COVID,TR_COVID,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,prob_reg34_new(t_march_2019),prob_reg24_new(t_march_2019),i);
updateWaitbar(); %#ok<PFBNS>
end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_nov8_mar = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));

%%
%-----------------------------------------------------------------
% Figure: COVID-19 October/2020
%-----------------------------------------------------------------
fig=figure;
[pdf_sep,xi_sep]=ksdensity(y_erg_barcovid_nov8_sep(1,:));
[pdf_oct,xi_oct]=ksdensity(y_erg_barcovid_nov8_oct(1,:));
[pdf_mar,xi_mar]=ksdensity(y_erg_barcovid_nov8_mar(1,:));

hold on
l1=plot(xi_oct,pdf_oct,'k-','LineWidth', 4,'DisplayName',['October-2020 (Nov-8 vintage)']);
l3=plot(xi_oct,pdf_sep,'LineWidth', 3,'Color',dcolor(1,:),'DisplayName',['September-2020 (Nov-8 vintage)']);
l2=plot(xi3,pdf3,'LineWidth', 3,'Color',dcolor(2,:),'DisplayName','March-2020 (Apr-2 vintage)');
l4=plot(xi_mar,pdf_mar,'k--','LineWidth', 3,'Color',dcolor(2,:),'DisplayName','March-2020 (Nov-8 vintage)');
hleg = legend([l1 l3 l2 l4],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('Predictive Density: March, September and October-2020','Interpreter','Latex','FontSize',9)
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

if saveit==1
    print('-dpdf',fig,[fig_folder 'Mar_Oct_2020_alternative'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end

