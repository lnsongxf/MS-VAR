%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all; clc;

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
saveit    = 0; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

% Data vintage, sample and country selection
datafilename = 'MacroRisk_November30.xlsx';
sheetuse     = 'DFM_73_Monthly';
start_date   = '1973-Jan';
end_date     = '2020-Oct';
ctry         = 'US';


% VAR configuration 
nlags=1;
exog={};
constant=true;
panel=[];

% Create date formats for plotting
inputformat = 'yyyy-MMM';
dataformat  = 'yyyy-mmm';
start_plot    = '1973-Feb';
end_plot      = '2019-Oct';
end_plot_full = '2020-Oct';

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Vector of dates for plotting
dates   = dates_full(sd:ed);

% Figure options
FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);
left_color = [0 0 0];
right_color = colors(2,:);

%% Load RISE
rise_startup()

%% Prepare Data

%load data
[db, db_full,startdb,enddb,tex] = data(datafilename,sheetuse,start_date,end_date,ctry);

% Collect MF, FF and trend (estimation sample)
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);
TR = db.TRENDH.data(nlags+1:end);

% Collect MF, FF and trend (full sample)
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);
TR_full = db_full.TRENDH.data(nlags+1:end);


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
model = 'SVAR_StateDependent_2Regimes';

% Model label
model_label = [datafilename(1:end-5) '_' ctry '_' startdb '_' enddb '_const' num2str(const) '_normal' num2str(normal)];

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

slides_folder = ['Results/Figures/' model '/'];
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
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'
            }};
    end
    prob_params = {{'a12','a21','b12','c12','b21','c21'}};
elseif const==0
    if normal==1
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)-c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(-b21*(FF)-c21*(MF)))'
            }};
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(+b21*(FF)-c21*(MF)))'
            }};
    end
    prob_params = {{'b12','c12','b21','c21'}};
end

markov_chains=struct('name','sync',...
    'number_of_states',2,...
    'controlled_parameters',{switches},...   % these correspond to the parameters that are switching
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
        %         switch_prior.a12={0.5,0.5,0.5,'normal'};
        %         switch_prior.a21={0.5,0.5,0.5,'normal'};
        %         switch_prior.b12={-1,-1,0.5,'normal'};
        %         switch_prior.c12={1,1,0.5,'normal'};
        %         switch_prior.b21={1,1,0.5,'normal'};
        %         switch_prior.c21={-1,-1,0.5,'normal'};
        switch_prior.a12={0,0,2,'normal'};
        switch_prior.b12={0,0,2,'normal'};
        switch_prior.c12={0,0,2,'normal'};
        switch_prior.a21={0,0,2,'normal'};
        switch_prior.b21={0,0,2,'normal'};
        switch_prior.c21={0,0,2,'normal'};
    elseif normal==0
        switch_prior.a12={0.5,0.5,0.5,'normal'};
        switch_prior.a21={0.5,0.5,0.5,'normal'};
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
    end
elseif const==0
    if normal==1
        %         switch_prior.b12={0,0,2,'normal'};
        %         switch_prior.c12={0,0,2,'normal'};
        %         switch_prior.b21={0,0,2,'normal'};
        %         switch_prior.c21={0,0,2,'normal'};
        switch_prior.b12={-0.25,-0.25,0.1,'normal'};
        switch_prior.c12={0.25,0.25,0.1,'normal'};
        switch_prior.b21={0.25,0.25,0.1,'normal'};
        switch_prior.c21={-0.25,-0.25,0.1,'normal'};
    elseif normal==0
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
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

% Smoothed and Filtered probabilities
[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;

p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;


%----------------------------------------------------
% Figure: Fitted value and realization of dependent variable
%----------------------------------------------------

% inputformat = 'yyyy-MMM';
% dates = datenum((datetime('1973-Jan','InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime('2019-Oct','InputFormat',inputformat)))';
% dataformat  = 'yyyy-mmm';
% 
% FontSize = 16;
% numticks = 48;
% figSize = [12 6];
% linestyle = {'-','--',':'};
% 
% start_plot = '1973-Feb';
% end_plot   = '2019-Oct';
% sd = find(datenum(start_plot,dataformat)==dates);
% ed = find(datenum(end_plot,dataformat)==dates);

y_fit = fit1.*p_reg1 + fit2.*p_reg2;

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
% colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
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
title('Smoothed Regime Probabilies and Realized GDP','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end

%----------------------------------------------------
% Figure: Regime Probabilities Against Data and Factors
%----------------------------------------------------
figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
% colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
close

numticks = 48;
left_color = [0 0 0];
right_color = colors(2,:);

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
plot(dates(sd:ed),db.MF.data(sd:ed),'c:','LineWidth', 1.5)
% plot(dates(sd:ed),realized(sd:ed),'k','LineWidth', 1.5)
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
legend('$f_t$','$m_t$','interpreter','Latex','Orientation','Horizontal')
legend boxoff
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
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs_X'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs_X']));
end



%% Time Series of Quantiles - QR Comparison
%----------------------------------------------------
% Figure: QR quantiles vs MS regime fitted values
%----------------------------------------------------

figure;
set(gca,'ColorOrder','factory');
colr = get(gca,'ColorOrder');
% colors = [colr(2,:);colr(1,:);colr(5,:)];
close

% load(['quantiles_QR_FN_1973']);
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
prob_reg2_new = zeros(1,length(p_reg2));
prob_reg2_new(p_reg2(sd:ed)>0.9) = 1;

% Draw shocks for simulation
rng(123);
shocks_sim = randn(length(prob_reg2_new),nDraws);

% Collect MF, FF and trend
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);
TR = db.TRENDH.data(nlags+1:end);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws

    for tt= 1:length(prob_reg2_new)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

        if tt<12
            st_lag = 1; % initialize in bad regime (Depends on when the sample starts)

            % Set st
            st_use = 1;

        else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag = prob_reg2_new(tt-11);

            % Get transition probabilities from t-11:tt
            if const==1
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF(tt-11:tt))-pmode.c12*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF(tt-11:tt))-pmode.c21*(MF(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
                end
            elseif const==0
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(-pmode.b12*(FF(tt-11:tt))-pmode.c12*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(-pmode.b21*(FF(tt-11:tt))-pmode.c21*(MF(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
                end
            end
            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;


            % Draw regimes from t-12:t 
            for tt2 = 1:12

                % Draw the Markov Chain for period t-12:t
                udraw = rand(1);

                if st_lag == 0 % started in good regime
                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end

                else % start in bad regime

                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end

                end

                st_lag = st(tt2,1);
            end

            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use = st(tt2);
        end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;


        % Good regime in period t
        dY_sim_1(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF(tt)...
            - pmode.a0_3_2_sync_1*MF(tt) + pmode.s_3_3_sync_1*shocks_sim(tt,nsim))...
            + TR(tt);

        % Bad regime in period t
        dY_sim_2(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF(tt)...
            - pmode.a0_3_2_sync_2*MF(tt)+ pmode.s_3_3_sync_2*shocks_sim(tt,nsim))...
            + TR(tt);

        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_good + dY_sim_2(tt,nsim)*(1-IND_good);

        St_sim(tt,nsim) = IND_good;

        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt);

    end

    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end

% Ergodic probability of regime 1;
p_reg1_sim(:,1) =  sum(St_sim,2)/size(sum(St_sim),2);
p_reg2_sim(:,1) =  sum((1-St_sim),2)/size(sum(St_sim),2);


% Compute percentiles
dY_sim_25 = prctile(dY_sim',25)';
dY_sim_75 = prctile(dY_sim',75)';
dY_sim_10 = prctile(dY_sim',10)';
dY_sim_90 = prctile(dY_sim',90)';

% Compute Regime specific means
dY_sim1_mean = mean(dY_sim_1,2);
dY_sim2_mean = mean(dY_sim_2,2);


% COLLECT BEGIN/END DATE OF TRANSITIONS INTO THE BAD REGIME FOR SHADING
count1 = 1;
count2 = 1;

regime2_shade=[NaN, NaN];

if prob_reg2_new(1)==1
    % beggining of regime switch
    regime2_shade(1,1) = dates(1+12);
    count1 = 2;
end

for tt=1:length(prob_reg2_new)-12


    if tt>2

        if prob_reg2_new(tt-1)==0 && prob_reg2_new(tt)==1
            % beggining of regime switch
            regime2_shade(count1,1) = dates(tt+12);
            count1 = count1+1;

        elseif prob_reg2_new(tt-1)==1 && prob_reg2_new(tt)==0
            % end of regime switch
            regime2_shade(count2,2) = dates(tt+12);
            count2 = count2+1;

        end



    end

end

% % Check last period
% if prob_reg2_new(end)==1
%     % beggining of regime switch
%     regime2_shade(end,2) = dates(end);
% elseif prob_reg2_new(end)==0
%     % beggining of regime switch
%     regime2_shade(end,1) = dates(end);    
% end

%% Recursive updating of s(t) from simulations of Markov-chain (FULL SET)



% Number of draws for simulting predictive density
nDraws = 20000;

% Set a 0.9 threshold to switch from good to bad regime
% prob_reg2_new_full = zeros(1,length(p_reg2));
prob_reg2_new_full = zeros(1,length(dates_full));
prob_reg2_new_full(p_reg2(sd:ed)>0.9) = 1;
% prob_reg2_new_full(p_reg2(sd:ed)>0.7) = 1;

% Draw shocks for simulation
rng(123);
shocks_sim_full = randn(length(prob_reg2_new_full),nDraws);

% Collect MF, FF and trend
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);
TR_full = db_full.TRENDH.data(nlags+1:end);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws

    for tt= 1:length(prob_reg2_new_full)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

        if tt<12
            st_lag = 1; % initialize in bad regime (Depends on when the sample starts)

            % Set st
            st_use = 1;

        else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag = prob_reg2_new_full(tt-11);

            % Get transition probabilities from t-11:tt
            if const==1
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF_full(tt-11:tt))-pmode.c12*(MF_full(tt-11:tt))));
                    p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF_full(tt-11:tt))-pmode.c21*(MF_full(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF_full(tt-11:tt))+pmode.c12.*(MF_full(tt-11:tt))));
                    p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF_full(tt-11:tt))-pmode.c21.*(MF_full(tt-11:tt))));
                end
            elseif const==0
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(-pmode.b12*(FF_full(tt-11:tt))-pmode.c12*(MF_full(tt-11:tt))));
                    p21 = 1./(1+exp(-pmode.b21*(FF_full(tt-11:tt))-pmode.c21*(MF_full(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(-pmode.b12.*(FF_full(tt-11:tt))+pmode.c12.*(MF_full(tt-11:tt))));
                    p21 = 1./(1+exp(+pmode.b21.*(FF_full(tt-11:tt))-pmode.c21.*(MF_full(tt-11:tt))));
                end
            end


            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;


            % Draw regimes from t-12:t 
            for tt2 = 1:12

                % Draw the Markov Chain for period t-12:t
                udraw = rand(1);

                if st_lag == 0 % started in good regime
                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end

                else % start in bad regime

                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end

                end

                st_lag = st(tt2,1);
            end

            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use = st(tt2);
        end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good_full = 1 - st_use;


        % Good regime in period t
        dY_sim_1_full(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF_full(tt)...
            - pmode.a0_3_2_sync_1*MF_full(tt) + pmode.s_3_3_sync_1*shocks_sim_full(tt,nsim))...
            + TR_full(tt);

        % Bad regime in period t
        dY_sim_2_full(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF_full(tt)...
            - pmode.a0_3_2_sync_2*MF_full(tt)+ pmode.s_3_3_sync_2*shocks_sim_full(tt,nsim))...
            + TR_full(tt);

        dY_sim_full(tt,nsim) = dY_sim_1_full(tt,nsim)*IND_good_full + dY_sim_2_full(tt,nsim)*(1-IND_good_full);

        St_sim_full(tt,nsim) = IND_good_full;

        % Weight regime simulation using filtered probabilities
%         dY_sim_weighted_full(tt,nsim) = dY_sim_1_full(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2_full(tt,nsim)*p_reg2_filtered(tt);

    end

    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end

% Ergodic probability of regime 1;
p_reg1_sim_full(:,1) =  sum(St_sim_full,2)/size(sum(St_sim_full),2);
p_reg2_sim_full(:,1) =  sum((1-St_sim_full),2)/size(sum(St_sim_full),2);


% Compute percentiles
dY_sim_25_full = prctile(dY_sim_full',25)';
dY_sim_75_full = prctile(dY_sim_full',75)';
dY_sim_10_full = prctile(dY_sim_full',10)';
dY_sim_90_full = prctile(dY_sim_full',90)';

% Compute Regime specific means
dY_sim1_mean_full = mean(dY_sim_1_full,2);
dY_sim2_mean_full = mean(dY_sim_2_full,2);

% COLLECT BEGIN/END DATE OF TRANSITIONS INTO THE BAD REGIME FOR SHADING
count1 = 1;
count2 = 1;

regime2_shade_full=[NaN, NaN];

if prob_reg2_new_full(1)==1
    % beggining of regime switch
    regime2_shade_full(1,1) = dates_full(1+12);
    count1 = 2;
end

for tt=1:length(prob_reg2_new_full)-12


    if tt>2

        if prob_reg2_new_full(tt-1)==0 && prob_reg2_new_full(tt)==1
            % beggining of regime switch
            regime2_shade_full(count1,1) = dates_full(tt+12);
            count1 = count1+1;

        elseif prob_reg2_new_full(tt-1)==1 && prob_reg2_new_full(tt)==0
            % end of regime switch
            regime2_shade_full(count2,2) = dates_full(tt+12);
            count2 = count2+1;

        end



    end

end

% % Check last period
% if prob_reg2_new_full(end)==1
%     % beggining of regime switch
%     regime2_shade_full(end,2) = dates_full(end);
% % elseif prob_reg2_new_full(end)==0
% %     % beggining of regime switch
% %     regime2_shade_full(end,1) = dates_full(end);    
% end

%% Figure: Compare Real-Time Estimate of Prob of Bad Regime 

fig=figure; clf;
l1=plot(dates(sd:ed),p_reg2,'r-','DisplayName','p(bad regime) filtered'); hold on;
l2=plot(dates(sd:ed),p_reg2_sim,'b--','DisplayName','p(bad regime) simulated');
rr = recessionplot;
hleg = legend([l1 l2],'Orientation','Vertical','Location','NE','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight; ylim([0 1]);
setfig(fig,dates,sd,ed,numticks,FontSize,figSize,'Probability of Bad Regime: Filtered vs Real-Time');

%% Figure: Compare Real-Time Estimate of Prob of Bad Regime (FULL SET)

fig=figure; clf;
l2=plot(dates_full(sd:ed_full),p_reg2_sim_full,'-','LineWidth', 1.5,'DisplayName','p(bad regime) simulated'); hold on;
% hleg = legend(l2,'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
rr = recessionplot;
axis tight
setfig(fig,dates_full,sd,ed_full,numticks,FontSize,figSize,'Probability of Bad Regime: Full Sample');
if saveit==1
    print('-dpdf',fig,[fig_folder 'Prob_bad_regime'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
end

%% Figure: Transition probabilities Normal-to-Bad and Bad-to-Normal


if const==1
p12_fitted = 1./(1+exp(pmode.a12-pmode.b12.*(FF_full)+pmode.c12.*(MF_full)));
p12_fitted_ffonly = 1./(1+exp(pmode.a12-pmode.b12.*(FF_full)+pmode.c12.*0*(MF_full)));

p21_fitted = 1./(1+exp(pmode.a21+pmode.b21.*(FF_full)-pmode.c21.*(MF_full)));
p21_fitted_ffonly = 1./(1+exp(pmode.a21+pmode.b21.*(FF_full)-pmode.c21.*0*(MF_full)));
else
p12_fitted = 1./(1+exp(1-pmode.b12.*(FF_full)+pmode.c12.*(MF_full)));
p12_fitted_ffonly = 1./(1+exp(-pmode.b12.*(FF_full)+pmode.c12.*0*(MF_full)));

p21_fitted = 1./(1+exp(1+pmode.b21.*(FF_full)-pmode.c21.*(MF_full)));
p21_fitted_ffonly = 1./(1+exp(+pmode.b21.*(FF_full)-pmode.c21.*0*(MF_full)));
end

fig=figure
set(fig,'defaultAxesColorOrder',[right_color]);
l1=plot(dates_full(sd:ed_full),p12_fitted,'color',colors2(1,:),'DisplayName','$\hat{p}_{12}$: normal-to-bad','LineWidth',3); hold on;
rr = recessionplot;
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
%yticks([.2 .4 .6 .8]); yticklabels({'0.2','0.4','0.6', '0.8'})
setfig(fig,dates_full,sd,ed_full,numticks,FontSize,figSize,[]);

fig=figure
l1=plot(dates_full(sd:ed_full),p21_fitted,'color',colors2(3,:),'DisplayName','$\hat{p}_{21}$: bad-to-normal','LineWidth',3); hold on;
ylim([0, 1]);
rr = recessionplot;
hleg = legend([l1],'Orientation','vertical','Location','North','interpreter','Latex');legend boxoff;
setfig(fig,dates_full,sd,ed_full,numticks,FontSize,figSize,[]);


%% Figure: 10. 25. 75. 90. percentiles

fig=figure; clf;
hold on
% Plot percentiles from MS
l1=plot(dates(sd:ed), dY_sim_10(sd:ed),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(dates(sd:ed), dY_sim_90(sd:ed),'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
rr = recessionplot;
% plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
%l5=plot(dates(sd:ed),YQ_dist_fut(sd:ed,4),'--','Color',colors(20,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(1)*100) 'th Quantile']);
%l6=plot(dates(sd:ed),YQ_dist_fut(sd:ed,14),'--','Color',colors(50,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(3)*100) 'th Quantile']);
% Plot shades of bad regime
hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
setfig(fig,dates,sd,ed,numticks,FontSize,figSize,'Quantilees MS Model');

%% Figure: Quantiles: shorter Sample

start_plot2 = '2007-Jan';
sd2 = find(datenum(start_plot2,dataformat)==dates);

fig=figure; clf;
hold on
l1=plot(dates_full(sd2:ed_full), dY_sim_10_full(sd2:ed_full),'Color',colors(5,:),'LineWidth', 2.5,'DisplayName','10th');
l2=plot(dates_full(sd2:ed_full), dY_sim_25_full(sd2:ed_full),'--','Color',colors(15,:),'LineWidth', 2.5,'DisplayName','25th');
l3=plot(dates_full(sd2:ed_full), dY_sim_75_full(sd2:ed_full),'--','Color',colors(45,:),'LineWidth', 2.5,'DisplayName','75th');
l4=plot(dates_full(sd2:ed_full), dY_sim_90_full(sd2:ed_full),'Color',colors(55,:),'LineWidth', 2.5,'DisplayName','90th');
rr = recessionplot; 
hleg = legend([l1 l2 l3 l4],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
setfig(fig,dates_full,sd2,ed_full,numticks,FontSize,figSize,'Historical Quantiles 2007-2020');
% Save file
print('-dpdf',fig,[slides_folder model '_Quantiles_MS_recent'],'-bestfit');


%% Figure 25. 75. percentiles vs Expected value in each regime

fig=figure; clf;
%yyaxis left
hold on
l1=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l2=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(55,:),'LineWidth', 2,'DisplayName','MS 75th');
l3=plot(dates(sd:ed), dY_sim1_mean(sd:ed),'-.','Color',colors(45,:),'LineWidth', 1.5,'DisplayName','$E(Good Regime)$');
l4=plot(dates(sd:ed), dY_sim2_mean(sd:ed),'--','Color',colors(12,:),'LineWidth', 1.5,'DisplayName','$E(Bad Regime)$');
% Format plot
hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
ylabel('Percent','interpreter','Latex','fontsize',10)
setfig(fig,dates,sd,ed,numticks,FontSize,figSize,'Expected Value in Each Regime');


%% Percentile of the expected value in regime 1 and 2

for tt=1:length(prob_reg2_new)
    %[,Ydens] = ksdensity(dY_sim(tt,:));

    % Fit the period-t predictive distribution
    pdY = fitdist(dY_sim(tt,:)','Kernel','Kernel','epanechnikov');

    percentile_reg1(tt,1) = cdf(pdY,dY_sim1_mean(tt,1));
    percentile_reg2(tt,1) = cdf(pdY,dY_sim2_mean(tt,1));
end

fig=figure;
hold on
l1=plot(dates(sd:ed), percentile_reg1(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','prctl Mean of Good regime');
l2=plot(dates(sd:ed), percentile_reg2(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','prctl Mean of Bad regime');
l3=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg1(sd:ed)),'--','Color',colors(45,:),'LineWidth', 2);
l4=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg2(sd:ed)),'--','Color',colors(15,:),'LineWidth', 2);
% Format plot
hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
setfig(fig,dates,sd,ed,numticks,FontSize,figSize,'Percentiles of Regime Expected Values');

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
options.burnin=10^3;
options.N=2*10^4;
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

%% Posterior Estimates and Confidence Bands

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
% %==========================================================================
% % COMPUTE HISTORICAL QUANTILES FROM POSTERIOR DISTRIBUTION
% %==========================================================================
% 
% % Collect parameters
% params=[results.pop.x];
% 
% myparams=a2tilde_to_a(params);
% 
% p_reg1_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
% p_reg2_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
% y_erg = zeros(length(db.GDPGH.data)-nlags,options.N);
% 
% 
% nDraws2 = 1;
% nParamDraws = options.N;
% 
% % Allocate matrices
% dY_erg_1_mat = NaN(length(MF),nDraws2,nParamDraws);
% dY_erg_2_mat = NaN(length(MF),nDraws2,nParamDraws);
% dY_erg_mat = NaN(length(MF),nDraws2,nParamDraws);
% dY_erg_weighted_mat = NaN(length(MF),nDraws2,nParamDraws);
% 
% % Start waitbar
% updateWaitbar = waitbarParfor(nParamDraws, "Simulating model at posterior draws...");
% 
% % COMPUTE PREDICTIVE DENSITY
% tic
% for i=1:options.N
%     %sol = solve(sv,myparams(:,i));
%     % sv.estim_
% 
%     paramas_use = myparams(:,i);
% 
%     % Compute regime specific densities
%     [Resids,Fits]=residuals(sv,paramas_use);
%     [~,~,~,f]=filter(sv,paramas_use);
% 
%     fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
%     fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
%     resid1 = Resids.shock_3.data(:,1);
%     resid2 = Resids.shock_3.data(:,2);
% 
%     y1 = fit1 + (paramas_use(end-1))*randn;
%     y2 = fit2 + (paramas_use(end))*randn;
% 
%     p_reg1_mat(:,i) = f.smoothed_regime_probabilities.regime1.data;
%     p_reg2_mat(:,i) = f.smoothed_regime_probabilities.regime2.data;
%     y_erg(:,i) = y1.*p_reg1_mat(:,i) + y2.*p_reg2_mat(:,i);
% 
%     %fprintf('\n *** Processing parameter draw = %i of %i ****', i, nParamDraws);
% 
% 
%     % Compute predictive density for each period (t, nSims,nParamDraws)
%     [dY_erg_1_mat(:,:,i),dY_erg_2_mat(:,:,i), dY_erg_mat(:,:,i),dY_erg_weighted_mat(:,:,i)] = fPredictiveDist_2reg_endoprob(MF,FF,TR,myparams,pnames, p_reg1_filtered, p_reg2_filtered,prob_reg2_new, nDraws2,normal);
% 
%     updateWaitbar(); %#ok<PFBNS>
% 
% 
% end
% time_erg = toc;
% 
% fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);
% 
% 
% 
% % Percentiles of Predictive Density
% p_reg1_LB = prctile(p_reg1_mat,5,2);
% p_reg1_M = prctile(p_reg1_mat,50,2);
% p_reg1_UB = prctile(p_reg1_mat,95,2);
% 
% p_reg2_LB = prctile(p_reg2_mat,5,2);
% p_reg2_M = prctile(p_reg2_mat,50,2);
% p_reg2_UB = prctile(p_reg2_mat,95,2);
% 
% q10 = prctile(y_erg,10,2);
% q50 = prctile(y_erg,50,2);
% q90 = prctile(y_erg,90,2);
% 
% % Average out simulations -> returns a TxmParamDraws matrix
% y_erg_bar = (reshape(dY_erg_mat,length(MF),nDraws2*nParamDraws));
% y_erg_1_bar = (reshape(dY_erg_1_mat,length(MF),nDraws2*nParamDraws));
% y_erg_2_bar = (reshape(dY_erg_2_mat,length(MF),nDraws2*nParamDraws));



%%

%==========================================================================
% Density plots
%==========================================================================
return
% Options
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

% Define target periods
outMar2020 = fGetDensityMS('2020-Mar',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);
outApr2020 = fGetDensityMS('2020-Apr',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);
outAug2020 = fGetDensityMS('2020-Aug',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);
outSep2020 = fGetDensityMS('2020-Sep',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);
outOct2020 = fGetDensityMS('2020-Oct',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);
% outDec2020 = fGetDensityMS('2020-Dec',data_use,a2tilde_to_a,results,pnames,nDraws,nParamDraws);

% %% Figure: December 2020 and Good and Bad Regimes
% [pdf_dec,xi_dec]=ksdensity(outDec2020.y_erg_bar(1,:));
% [cdf_dec_base,xic_dec_base]=ksdensity(outDec2020.y_erg_bar(1,:),'function','cdf');
% 
% [pdf_dec_1,xi_dec_1]=ksdensity(outDec2020.y_erg_bar_1(1,:));
% [pdf_dec_2,xi_dec_2]=ksdensity(outDec2020.y_erg_bar_2(1,:));
% 
% fig=figure;
% hold on
% l1=plot(xi_dec,pdf_dec,'k-','LineWidth', 4,'DisplayName',['December-2020 ']);
% l2=plot(xi_dec_1,pdf_dec_1,'r--','LineWidth', 4,'DisplayName',['December-2020 (Good Regime)']);
% l3=plot(xi_dec_2,pdf_dec_2,'g--','LineWidth', 4,'DisplayName',['December-2020  (Bad Regime)']);
% hleg = legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
% title('Predictive Density: December-2020','Interpreter','Latex','FontSize',9)
% ylabel('PDF','fontsize',10,'interpreter','Latex')
% xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% axis tight
% xlim([-8 8])
% set(fig,'PaperOrientation','portrait');
% set(fig, 'PaperSize', figSize);
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'Units','inches');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
% tightfig;
% %%
% if saveit==1
%     print('-dpdf',fig,[fig_folder 'Dec_2020_foreign'],'-bestfit');
%     %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
% end

%% Figure: March, April, August, October 2020 
% [pdf_dec,xi_dec]=ksdensity(outDec2020.y_erg_bar(1,:));
[pdf_mar,xi_mar]=ksdensity(outMar2020.y_erg_bar(1,:));
[pdf_apr,xi_apr]=ksdensity(outApr2020.y_erg_bar(1,:));
%[pdf_aug,xi_aug]=ksdensity(outAug2020.y_erg_bar(1,:));
%[pdf_sep,xi_sep]=ksdensity(outSep2020.y_erg_bar(1,:));
[pdf_oct,xi_oct]=ksdensity(outOct2020.y_erg_bar(1,:));
[cdf_oct,xic_oct]=ksdensity(outOct2020.y_erg_bar(1,:),'function','cdf');
% [pdf_dec,xi_dec]=ksdensity(outDec2020.y_erg_bar(1,:));

[pdf_oct_1,xi_oct_1]=ksdensity(outOct2020.y_erg_bar_1(1,:));
[pdf_oct_2,xi_oct_2]=ksdensity(outOct2020.y_erg_bar_2(1,:));
[pdf_mar_1,xi_mar_1]=ksdensity(outMar2020.y_erg_bar_1(1,:));
[pdf_mar_2,xi_mar_2]=ksdensity(outMar2020.y_erg_bar_2(1,:));


save data_plot_us_bb_1130.mat xi_mar pdf_mar xi_apr pdf_apr xi_oct pdf_oct xi_oct_1 pdf_oct_1 xi_oct_2 pdf_oct_2 xi_mar_1 pdf_mar_1 xi_mar_2 pdf_mar_2
return

%save data_plot_us_iedo_1130.mat xi_mar pdf_mar xi_aug pdf_aug xi_sep pdf_sep xi_oct pdf_oct


colors3 = cbrewer('seq', 'BuPu', 8);
% colors = cbrewer('div', 'RdYlGn', 64);
staff = 2.14; % Q3/Q3
s_wave = -2.45; % Q3/Q3

fig=figure;
hold on
l1=plot(xi_mar,pdf_mar,'Color',colors(1,:),'LineWidth', 3,'DisplayName',['March-2020 ']); hold on;
% l2=plot(xi_apr,pdf_apr,'Color',colors3(8,:),'LineWidth', 3,'DisplayName',['April-2020 ']); hold on;
% l3=plot(xi_aug,pdf_aug,'Color',colors(22,:),'LineWidth', 3,'DisplayName',['August-2020 ']);
%l3=plot(xi_sep,pdf_sep,'Color',colors(22,:),'LineWidth', 3,'DisplayName',['September-2020 ']);
l4=plot(xi_oct,pdf_oct,'Color',colors(60,:),'LineWidth', 3,'DisplayName',['October-2020 ']);
% l5=plot(xi_dec,pdf_dec,'Color',colors(45,:),'LineWidth', 3,'DisplayName',['December-2020 ']);
ylimits = ylim; xlimtis = xlim;
plot(zeros(10,1),linspace(0,ylimits(2),10),'k--');
%l6=plot([staff staff],[0 ylimits(2)],'b','LineWidth',1,'DisplayName',['Staff forecast = ' num2str(round(staff,1)) '\%']);
l7=plot([s_wave s_wave],[0 ylimits(2)],'r','LineWidth',1,'DisplayName',['Second Waves = ' num2str(round(s_wave,1)) '\%']);
% hleg = legend([l1 l2 l3 l4 l5 l6 l7],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
hleg = legend([l1 l4 l7],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('MS-VAR: Densitites for GDP growth over the next 12 months (U.S.)','FontSize',16','Interpreter','Latex');
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'XTick',-16:2:6)
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
if saveit==1
    print('-dpdf',fig,[fig_folder 'Mar_Sep_Oct_2020_US'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end
%%
fig=figure;
hold on
l1=plot(xi_oct,pdf_oct,'k-','LineWidth', 4,'DisplayName',['December-2020 ']);
l2=plot(xi_oct_1,pdf_oct_1,'r--','LineWidth', 4,'DisplayName',['Good Regime']);
l3=plot(xi_oct_2,pdf_oct_2,'g--','LineWidth', 4,'DisplayName',['Bad Regime']);
hleg = legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
ylimits = ylim; % xlimtis = xlim;
plot(zeros(10,1),linspace(0,ylimits(2),10),'k--');
hleg = legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('MS-VAR: Densitites for GDP growth over the next 12 months (U.S.)','FontSize',16','Interpreter','Latex');
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'XTick',-12:2:8)
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
