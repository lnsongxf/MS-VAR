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

switches = {'c(3)','a0(3)','s(3)'};
model = 'SVAR_StateDependent_2Regimes_alternative';

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
colors = [colr(1,:);colr(2,:);colr(3,:);colr(4,:);colr(5,:)];
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
plot(dates(sd:ed), fit2(sd:ed),'Color',colors(1,:),'LineWidth', 1.5)
plot(dates(sd:ed), fit1(sd:ed),'Color',colors(2,:),'LineWidth', 1.5)
j=0;
for iquantile=indq
    j=j+1;
    plot(dates(sd:ed),YQ_dist_fut(sd:ed,iquantile),'-.','Color',colors(j,:),'LineWidth',1.5);
    hold on
end
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {'Bad Regime','Normal Regime',[num2str(Qplot(1)*100) 'th Quantile'],'Median',[num2str(Qplot(3)*100) 'th Quantile']};
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
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF(tt-11:tt))-pmode.c12*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF(tt-11:tt))-pmode.c21*(MF(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
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
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF_full(tt-11:tt))-pmode.c12*(MF_full(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF_full(tt-11:tt))-pmode.c21*(MF_full(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF_full(tt-11:tt))+pmode.c12.*(MF_full(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF_full(tt-11:tt))-pmode.c21.*(MF_full(tt-11:tt))));
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

%%
%----------------------------------------------------
% Figure: Compare Real-Time Estimate of Prob of Bad Regime 
%----------------------------------------------------

fig=figure; clf;
l1=plot(dates(sd:ed),p_reg2,'r-','DisplayName','p(bad regime) filtered'); hold on;
l2=plot(dates(sd:ed),p_reg2_sim,'--','DisplayName','p(bad regime) simulated');
for row=1:length(regime2_shade)
 patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [0 0, 1 1], [0.8 0.8 0.8],'EdgeColor','none'); hold on;
end
hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
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
title('Probability of Bad Regime: Filtered vs Real-Time','FontSize',16','Interpreter','Latex');
tightfig;

%%
%----------------------------------------------------
% Figure: Compare Real-Time Estimate of Prob of Bad Regime (FULL SET)
%----------------------------------------------------


fig=figure; clf;
l2=plot(dates_full(sd:ed_full),p_reg2_sim_full,'-','LineWidth', 1.5,'DisplayName','p(bad regime) simulated'); hold on;
for row=1:length(regime2_shade_full)
 patch([regime2_shade_full(row,1) regime2_shade_full(row,2) regime2_shade_full(row,2) regime2_shade_full(row,1)], [0 0, 1 1], [0.8 0.8 0.8],'EdgeColor','none'); hold on;
end
% hleg = legend(l2,'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
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
% title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'Prob_bad_regime'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_MSFitted_QR']));
end

%% ========================================================================
% PLOT FIGURE 10. 25. 75. 90. percentiles
% =========================================================================

% Set colors
%colors = cbrewer('qual', 'Set2', 8);
colors = cbrewer('div', 'RdYlBu', 64);

fig=figure; clf;
%yyaxis left
hold on
% Plot percentiles from MS
l1=plot(dates(sd:ed), dY_sim_10(sd:ed),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(dates(sd:ed), dY_sim_90(sd:ed),'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');

% plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
l5=plot(dates(sd:ed),YQ_dist_fut(sd:ed,4),'--','Color',colors(20,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(1)*100) 'th Quantile']);
l6=plot(dates(sd:ed),YQ_dist_fut(sd:ed,14),'--','Color',colors(50,:),'LineWidth',1.6,'DisplayName',[num2str(Qplot(3)*100) 'th Quantile']);
% Plot shades of bad regime
for row=1:length(regime2_shade)
patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [-8 -8, 8 8], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
end
hleg = legend([l1 l2 l3 l4 l5 l6],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
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


%% ========================================================================
% PLOT FIGURE 25. 75. percentiles vs Expected value in each regime
% =========================================================================

% Set colors
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);

fig=figure; clf;
%yyaxis left
hold on
l1=plot(dates(sd:ed), dY_sim_25(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l2=plot(dates(sd:ed), dY_sim_75(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l3=plot(dates(sd:ed), dY_sim1_mean(sd:ed),'-.','Color',colors2(2,:),'LineWidth', 1.5,'DisplayName','$E(Good Regime)$');
l4=plot(dates(sd:ed), dY_sim2_mean(sd:ed),'--','Color',colors2(4,:),'LineWidth', 1.5,'DisplayName','$E(Bad Regime)$');

% Plot shades of bad regime
for row=1:length(regime2_shade)
patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [-8 -8, 8 8], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
end

% Format plot
hleg = legend([l1 l2 l3 l4],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
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

%%
%========================================================================
% Percentile of the expected value in regime 1 and 2
%=========================================================================
for tt=1:length(prob_reg2_new)
    %[,Ydens] = ksdensity(dY_sim(tt,:));

    % Fit the period-t predictive distribution
    pdY = fitdist(dY_sim(tt,:)','Kernel','Kernel','epanechnikov');

    percentile_reg1(tt,1) = cdf(pdY,dY_sim1_mean(tt,1));
    percentile_reg2(tt,1) = cdf(pdY,dY_sim2_mean(tt,1));
end



colors = cbrewer('div', 'RdYlBu', 64);

fig=figure;
%yyaxis left
hold on
l1=plot(dates(sd:ed), percentile_reg1(sd:ed),'Color',colors(45,:),'LineWidth', 2,'DisplayName','prctl Mean of Good regime');
l2=plot(dates(sd:ed), percentile_reg2(sd:ed),'Color',colors(15,:),'LineWidth', 2,'DisplayName','prctl Mean of Bad regime');
l3=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg1(sd:ed)),'--','Color',colors(45,:),'LineWidth', 2);
l4=plot(dates(sd:ed), ones(length(dates(sd:ed)),1)*mean(percentile_reg2(sd:ed)),'--','Color',colors(15,:),'LineWidth', 2);

% Plot shades of bad regime
for row=1:length(regime2_shade)
patch([regime2_shade(row,1) regime2_shade(row,2) regime2_shade(row,2) regime2_shade(row,1)], [0 0, 1 1], [0.9 0.9 0.9],'EdgeColor','none'); hold on;
end

% Format plot
hleg = legend([l1 l2],'Orientation','Vertical','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
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

%==========================================================================
% COMPUTE OBJECTS FROM POSTERIOR DISTRIBUTION
%==========================================================================

% Collect parameters
params=[results.pop.x];

myparams=a2tilde_to_a(params);

p_reg1_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
p_reg2_mat = zeros(length(db.GDPGH.data)-nlags,options.N);
y_erg = zeros(length(db.GDPGH.data)-nlags,options.N);


nDraws2 = 1;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_mat = NaN(length(MF),nDraws2,nParamDraws);
dY_erg_2_mat = NaN(length(MF),nDraws2,nParamDraws);
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
    resid1 = Resids.shock_3.data(:,1);
    resid2 = Resids.shock_3.data(:,2);

    y1 = fit1 + (paramas_use(end-1))*randn;
    y2 = fit2 + (paramas_use(end))*randn;

    p_reg1_mat(:,i) = f.smoothed_regime_probabilities.regime1.data;
    p_reg2_mat(:,i) = f.smoothed_regime_probabilities.regime2.data;
    y_erg(:,i) = y1.*p_reg1_mat(:,i) + y2.*p_reg2_mat(:,i);

    %fprintf('\n *** Processing parameter draw = %i of %i ****', i, nParamDraws);


    % Compute predictive density for each period (t, nSims,nParamDraws)
    [dY_erg_1_mat(:,:,i),dY_erg_2_mat(:,:,i), dY_erg_mat(:,:,i),dY_erg_weighted_mat(:,:,i)] = fPredictiveDist_2reg_endoprob(MF,FF,TR,myparams,pnames, p_reg1_filtered, p_reg2_filtered,prob_reg2_new, nDraws2,normal,const);

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
[pdf3,xi3]=ksdensity(dY_sim_1(mydate,:));
[pdf4,xi4]=ksdensity(dY_sim_2(mydate,:));
[pdf5,xi5]=ksdensity(dY_sim(mydate,:));
hold on
plot(xi3,pdf3,'--','Color',[0 0.9 0.1],'LineWidth', 1.5,'DisplayName',['Normal regime']);
plot(xi4,pdf4,'--','Color',[1 0.0 0.0],'LineWidth', 1.5,'DisplayName',['Bad regime']);
plot(xi,pdf,'k-','LineWidth', 2,'DisplayName',['Full']);
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
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);


% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for Mar13...");

tic
for i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
    
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for March 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in March 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_march_2019);
        
        % Get transition probabilities from Apri-2019 to March-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Apr-2019 to March-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);
        
        St_sim_COV_temp(1,nsim) = IND_good;
        
    end

    St_sim_COV(:,i) = mean(St_sim_COV_temp);
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_mar13 =  mean(St_sim_COV);
p_reg2_sim_mar13 =  1-p_reg1_sim_mar13;


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
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for Apr2...");

tic
for i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
    
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for March 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in March 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_march_2019);
        
        % Get transition probabilities from Apri-2019 to March-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Apr-2019 to March-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);
        
        St_sim_COV_temp(1,nsim) = IND_good;
        
    end

    St_sim_COV(:,i) = mean(St_sim_COV_temp);
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_apr2 =  mean(St_sim_COV);
p_reg2_sim_apr2 =  1-p_reg1_sim_apr2;


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
xlim([-12 8])
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

% Index for August 2019
t_aug_2019 = find(datenum('2019-Aug',dataformat)==dates);


% Read FF AND MF FROM VINTAGE FILES
% DATA_APR2 = readtable('Data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
% FF_COVID = DATA_APR2.FF(556:567);
% MF_COVID = DATA_APR2.MF_US(556:567);
% TR_COVID = DATA_APR2.TREND_US(556:567);
FF_COVID = db_full.FF.data(end-11-2:end-2);
MF_COVID = db_full.MF.data(end-11-2:end-2);
TR_COVID = db_full.TRENDH.data(end-11-2:end-2);

% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for August (Nov 8)...");

tic
for i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
    
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for October 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in August 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_aug_2019);
        
        % Get transition probabilities from Sep-2019 to Aug-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Sep-2019 to Aug-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);
        
        St_sim_COV_temp(1,nsim) = IND_good;
        
    end
    
     St_sim_COV(:,i) = mean(St_sim_COV_temp);
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_aug =  mean(St_sim_COV);
p_reg2_sim_aug =  1-p_reg1_sim_aug;

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_nov8_aug = (reshape(dY_erg_matcovid,1,nDraws*nParamDraws));

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
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
    
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for October 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in September 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_sep_2019);
        
        % Get transition probabilities from Oct-2019 to Sep-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Oct-2019 to Sep-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);

        St_sim_COV_temp(1,nsim) = IND_good;
        
    end        
     St_sim_COV(:,i) = mean(St_sim_COV_temp);
     
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_sep =  mean(St_sim_COV);
p_reg2_sim_sep =  1-p_reg1_sim_sep;

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
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
%     ax.XTick = datenum(dates_full(sd:numticks:ed_full));
% datetick('x','yyyy','keepticks')
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for October 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in October 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_oct_2019);
        
        % Get transition probabilities from Nov-2019 to Oct-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Nov-2019 to Oct-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);

        St_sim_COV_temp(1,nsim) = IND_good;
        
    end
    
     St_sim_COV(:,i) = mean(St_sim_COV_temp);
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_oct =  mean(St_sim_COV);
p_reg2_sim_oct =  1-p_reg1_sim_oct;

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
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
    
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for March 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in March 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_march_2019);
        
        % Get transition probabilities from Apri-2019 to March-2020
        if normal
            % Transition probabilities when coefficients have normal prior
            p12 = 1./(1+exp(a12-b12*(FF_COVID)-c12*(MF_COVID)));
            p21 = 1./(1+exp(a21-b21*(FF_COVID)-c21*(MF_COVID)));
        else
            % Transition probabilities when coefficients have gamma prior
            p12 = 1./(1+exp(a12-b12.*(FF_COVID)+c12.*(MF_COVID)));
            p21 = 1./(1+exp(a21+b21.*(FF_COVID)-c21.*(MF_COVID)));
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Apr-2019 to March-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matcovid(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID(tt2)...
            - a0_3_2_sync_1*MF_COVID(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID(tt2)...
            - a0_3_2_sync_2*MF_COVID(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid(1,nsim,i) = dY_erg_1_matcovid(1,nsim,i)*IND_good + dY_erg_2_matcovid(1,nsim,i)*(1-IND_good);
        
        %St_sim_COV(1,nsim) = IND_good;
        
    end
    

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
[pdf_aug,xi_aug]=ksdensity(y_erg_barcovid_nov8_aug(1,:));
[cdf_aug,xic_aug]=ksdensity(y_erg_barcovid_nov8_aug(1,:),'function','cdf');
[pdf_sep,xi_sep]=ksdensity(y_erg_barcovid_nov8_sep(1,:));
[cdf_sep,xic_sep]=ksdensity(y_erg_barcovid_nov8_sep(1,:),'function','cdf');
[pdf_oct,xi_oct]=ksdensity(y_erg_barcovid_nov8_oct(1,:));
[cdf_oct,xic_oct]=ksdensity(y_erg_barcovid_nov8_oct(1,:),'function','cdf');
[pdf_mar,xi_mar]=ksdensity(y_erg_barcovid_nov8_mar(1,:));
[cdf_mar,xic_mar]=ksdensity(y_erg_barcovid_nov8_mar(1,:),'function','cdf');

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

pdf_sep_one_ms = pdf_sep; xi_sep_one_ms = xi_sep;
pdf_oct_one_ms = pdf_oct; xi_oct_one_ms = xi_oct;
save one_ms_chain pdf_sep_one_ms xi_sep_one_ms pdf_oct_one_ms xi_oct_one_ms

%% SEPTEMBER DENSITIES with QR
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

md = '2020-Sep';
mydate = find(datenum(md,dataformat)==dates_full);

mcolor = colors(1,:);
ocolor = colors(2,:);

load('results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=November11_Smooth___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-Oct.mat')

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};
clear tind

tind = find(datenum(char(period_pdf),dataformat)==dates_full)+1;


% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST;

[~,i10] = min(abs(CDF(tind,:)'-0.1));
[~,i25] = min(abs(CDF(tind,:)'-0.25));
[~,imedian] = min(abs(CDF(tind,:)'-0.5));
[~,i10_ms] = min(abs(cdf_sep-0.1));
[~,i25_ms_oct] = min(abs(cdf_sep-0.25));
[~,imedian_ms_oct] = min(abs(cdf_sep-0.5));
ym(t) = max(PST(tind,:)');
ymax = max(ym)*1.05;
ymax = max([(max(ym)*1.05);max(pdf_sep)]);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -8;
xmax = 8;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
plot(ax1,ResMatch.YY,PST(tind,:)','-','Color',dcolor(1,:),'LineWidth',2,'DisplayName',[char(period_pdf) ' - QR']);
plot(ax1,xi_sep,pdf_sep,'LineWidth', 2,'Color',dcolor(2,:),'DisplayName',[char(period_pdf) ' - MS-VAR']);
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
plot(ax1,ResMatch.YY(:,imedian),PST(tind,imedian)','o','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(ResMatch.YY(:,imedian),1)) '\% (QR)']);
plot(ax1,xi_sep(imedian_ms_oct),pdf_sep(imedian_ms_oct),'o','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(xi_sep(imedian_ms_oct),1)) '\% (MS)']);
plot(ax1,ResMatch.YY(:,i25),PST(tind,i25)','d','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25),1)) '\% (QR)']);
plot(ax1,xi_sep(i25_ms_oct),pdf_sep(i25_ms_oct),'d','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_sep(i25_ms_oct),1)) '\% (MS)']);
plot(ax1,ResMatch.YY(:,i10),PST(tind,i10)','s','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(ResMatch.YY(:,i10),1)) '\% (QR)']);
plot(ax1,xi_sep(i10_ms),pdf_sep(i10_ms),'s','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(xi_sep(i10_ms),1)) '\% (MS)']);
axis tight
yl = ylim;
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('Densitites for GDP growth over the next 12 months (September 2020)','FontSize',16','Interpreter','Latex');
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
    print('-dpdf',fig,[fig_folder 'DensityQR_Sep'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end


%% OCTOBER DENSITIES with QR
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

md = '2020-Oct';
mydate = find(datenum(md,dataformat)==dates_full);

mcolor = colors(1,:);
ocolor = colors(2,:);

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};
clear tind

tind = find(datenum(char(period_pdf),dataformat)==dates_full)+1;


% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST;

[~,i10] = min(abs(CDF(tind,:)'-0.1));
[~,i25] = min(abs(CDF(tind,:)'-0.25));
[~,imedian] = min(abs(CDF(tind,:)'-0.5));
[~,i10_ms] = min(abs(cdf_oct-0.1));
[~,i25_ms_oct] = min(abs(cdf_oct-0.25));
[~,imedian_ms_oct] = min(abs(cdf_oct-0.5));
ym(t) = max(PST(tind,:)');
ymax = max(ym)*1.05;
ymax = max([(max(ym)*1.05);max(pdf_oct)]);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -8;
xmax = 8;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
plot(ax1,ResMatch.YY,PST(tind,:)','-','Color',dcolor(1,:),'LineWidth',2,'DisplayName',[char(period_pdf) ' - QR']);
plot(ax1,xi_oct,pdf_oct,'LineWidth', 2,'Color',dcolor(2,:),'DisplayName',[char(period_pdf) ' - MS-VAR']);
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
plot(ax1,ResMatch.YY(:,imedian),PST(tind,imedian)','o','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(ResMatch.YY(:,imedian),1)) '\% (QR)']);
plot(ax1,xi_oct(imedian_ms_oct),pdf_oct(imedian_ms_oct),'o','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(xi_oct(imedian_ms_oct),1)) '\% (MS)']);
plot(ax1,ResMatch.YY(:,i25),PST(tind,i25)','d','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25),1)) '\% (QR)']);
plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'d','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_oct(i25_ms_oct),1)) '\% (MS)']);
plot(ax1,ResMatch.YY(:,i10),PST(tind,i10)','s','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(ResMatch.YY(:,i10),1)) '\% (QR)']);
plot(ax1,xi_oct(i10_ms),pdf_oct(i10_ms),'s','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(xi_oct(i10_ms),1)) '\% (MS)']);
axis tight
yl = ylim;
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('Densitites for GDP growth over the next 12 months (October 2020)','FontSize',16','Interpreter','Latex');
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
    print('-dpdf',fig,[fig_folder 'DensityQR_Oct'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end


%% MARCH, AUGUST AND OCTOBER DENSITIES with MS
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

colors2 = cbrewer('div', 'RdYlGn', 18);

md = '2020-Oct';
md2 = '2020-Aug';
md3 = '2020-Mar';
mydate = find(datenum(md,dataformat)==dates_full);

mcolor = colors(1,:);
ocolor = colors(2,:);

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};
period_pdf2 = {md2};
period_pdf3 = {md3};
clear tind


tind = find(datenum(char(period_pdf),dataformat)==dates_full)+1;
tind2 = find(datenum(char(period_pdf2),dataformat)==dates_full)+1;
tind3 = find(datenum(char(period_pdf3),dataformat)==dates_full)+1;


% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST;

[~,i10_oct] = min(abs(CDF(tind,:)'-0.1));
[~,i25_oct] = min(abs(CDF(tind,:)'-0.25));
[~,imedian_oct] = min(abs(CDF(tind,:)'-0.5));
[~,i10_sep] = min(abs(CDF(tind2,:)'-0.1));
[~,i25_sep] = min(abs(CDF(tind2,:)'-0.25));
[~,imedian_sep] = min(abs(CDF(tind2,:)'-0.5));
[~,i10_mar] = min(abs(CDF(tind3,:)'-0.1));
[~,i25_mar] = min(abs(CDF(tind3,:)'-0.25));
[~,imedian_mar] = min(abs(CDF(tind3,:)'-0.5));
[~,i10_ms] = min(abs(cdf_oct-0.1));
[~,i25_ms_oct] = min(abs(cdf_oct-0.25));
[~,i25_ms_sep] = min(abs(cdf_sep-0.25));
[~,i25_ms_mar] = min(abs(cdf_mar-0.25));
[~,imedian_ms_oct] = min(abs(cdf_oct-0.5));

ym(t) = max(PST(tind2,:)');
ymax = max(ym)*1.05;
ymax = max([(max(ym)*1.05);max(pdf_aug)]);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -12;
xmax = 8;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
% plot(ax1,ResMatch.YY,PST(tind,:)','-','Color',dcolor(1,:),'LineWidth',2,'DisplayName',[char(period_pdf) ' - QR']);
plot(ax1,xi_oct,pdf_oct,'LineWidth', 2,'Color',colors2(16,:),'DisplayName',['2020-Oct']);
% plot(ax1,xi_sep,pdf_sep,'LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Sep - MS-VAR']);
plot(ax1,xi_aug,pdf_aug,'LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Aug']);
plot(ax1,xi_mar,pdf_mar,'LineWidth', 2,'Color',colors2(3,:),'DisplayName',['2020-Mar']);
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
% plot(ax1,ResMatch.YY(:,imedian),PST(tind,imedian)','o','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(ResMatch.YY(:,imedian),1)) '\% (QR)']);
% plot(ax1,xi_oct(imedian_ms_oct),pdf_oct(imedian_ms_oct),'o','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(xi_oct(imedian_ms_oct),1)) '\% (MS)']);
% plot(ax1,ResMatch.YY(:,i25),PST(tind,i25)','d','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25),1)) '\% (QR)']);
plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_oct(i25_ms_oct),1)) '\% (October)']);
% plot(ax1,xi_sep(i25_ms_sep),pdf_sep(i25_ms_sep),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_sep(i25_ms_sep),1)) '\% (September)']);
plot(ax1,xi_aug(i25_ms_aug),pdf_aug(i25_ms_aug),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_aug(i25_ms_aug),1)) '\% (August)']);
plot(ax1,xi_mar(i25_ms_mar),pdf_mar(i25_ms_mar),'d','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_mar(i25_ms_mar),1)) '\% (March)']);
% plot(ax1,ResMatch.YY(:,i10),PST(tind,i10)','s','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(ResMatch.YY(:,i10),1)) '\% (QR)']);
% plot(ax1,xi_oct(i10_ms),pdf_oct(i10_ms),'s','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(xi_oct(i10_ms),1)) '\% (MS)']);
axis tight
yl = ylim;
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('MS-VAR: Densitites for GDP growth over the next 12 months','FontSize',16','Interpreter','Latex');
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
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Sep_Oct'],'-bestfit');
    print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Aug_Oct'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Oct'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end

%% MARCH, AUGUST AND OCTOBER DENSITIES with QR
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

colors2 = cbrewer('div', 'RdYlGn', 18);

md = '2020-Oct';
md2 = '2020-Aug';
md3 = '2020-Mar';
mydate = find(datenum(md,dataformat)==dates_full);

mcolor = colors(1,:);
ocolor = colors(2,:);

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};
period_pdf2 = {md2};
period_pdf3 = {md3};
clear tind


tind = find(datenum(char(period_pdf),dataformat)==dates_full)+1;
tind2 = find(datenum(char(period_pdf2),dataformat)==dates_full)+1;
tind3 = find(datenum(char(period_pdf3),dataformat)==dates_full)+1;


% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST;

[~,i10_oct] = min(abs(CDF(tind,:)'-0.1));
[~,i25_oct] = min(abs(CDF(tind,:)'-0.25));
[~,imedian_oct] = min(abs(CDF(tind,:)'-0.5));
[~,i10_sep] = min(abs(CDF(tind2,:)'-0.1));
[~,i25_sep] = min(abs(CDF(tind2,:)'-0.25));
[~,imedian_sep] = min(abs(CDF(tind2,:)'-0.5));
[~,i10_mar] = min(abs(CDF(tind3,:)'-0.1));
[~,i25_mar] = min(abs(CDF(tind3,:)'-0.25));
[~,imedian_mar] = min(abs(CDF(tind3,:)'-0.5));
[~,i10_ms] = min(abs(cdf_oct-0.1));
[~,i25_ms_oct] = min(abs(cdf_oct-0.25));
[~,i25_ms_sep] = min(abs(cdf_sep-0.25));
[~,i25_ms_mar] = min(abs(cdf_mar-0.25));
[~,imedian_ms_oct] = min(abs(cdf_oct-0.5));

ym(t) = max(PST(tind2,:)');
ymax = max(ym)*1.05;
ymax = max([(max(ym)*1.05);max(pdf_aug)]);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -12;
xmax = 8;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
l1=plot(ax1,ResMatch.YY,PST(tind,:)','-.','Color',colors2(16,:),'LineWidth',2,'DisplayName',['2020-Oct (QR)']);
l11=plot(ax1,xi_oct,pdf_oct,'-','LineWidth', 2,'Color',colors2(16,:),'DisplayName',['2020-Oct (MS-VAR)']);
% l2=plot(ax1,ResMatch.YY,PST(tind2,:)','-.','Color',colors2(6,:),'LineWidth',2,'DisplayName',['2020-Aug (QR)']);
% l22=plot(ax1,xi_aug,pdf_aug,'-','LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Aug (MS-VAR)']);
l3=plot(ax1,ResMatch.YY,PST(tind3,:)','-.','Color',colors2(3,:),'LineWidth',2,'DisplayName',['2020-Mar (QR)']);
l33=plot(ax1,xi_mar,pdf_mar,'-','LineWidth', 2,'Color',colors2(3,:),'DisplayName',['2020-Mar (MS-VAR)']);
l1.Color = [l1.Color 0.5];
% plot(ax1,xi_sep,pdf_sep,'LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Sep - MS-VAR']);
% l2.Color = [l2.Color 0.5];
l3.Color = [l3.Color 0.5];
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
% plot(ax1,ResMatch.YY(:,imedian),PST(tind,imedian)','o','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(ResMatch.YY(:,imedian),1)) '\% (QR)']);
% plot(ax1,xi_oct(imedian_ms_oct),pdf_oct(imedian_ms_oct),'o','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(xi_oct(imedian_ms_oct),1)) '\% (MS)']);
% l6=plot(ax1,ResMatch.YY(:,i25_oct),PST(tind,i25_oct)','d','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25_oct),1)) '\% (October, QR)']);
% plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'s','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_oct(i25_ms_oct),1)) '\% (October, MS-VAR)']);
% l7=plot(ax1,ResMatch.YY(:,i25_sep),PST(tind2,i25_sep)','d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25_sep),1)) '\% (August, QR)']);
% plot(ax1,xi_aug(i25_ms_aug),pdf_aug(i25_ms_aug),'s','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_aug(i25_ms_aug),1)) '\% (August, MS-VAR)']);
% l8=plot(ax1,ResMatch.YY(:,i25_mar),PST(tind3,i25_mar)','d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(ResMatch.YY(:,i25_mar),1)) '\% (March, QR)']);
% plot(ax1,xi_mar(i25_ms_mar),pdf_mar(i25_ms_mar),'s','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_mar(i25_ms_mar),1)) '\% (March, MS-VAR)']);
l6=plot(ax1,ResMatch.YY(:,i25_oct),PST(tind,i25_oct)','d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10);
plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'s','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10);
% l7=plot(ax1,ResMatch.YY(:,i25_sep),PST(tind2,i25_sep)','d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10);
% plot(ax1,xi_aug(i25_ms_aug),pdf_aug(i25_ms_aug),'s','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10);
l8=plot(ax1,ResMatch.YY(:,i25_mar),PST(tind3,i25_mar)','d','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10);
plot(ax1,xi_mar(i25_ms_mar),pdf_mar(i25_ms_mar),'s','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10);
l6.Color = [l6.Color 0.5];
% l7.Color = [l7.Color 0.5];
l8.Color = [l8.Color 0.5];
% plot(ax1,xi_sep(i25_ms_sep),pdf_sep(i25_ms_sep),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_sep(i25_ms_sep),1)) '\% (September)']);
% plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'d','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_oct(i25_ms_oct),1)) '\% (October)']);
% plot(ax1,xi_sep(i25_ms_sep),pdf_sep(i25_ms_sep),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_sep(i25_ms_sep),1)) '\% (September)']);
% plot(ax1,xi_mar(i25_ms_mar),pdf_mar(i25_ms_mar),'d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_mar(i25_ms_mar),1)) '\% (March)']);
% plot(ax1,ResMatch.YY(:,i10),PST(tind,i10)','s','MarkerFaceColor',dcolor(1,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(ResMatch.YY(:,i10),1)) '\% (QR)']);
% plot(ax1,xi_oct(i10_ms),pdf_oct(i10_ms),'s','MarkerFaceColor',dcolor(2,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['10th Quantile = ' num2str(round(xi_oct(i10_ms),1)) '\% (MS)']);
axis tight
yl = ylim;
% hL=legend([l1 l11 l2 l22 l3 l33],'Location','NorthWest');
hL=legend([l1 l11 l3 l33],'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('QR and MS-VAR: Densitites for GDP growth over the next 12 months','FontSize',16','Interpreter','Latex');
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
%     print('-dpdf',fig,[fig_folder 'DensityQR_Mar_Sep_Oct'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityQR_Mar_Aug_Oct'],'-bestfit');
    print('-dpdf',fig,[fig_folder 'DensityQR_Mar_Oct'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end

%%
%==========================================================================
% Historical Quantiles QR and MS
%==========================================================================
colors = cbrewer('div', 'RdYlBu', 64);


fig=figure; clf;
%yyaxis left
hold on
% Plot percentiles from MS
l2=plot(dates_full(sd:ed_full), dY_sim_25_full(sd:ed_full),'--','Color',colors(15,:),'LineWidth', 3,'DisplayName','MS 25th');
l3=plot(dates_full(sd:ed_full), dY_sim_75_full(sd:ed_full),'--','Color',colors(55,:),'LineWidth', 3,'DisplayName','MS 75th');

% plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
l5=plot(dates_full(sd:ed_full),YQ_dist_fut(sd+1:ed_full+1,4),'-','Color',colors(20,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(1)*100) 'th ']);
l6=plot(dates_full(sd:ed_full),YQ_dist_fut(sd+1:ed_full+1,14),'-','Color',colors(50,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(3)*100) 'th ']);
ylim([-10 8]);
rr=recessionplot;
hleg = legend([l2 l3 l5 l6],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))

hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
axis tight
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
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
print('-dpdf',fig,[slides_folder model '_Quantiles_MS_QR'],'-bestfit');


fig=figure; clf;
%yyaxis left
hold on
% Plot percentiles from MS
l2=plot(dates_full(sd2:ed_full), dY_sim_25_full(sd2:ed_full),'--','Color',colors(15,:),'LineWidth', 3,'DisplayName','MS 25th');
l3=plot(dates_full(sd2:ed_full), dY_sim_75_full(sd2:ed_full),'--','Color',colors(55,:),'LineWidth', 3,'DisplayName','MS 75th');

% plot Quantiles YQ_dist_fut(:,4) is the 25th quantile, YQ_dist_fut(:,14)
% is the 75th Quantile. Other quantiles are in quantiles_dist.
l5=plot(dates_full(sd2:ed_full),YQ_dist_fut(sd2+1:ed_full+1,4),'-','Color',colors(20,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(1)*100) 'th ']);
l6=plot(dates_full(sd2:ed_full),YQ_dist_fut(sd2+1:ed_full+1,14),'-','Color',colors(50,:),'LineWidth',2.5,'DisplayName',['QR ' num2str(Qplot(3)*100) 'th ']);
ylim([-10 8]);
rr=recessionplot;
hleg = legend([l2 l3 l5 l6],'Orientation','Horizontal','Location','South','interpreter','Latex');legend boxoff;
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
print('-dpdf',fig,[slides_folder model '_Quantiles_MS_QR_recent'],'-bestfit');

%% SYNTHETIC DECEMBER 2020
t_oct_2019 = find(datenum('2019-Oct',dataformat)==dates);

for_s = 2; % how many months ahead


% Read FF AND MF FROM VINTAGE FILES
DATA_SYN = readtable('Syn_MF_FF.xlsx','Sheet','data');
FF_COVID_short = DATA_SYN.FF(end-11-2:end-2);
MF_COVID_short = DATA_SYN.MF_US(end-11-2:end-2);
FF_COVID_add = DATA_SYN.FF(end-1:end);
MF_COVID_add = DATA_SYN.MF_US(end-1:end);
% TR_COVID = DATA_APR2.TREND_US(556:567);
% FF_COVID = db_full.FF.data(end-11:end);
% MF_COVID = db_full.MF.data(end-11:end);
TR_COVID = db_full.TRENDH.data(end-11:end);
TR_COVID = [TR_COVID;TR_COVID(end);TR_COVID(end)];


% Simulation options
nDraws = 10;
nParamDraws = options.N;

% Allocate matrices
dY_erg_1_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_2_matcovid = NaN(1,nDraws,nParamDraws);
dY_erg_matcovid = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, "Covid-19 distribution: Simulating for December (Synthetic)...");

tic
for i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
    %===================
    % MAP PARAMETERS
    %===================
    
    % Transition prob parameters
    a12 = myparams(strcmp('a12',pnames),i);
    b12 = myparams(strcmp('b12',pnames),i);
    c12 = myparams(strcmp('c12',pnames),i);
    a21 = myparams(strcmp('a21',pnames),i);
    b21 = myparams(strcmp('b21',pnames),i);
    c21 = myparams(strcmp('c21',pnames),i);
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);
%     ax.XTick = datenum(dates_full(sd:numticks:ed_full));
% datetick('x','yyyy','keepticks')
    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    % VAR equation, FF
    c_1_1  = myparams(strcmp('c_1_1',pnames),i); % constant
    a1_1_1  = myparams(strcmp('a1_1_1',pnames),i); % FF(-1)
    a1_1_2  = myparams(strcmp('a1_1_2',pnames),i); % MF(-1)
    a0_1_2  = myparams(strcmp('a0_1_2',pnames),i)*(-1); % MF
    s_1_1  = myparams(strcmp('s_1_1',pnames),i); % std dev
    % VAR equation, MF
    c_2_1  = myparams(strcmp('c_2_1',pnames),i); % constant
    a1_2_1  = myparams(strcmp('a1_2_1',pnames),i); % FF(-1)
    a1_2_2  = myparams(strcmp('a1_2_2',pnames),i); % MF(-1) 
    s_2_2  = myparams(strcmp('s_2_2',pnames),i); % std dev
    % Draw shocks for October 2020
    shocks_sim = randn(1,nDraws);

    FF_COVID_base = FF_COVID_short;
    MF_COVID_base = MF_COVID_short;
    FF_COVID_syn = FF_COVID_short;
    MF_COVID_syn = MF_COVID_short;
    
    for ww=1:for_s
        shock_mf = randn(1); shock_ff = randn(1);
        MF_temp = c_2_1 + a1_2_1*FF_COVID_base(end) + a1_2_2*MF_COVID_base(end) + s_2_2*shock_mf;
        FF_temp = c_1_1 - a0_1_2*MF_temp + a1_1_1*FF_COVID_base(end) + a1_1_2*MF_COVID_base(end) + s_1_1*shock_ff;
        MF_COVID_base = [MF_COVID_base;MF_temp];
        FF_COVID_base = [FF_COVID_base;FF_temp];
%         MF_temp_syn = min(MF_COVID_syn)*ww*0.25;
%         FF_temp_syn = max(FF_COVID_syn)*ww*0.25;
        MF_temp_syn = MF_COVID_add(ww);
        FF_temp_syn = FF_COVID_add(ww);
        MF_COVID_syn = [MF_COVID_syn;MF_temp_syn];
        FF_COVID_syn = [FF_COVID_syn;FF_temp_syn];
%         MF_COVID_syn = [MF_COVID_syn;MF_COVID_syn(end)-(abs(min(MF_COVID_syn))*0.25)];
%         FF_COVID_syn = [FF_COVID_syn;FF_COVID_syn(end)+(abs(max(FF_COVID_syn))*0.25)];
    end
        
   
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in October 2019
        
        % regime-2 ***
        st_lag_base = prob_reg2_new(t_oct_2019);
        st_lag_syn = prob_reg2_new(t_oct_2019);
        
        % Get transition probabilities from Nov-2019 to Oct-2020
        if const==1
            if normal
                % Transition probabilities when coefficients have normal prior
                p12_base = 1./(1+exp(a12-b12*(FF_COVID_base)-c12*(MF_COVID_base)));
                p21_base = 1./(1+exp(a21-b21*(FF_COVID_base)-c21*(MF_COVID_base)));
                p12_syn = 1./(1+exp(a12-b12*(FF_COVID_syn)-c12*(MF_COVID_syn)));
                p21_syn = 1./(1+exp(a21-b21*(FF_COVID_syn)-c21*(MF_COVID_syn)));
            else
                % Transition probabilities when coefficients have gamma prior
                p12_base = 1./(1+exp(a12-b12.*(FF_COVID_base)+c12.*(MF_COVID_base)));
                p21_base = 1./(1+exp(a21+b21.*(FF_COVID_base)-c21.*(MF_COVID_base)));
                p12_syn = 1./(1+exp(a12-b12.*(FF_COVID_syn)+c12.*(MF_COVID_syn)));
                p21_syn = 1./(1+exp(a21+b21.*(FF_COVID_syn)-c21.*(MF_COVID_syn)));
            end
        elseif const==0
            if normal
                % Transition probabilities when coefficients have normal prior
                p12_base = 1./(1+exp(-b12*(FF_COVID_base)-c12*(MF_COVID_base)));
                p21_base = 1./(1+exp(-b21*(FF_COVID_base)-c21*(MF_COVID_base)));
                p12_syn = 1./(1+exp(-b12*(FF_COVID_syn)-c12*(MF_COVID_syn)));
                p21_syn = 1./(1+exp(-b21*(FF_COVID_syn)-c21*(MF_COVID_syn)));
            else
                % Transition probabilities when coefficients have gamma prior
                p12_base = 1./(1+exp(-b12.*(FF_COVID_base)+c12.*(MF_COVID_base)));
                p21_base = 1./(1+exp(+b21.*(FF_COVID_base)-c21.*(MF_COVID_base)));
                p12_syn = 1./(1+exp(-b12.*(FF_COVID_syn)+c12.*(MF_COVID_syn)));
                p21_syn = 1./(1+exp(+b21.*(FF_COVID_syn)-c21.*(MF_COVID_syn)));
            end
        end
        
        % Compute the probability of remaining in a given regime
        p11_base = ones(12+for_s,1) - p12_base;
        p22_base = ones(12+for_s,1) - p21_base;
        p11_syn = ones(12+for_s,1) - p12_syn;
        p22_syn = ones(12+for_s,1) - p21_syn;
        
        % Draw Markov Chain from Nov-2019 to Oct-2020
        for tt2 = 1:12+for_s
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag_base == 0 % started in good regime
                if udraw > p11_base(tt2)
                    st_base = 1; % switch from good to bad
                else
                    st_base = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22_base(tt2)
                    st_base = 0; % switch from bad to good
                else
                    st_base = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag_base = st_base;
            
            if st_lag_syn == 0 % started in good regime
                if udraw > p11_syn(tt2)
                    st_syn = 1; % switch from good to bad
                else
                    st_syn = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22_syn(tt2)
                    st_syn = 0; % switch from bad to good
                else
                    st_syn = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag_syn = st_syn;            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use_base = st_base;
        st_use_syn = st_syn;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good_base = 1 - st_use_base;
        IND_good_syn = 1 - st_use_syn;
        
        % Good regime in period t
        dY_erg_1_matcovid_base(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID_base(tt2)...
            - a0_3_2_sync_1*MF_COVID_base(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid_base(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID_base(tt2)...
            - a0_3_2_sync_2*MF_COVID_base(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid_base(1,nsim,i) = dY_erg_1_matcovid_base(1,nsim,i)*IND_good_base + dY_erg_2_matcovid_base(1,nsim,i)*(1-IND_good_base);
        
        St_sim_COV_temp_base(1,nsim) = IND_good_base;
        
        % Good regime in period t
        dY_erg_1_matcovid_syn(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_COVID_syn(tt2)...
            - a0_3_2_sync_1*MF_COVID_syn(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        % Bad regime in period t
        dY_erg_2_matcovid_syn(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_COVID_syn(tt2)...
            - a0_3_2_sync_2*MF_COVID_syn(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_COVID(tt2);
        
        dY_erg_matcovid_syn(1,nsim,i) = dY_erg_1_matcovid_syn(1,nsim,i)*IND_good_syn + dY_erg_2_matcovid_syn(1,nsim,i)*(1-IND_good_syn);
        
        St_sim_COV_temp_syn(1,nsim) = IND_good_syn;        
        
    end
    
     St_sim_COV_base(:,i) = mean(St_sim_COV_temp_base);
     St_sim_COV_syn(:,i) = mean(St_sim_COV_temp_syn);    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);

% Ergodic probability of regime 1;
p_reg1_sim_dec_base =  mean(St_sim_COV_base);
p_reg2_sim_dec_base =  1-p_reg1_sim_dec_base;
p_reg1_sim_dec_syn =  mean(St_sim_COV_syn);
p_reg2_sim_dec_syn =  1-p_reg1_sim_dec_syn;

% Average out simulations -> returns a TxmParamDraws matrix
y_erg_barcovid_dec_base = (reshape(dY_erg_matcovid_base,1,nDraws*nParamDraws));
y_erg_barcovid_dec_syn = (reshape(dY_erg_matcovid_syn,1,nDraws*nParamDraws));
y_erg_1_barcovid_dec_syn = (reshape(dY_erg_1_matcovid_syn,1,nDraws*nParamDraws));
y_erg_2_barcovid_dec_syn = (reshape(dY_erg_2_matcovid_syn,1,nDraws*nParamDraws));


%% SYNTHETIC DECEMBER DENSITIES with MS
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------
colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

colors2 = cbrewer('div', 'RdYlGn', 18);

[pdf_dec_base,xi_dec_base]=ksdensity(y_erg_barcovid_dec_base(1,:));
[cdf_dec_base,xic_dec_base]=ksdensity(y_erg_barcovid_dec_base(1,:),'function','cdf');
[pdf_dec_reg1,xi_dec_reg1]=ksdensity(y_erg_1_barcovid_dec_syn(1,:));
[cdf_dec_reg1,xic_dec_reg1]=ksdensity(y_erg_1_barcovid_dec_syn(1,:),'function','cdf');
[pdf_dec_reg2,xi_dec_reg2]=ksdensity(y_erg_2_barcovid_dec_syn(1,:));
[cdf_dec_reg2,xic_dec_reg2]=ksdensity(y_erg_2_barcovid_dec_syn(1,:),'function','cdf');

[~,i25_ms_dec_base] = min(abs(cdf_dec_base-0.25));
[~,i25_ms_dec_syn] = min(abs(cdf_dec_syn-0.25));
[~,i25_ms_dec_reg1] = min(abs(cdf_dec_reg1-0.25));
[~,i25_ms_dec_reg2] = min(abs(cdf_dec_reg2-0.25));
staff = 3.13;
s_wave = -0.69;

ymax = max([max(pdf_dec_base)*1.05,max(pdf_dec_syn)*1.05,max(pdf_dec_reg1)*1.05,max(pdf_dec_reg2)*1.05]);

dcolor = [mcolor;ocolor];
left_color = mcolor;
right_color = ocolor;
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -12;
xmax = 12;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
% plot(ax1,xi_dec_base,pdf_dec_base,'LineWidth', 2,'Color',colors2(16,:),'DisplayName',['2020-Dec (Baseline)']);
% plot(ax1,xi_dec_syn,pdf_dec_syn,'LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Dec (Synthetic)']);
plot(ax1,xi_dec_syn,pdf_dec_syn,'LineWidth', 3,'Color',colors2(6,:),'DisplayName',['2020-Dec']);
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
% plot(ax1,xi_dec_base(i25_ms_dec_base),pdf_dec_base(i25_ms_dec_base),'d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_base(i25_ms_dec_base),1)) '\% (Dec. Baseline)']);
% plot(ax1,xi_dec_syn(i25_ms_dec_syn),pdf_dec_syn(i25_ms_dec_syn),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_syn(i25_ms_dec_syn),1)) '\% (Dec. Synthetic)']);
plot(ax1,xi_dec_syn(i25_ms_dec_syn),pdf_dec_syn(i25_ms_dec_syn),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_syn(i25_ms_dec_syn),1)) '\%']);
plot([staff staff],[0 ymax],'Color',colors2(14,:),'LineWidth',1,'DisplayName',['Staff forecast = ' num2str(round(staff,1)) '\%']);
plot([s_wave s_wave],[0 ymax],'Color',colors2(2,:),'LineWidth',1,'DisplayName',['Second Waves = ' num2str(round(s_wave,1)) '\%']);
axis tight
yl = ylim;
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('MS-VAR: Densitites for GDP growth over the next 12 months (U.S.)','FontSize',16','Interpreter','Latex');
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
    print('-dpdf',fig,[fig_folder 'DensityMS_Dec_Syn'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Aug_Oct_foreign'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Oct_foreign'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end

%% BAD AND NORMAL REGIME
fig=figure;
set(fig,'defaultAxesColorOrder',[right_color]);
xmin = -12;
xmax = 12;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
% PDF Plot
%yyaxis left
hold on;
% plot(ax1,xi_dec_base,pdf_dec_base,'LineWidth', 2,'Color',colors2(16,:),'DisplayName',['2020-Dec (Baseline)']);
% plot(ax1,xi_dec_syn,pdf_dec_syn,'LineWidth', 2,'Color',colors2(6,:),'DisplayName',['2020-Dec (Synthetic)']);
plot(ax1,xi_dec_reg1,pdf_dec_reg1,'LineWidth', 3,'Color',colors2(6,:),'DisplayName',['2020-Dec (normal regime)']);
plot(ax1,xi_dec_reg2,pdf_dec_reg2,'LineWidth', 3,'Color',colors2(16,:),'DisplayName',['2020-Dec (bad regime)']);
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
% plot(ax1,xi_dec_base(i25_ms_dec_base),pdf_dec_base(i25_ms_dec_base),'d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_base(i25_ms_dec_base),1)) '\% (Dec. Baseline)']);
% plot(ax1,xi_dec_syn(i25_ms_dec_syn),pdf_dec_syn(i25_ms_dec_syn),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_syn(i25_ms_dec_syn),1)) '\% (Dec. Synthetic)']);
plot(ax1,xi_dec_reg1(i25_ms_dec_reg1),pdf_dec_reg1(i25_ms_dec_reg1),'d','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_reg1(i25_ms_dec_reg1),1)) '\%  (normal regime)']);
plot(ax1,xi_dec_reg2(i25_ms_dec_reg2),pdf_dec_reg2(i25_ms_dec_reg2),'d','MarkerFaceColor',colors2(16,:),'MarkerEdgeColor','black','Markersize',10,'DisplayName',['25th Quantile = ' num2str(round(xi_dec_reg2(i25_ms_dec_reg2),1)) '\%  (bad regime)']);
plot([staff staff],[0 ymax],'Color',colors2(14,:),'LineWidth',1,'DisplayName',['Staff forecast = ' num2str(round(staff,1)) '\%']);
plot([s_wave s_wave],[0 ymax],'Color',colors2(2,:),'LineWidth',1,'DisplayName',['Second Waves = ' num2str(round(s_wave,1)) '\%']);
axis tight
yl = ylim;
hL=legend(ax1,'Location','NorthWest');
set(hL,'interpreter','Latex')
legend boxoff
title('MS-VAR: Densitites for GDP growth over the next 12 months (U.S.)','FontSize',16','Interpreter','Latex');
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
    print('-dpdf',fig,[fig_folder 'DensityMS_Dec_Regs'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Aug_Oct_foreign'],'-bestfit');
%     print('-dpdf',fig,[fig_folder 'DensityMS_Mar_Oct_foreign'],'-bestfit');
    %saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end
