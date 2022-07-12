%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
dir_or_it = 1; % 1 = direct, 2 = iterated
saveit    = 0; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

% Set DGP options
burn = 500; % 500 burn, 500 keep
Tsim = 1200;
tt=Tsim+burn; % total time periods

% Pick a particular simulation
draw =10;


% Date vector for simulated data (After burn-in)
dates_full = datenum(1900,1:Tsim,1);

% Data vintage, sample and country selection
start_date   = datestr(dates_full(1),'yyyy-mmm');
end_date     = datestr(dates_full(end-12),'yyyy-mmm');
end_date_db  = datestr(dates_full(end),'yyyy-mmm');

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

% Create date formats for plotting
inputformat = 'yyyy-MMM';
dataformat  = 'yyyy-mmm';
start_plot    = datestr(dates_full(1+nlags),'yyyy-mmm');
end_plot      = end_date_db;
end_plot_full = end_date;

% Vector of dates for the full sample
%dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

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

%% Simulate DGP
% FF = c_1_1 + a0_1_2*MF + + a1_1_1*FF(-1) + a1_1_2*MF(-1)
% FF s.d. = s_1_1
% MF = c_2_1 + + a1_2_1*FF(-1) + a1_2_2*MF(-1)
% MF s.d. = s_2_2
% GDP, normal = c_3_1_sync_1 + a0_3_1_sync_1*FF + a0_3_2_sync_1*MF
% GDP, normal s.d. = s_3_3_sync_1
% GDP, bad = c_3_1_sync_2 + a0_3_1_sync_2*FF + a0_3_2_sync_2*MF
% GDP, bad s.d. = s_3_3_sync_2
% 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
% 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'

dgp.c_1_1 = 0; dgp.a0_1_1 = 0; dgp.a0_1_2 = -0.15; dgp.a0_1_3 = 0; dgp.a1_1_1 =  0.75; dgp.a1_1_2 = -0.05; dgp.a1_1_3=0; dgp.s_1_1 = 1;
dgp.c_2_1 = 0; dgp.a0_2_1 = 0; dgp.a0_2_2 = 0;    dgp.a0_2_3 = 0; dgp.a1_2_1 = -0.1; dgp.a1_2_2 = 0.85;  dgp.a1_2_3=0; dgp.s_2_2 = 1;
dgp.c_3_1_sync_1 =  0.5; dgp.a0_3_1_sync_1 = -0.1; dgp.a0_3_2_sync_1 = 0.1; dgp.a0_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = 0.5;
dgp.c_3_1_sync_2 = -0.5; dgp.a0_3_1_sync_2 = -0.3; dgp.a0_3_2_sync_2 = 0.3; dgp.a0_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = 2;
dgp.a12 = 2; dgp.a21 = 2; dgp.b12 = 0.5; dgp.b21 = 0.5; dgp.c12 = 0.5; dgp.c21 = 0.5;


%% SIMULATE N DRAWS OF THE GDP


% 2 = BAD, 1 = NORMAL
burn = 500; % 500 burn, 500 keep
tt   = Tsim+burn; % time periods
rr   = 600; % draws
f_draw = zeros(rr,tt);m_draw = zeros(rr,tt);
p12_draw = NaN(tt-burn,1); p21_draw = NaN(tt-burn,1);
y_mat = NaN(rr,tt-burn);

st_mat = NaN(rr,tt); st_mat(:,1) = ones(rr,1);

for dd=1:rr
    for jj=2:tt
        eta1 = randn(2,1); % financial and macro shocks
        m_draw(dd,jj) = dgp.c_2_1                            + dgp.a1_2_1*f_draw(dd,jj-1) + dgp.a1_2_2*m_draw(dd,jj-1) + dgp.s_2_2*eta1(2,1);
        f_draw(dd,jj) = dgp.c_1_1 + dgp.a0_1_2*m_draw(dd,jj) + dgp.a1_1_1*f_draw(dd,jj-1) + dgp.a1_1_2*m_draw(dd,jj-1) + dgp.s_1_1*eta1(1,1);
        
        p12 = 1/(1+exp(dgp.a12-dgp.b12*(f_draw(dd,jj))+dgp.c12*(m_draw(dd,jj)))); % 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))' %USE THIS!
        p21 = 1/(1+exp(dgp.a21+dgp.b21*(f_draw(dd,jj))-dgp.c21*(m_draw(dd,jj)))); % 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))' %USE THIS!
        
        p11 = 1 - p12; % probability of remaining in normal
        p22 = 1 - p21; % probability of remaining in bad
        
        
        udraw = rand(1);
        
        if st_mat(dd,jj-1) == 1 % started in good regime
            if udraw > p11
                st_temp = 2; % switch from good to bad
            else
                st_temp = 1; % don't switch and remain in good
            end
            
        else % start in bad regime
            if udraw > p22
                st_temp = 1; % switch from bad to good
            else
                st_temp = 2; % don't switch and remain in bad
            end
            
        end
        
        eta2 = randn(1,1); % GDP shock
        
        if st_temp==2
            y_draw = dgp.c_3_1_sync_2 + dgp.a0_3_1_sync_2*f_draw(dd,jj) + dgp.a0_3_2_sync_2*m_draw(dd,jj) + dgp.s_3_3_sync_2*eta2;
        else
            y_draw = dgp.c_3_1_sync_1 + dgp.a0_3_1_sync_1*f_draw(dd,jj) + dgp.a0_3_2_sync_1*m_draw(dd,jj) + dgp.s_3_3_sync_1*eta2;
        end
        
        %         if dd==rr % randomly assign a state for time jj
        %             ind_st = round(rand(1)*(length(st_temp))+0.5,0);
        %             st_draw(jj) = st_temp(ind_st);
        %         end
        st_mat(dd,jj) = st_temp;
        if jj>burn
            y_mat(dd,jj-burn) = y_draw;
            %             f_mat(dd,jj-burn) = f_draw(dd,jj); m_mat(dd,jj-burn) = m_draw(dd,jj);
            p12_draw(dd,jj-burn) = p12; p21_draw(dd,jj-burn) = p21;
        end
    end
end

st_mat = st_mat(:,burn+1:end);
f_mat =  f_draw(:,burn+1:end);
m_mat =  m_draw(:,burn+1:end);

%
% 12-months ahead
y_mat_fut = NaN(size(y_mat,1),size(y_mat,2)-12);
for hh=1:size(y_mat,1)
    for ww=1:size(y_mat,2)-12
        y_mat_fut(hh,ww) = mean(y_mat(hh,ww+1:ww+12));
    end
end


% Compute percentiles
dY_25 = prctile(y_mat,25)'; dY_25_fut = prctile(y_mat_fut,25)';
dY_75 = prctile(y_mat,75)'; dY_75_fut = prctile(y_mat_fut,75)';
dY_10 = prctile(y_mat,10)'; dY_10_fut = prctile(y_mat_fut,10)';
dY_90 = prctile(y_mat,90)'; dY_90_fut = prctile(y_mat_fut,90)';

ds_50 = prctile(st_mat,50)';
ds_10 = prctile(st_mat,10)';
ds_90 = prctile(st_mat,90)';

Res_iterated.y_mat = y_mat; Res_iterated.FF = f_mat; Res_iterated.MF = m_mat;
Res_iterated.dY_25 = dY_25; Res_iterated.dY_75 = dY_75;
Res_iterated.dY_10 = dY_10; Res_iterated.dY_90 = dY_90;

Res_direct.y_mat = y_mat_fut; Res_direct.FF = f_mat(1:end-12); Res_direct.MF = m_mat(1:end-12);
Res_direct.dY_25 = dY_25_fut; Res_direct.dY_75 = dY_75_fut;
Res_direct.dY_10 = dY_10_fut; Res_direct.dY_90 = dY_90_fut;



%% Prepare Data

%load data
Res.Res_iterated = Res_iterated; Res.Res_direct = Res_direct;

% Number of simulated paths
n_draws = size(Res.Res_iterated.y_mat,1);


[db, db_full, tex] = data_sim(Res,draw,start_date);
nlags=1;

% Collect MF, FF and trend (estimation sample)
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);

% Collect MF, FF and trend (full sample)
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);

if dir_or_it ==1
    varlist={'FF','MF','GDPGH'}';
elseif dir_or_it ==2
    varlist={'FF','MF','GDPG'}';
end

%% set up Markov chains

switches = {'c(3)','a0(3)','s(3)'};
model = 'SVAR_StateDependent_2Regimes';

temp = {'direct_','iterated_'};


% Model label
model_label = ['MC' temp{dir_or_it} '_Draw_' num2str(draw) '_Tsim_' num2str(Tsim) '_const' num2str(const) '_normal' num2str(normal)];

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
if dir_or_it ==1 % Direct forecast
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
elseif dir_or_it ==2 % Iterated forecast
    lin_restr={
        % first equation or "financial factor" equation
        %----------------------------------
        'a0(1,GDPG)=0'
        'a1(1,GDPG)=0'
        %'a1(1,MF)=0'
        % second equation or "macroeconomic factor" equation
        %-----------------------------------
        'a0(2,GDPG)=0'
        'a0(2,FF)=0'
        'a1(2,GDPG)=0'
        % third equation or "GDP Growth" equation
        %----------------------------------
        'a1(3,GDPG)=0'
        'a1(3,FF)=0'
        'a1(3,MF)=0'
        };
end

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


%% Store estimates

%options = optimset('TolX',1e-12);
if exist([mcmc_folder 'sv' model_label '.mat'],'file')==2
    load([mcmc_folder 'sv' model_label '.mat'])
else
    %sv=estimate(sv0,db,{'1973M1','2019M5'},prior,restrictions,'fmincon',false,'optimset',options);
    sv=estimate(sv0,db,{[datestr(dates_full(1),'yyyy') 'M' datestr(dates_full(1),'mm')]...
        ,[datestr(dates_full(end-12),'yyyy') 'M' datestr(dates_full(end-12),'mm')]},prior,restrictions,'fmincon');
    if save_mcmc==1
        save([mcmc_folder 'sv' model_label '.mat'],'sv')
    end
end


pmode=posterior_mode(sv);
print_structural_form(sv);
% print_solution(sv)


%%
%==========================================================================
% Posterior sampling
%==========================================================================

pnames=fieldnames(pmode);

a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);

options=struct();
options.alpha=0.234;
options.thin=5;
options.burnin=10^3;
options.N=1*10^4;
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
%% Compute Ergodic Fitted Values

[Resids,Fits]=residuals(sv);
fit1 = Fits.GDPGH.data(:,1);
fit2 = Fits.GDPGH.data(:,2);

% Smoothed and Filtered probabilities
[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;

p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;


%% Recursive updating of s(t) from simulations of Markov-chain

% Number of draws for simulting predictive density
nDraws = 20000;

% Set a 0.9 threshold to switch from good to bad regime
prob_reg2_new = zeros(1,length(p_reg2));
prob_reg2_new(p_reg2(sd:end)>0.9) = 1;

Tobs = length(prob_reg2_new);
% Draw shocks for simulation
rng(123);
shocks_sim = randn(length(prob_reg2_new),nDraws);

% Collect MF, FF and trend
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);
%TR = db.TRENDH.data(nlags+1:end);
st0_dgp = zeros(1,nDraws);
st_mat_dgp = NaN(Tobs,nDraws);
DGP_sim    = NaN(Tobs,nDraws);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws
    
    for tt= 1:Tobs
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities
        
        % Initialize state in period-0
        if tt==1
            st_lag  = st0_dgp(1,nsim);
        end
        
                
        % Get transition probabilities from t-1:t
        if const==1
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF(tt))-pmode.c12*(MF(tt))));
                p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF(tt))-pmode.c21*(MF(tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                %                     p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                %                     p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
                
                p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt))+pmode.c12.*(MF(tt))));
                p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt))-pmode.c21.*(MF(tt))));
                
            end
        elseif const==0
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(-pmode.b12*(FF(tt))-pmode.c12*(MF(tt))));
                p21 = 1./(1+exp(-pmode.b21*(FF(tt))-pmode.c21*(MF(tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(-pmode.b12.*(FF(tt))+pmode.c12.*(MF(tt))));
                p21 = 1./(1+exp(+pmode.b21.*(FF(tt))-pmode.c21.*(MF(tt))));
            end
        end
        % Compute the probability of remaining in a given regime
        p11 = 1 - p12; %ones(12,1) - p12;
        p22 = 1 - p21; %ones(12,1) - p21;
        
        
        
        % Draw the Markov Chain for period t|t-1
        % Note: Assume that in the MS-VAR econometrician has access to
        % s(t-1) in real time. In our application we need to simulate
        % s(t)|s(t-11)
        
        udraw = rand(1);
        
        if st_lag == 0 % started in good regime
            if udraw > p11
                st = 1; % switch from good to bad
            else
                st = 0; % don't switch and remain in good
            end
            
        else % start in bad regime
            
            if udraw > p22
                st = 0; % switch from bad to good
            else
                st = 1; % don't switch and remain in bad
            end
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_lag = st;
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        %*************************************************************************************
        % SIMULATE ESTIMATED MODEL: y(t+12) = c(st) + a(st)*ff(t) + b(st)*mf(t) + d(st)*sigma(t)
        %*************************************************************************************
        
        % Good regime in period t
        dY_sim_1(tt,nsim) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF(tt)...
            - pmode.a0_3_2_sync_1*MF(tt) + pmode.s_3_3_sync_1*shocks_sim(tt,nsim));
        
        % Bad regime in period t
        dY_sim_2(tt,nsim) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF(tt)...
            - pmode.a0_3_2_sync_2*MF(tt)+ pmode.s_3_3_sync_2*shocks_sim(tt,nsim));
        
        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_good + dY_sim_2(tt,nsim)*(1-IND_good);
        
        St_sim(tt,nsim) = IND_good;
        
        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt);
        
        
        %*************************************************************************************
        % SIMULATE  DGP: y(t) = c(st) + a(st)*ff(t) + b(st)*mf(t) + d(st)*sigma(t)
        %*************************************************************************************
        
        if tt==1
            st_lag_dgp  = st0_dgp(1,nsim);
        else
            st_lag_dgp  = st_mat_dgp(tt-1,nsim);
        end
        


        p12_dgp = 1/(1+exp(dgp.a12-dgp.b12*(FF(tt))+dgp.c12*(MF(tt)))); % 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))' %USE THIS!
        p21_dgp = 1/(1+exp(dgp.a21+dgp.b21*(FF(tt))-dgp.c21*(MF(tt)))); % 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))' %USE THIS!
        
        p11_dgp = 1 - p12_dgp; % probability of remaining in normal
        p22_dgp = 1 - p21_dgp; % probability of remaining in bad
        
        if st_lag == 1 % started in good regime
            if udraw > p11_dgp
                st_temp = 2; % switch from good to bad
            else
                st_temp = 1; % don't switch and remain in good
            end
            
        else % start in bad regime
            if udraw > p22_dgp
                st_temp = 1; % switch from bad to good
            else
                st_temp = 2; % don't switch and remain in bad
            end
            
        end
        
        st_mat_dgp(tt,nsim) = st_temp;
        
        if st_temp==2
            DGP_sim(tt,nsim) = dgp.c_3_1_sync_2 + dgp.a0_3_1_sync_2*FF(tt) + dgp.a0_3_2_sync_2*MF(tt) + dgp.s_3_3_sync_2*shocks_sim(tt,nsim);
        else
            DGP_sim(tt,nsim) = dgp.c_3_1_sync_1 + dgp.a0_3_1_sync_1*FF(tt) + dgp.a0_3_2_sync_1*MF(tt) + dgp.s_3_3_sync_1*shocks_sim(tt,nsim);
        end
        
    end
    
    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end


% transform GDP to 12-month ahead averages
y_mat_fut2 = NaN(Tobs,nDraws);
for ww=1:Tobs-12
    y_mat_fut2(ww,:) = mean(DGP_sim(ww+1:ww+12,:));
end

%%
% Ergodic probability of regime 1;
p_reg1_sim(:,1) =  sum(St_sim,2)/size(sum(St_sim),2);
p_reg2_sim(:,1) =  sum((1-St_sim),2)/size(sum(St_sim),2);


% Compute percentiles of direct forecasting model
dY_sim_25 = prctile(dY_sim',25)';
dY_sim_75 = prctile(dY_sim',75)';
dY_sim_10 = prctile(dY_sim',10)';
dY_sim_90 = prctile(dY_sim',90)';

% Compute Regime specific means
dY_sim1_mean = mean(dY_sim_1,2);
dY_sim2_mean = mean(dY_sim_2,2);

% Compute percentiles of iterated model
dY_25_fut2 = prctile(y_mat_fut2',25)';
dY_75_fut2 = prctile(y_mat_fut2',75)';
dY_10_fut2 = prctile(y_mat_fut2',10)';
dY_90_fut2 = prctile(y_mat_fut2',90)';

%% Figure: Predictive Density DGP vs MSVAR

fig=figure; clf;
hold on
% Plot percentiles from MS
l1=plot(dY_10_fut2(1:200),':','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','DGP 10th');
l2=plot(dY_90_fut2(1:200),':','Color',colors(55,:),'LineWidth', 1.5,'DisplayName','DGP 90th');
l3=plot(dY_sim_10(1:200),'Color',colors(5,:),'LineWidth', 2,'DisplayName','MS-VAR 10th');
l4=plot(dY_sim_90(1:200),'Color',colors(55,:),'LineWidth', 2,'DisplayName','MS-VAR 90th');
l5=plot(y_mat_fut(draw,2:200),'k-.','LineWidth',2,'DisplayName','Data');
hleg = legend([l1 l2 l3 l4 l5],'Orientation','Horizontal','Location','S','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('1-Year-Ahead GDP Growth (\%)','interpreter','Latex','fontsize',10)
axis tight
ylim([-4,1.5]);
ax=gca;
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
title(['Predictive Density: DGP vs MS-VAR, draw = ' num2str(draw) ],'FontSize',16','Interpreter','Latex');
tightfig;


print('-dpdf',fig,['PredictiveDensity_DGP_MSVAR_draw' num2str(draw)],'-bestfit');

%% DEPRECATED


% %% SIMULATE 1 DRAW OF THE GDP x 1000 times
%
%
% % 2 = BAD, 1 = NORMAL
% burn = 0; % 500 burn, 500 keep
% tt   = Tsim; % time periods
% rr   = 20000; % draws
% %f_draw = zeros(rr,tt);m_draw = zeros(rr,tt);
% p12_draw = NaN(tt-burn,1); p21_draw = NaN(tt-burn,1);
% y_mat2 = NaN(rr,tt-burn);
%
%
% rng(123);
%
% st_mat = NaN(rr,tt); st_mat(:,1) = ones(rr,1);
%
% for dd=1:rr
%     for jj=2:tt
%
%     p12 = 1/(1+exp(dgp.a12-dgp.b12*(f_mat(draw,jj))+dgp.c12*(m_mat(draw,jj)))); % 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))' %USE THIS!
%     p21 = 1/(1+exp(dgp.a21+dgp.b21*(f_mat(draw,jj))-dgp.c21*(m_mat(draw,jj)))); % 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))' %USE THIS!
%
%     p11 = 1 - p12; % probability of remaining in normal
%     p22 = 1 - p21; % probability of remaining in bad
%
%
%         udraw = rand(1);
%
%         if st_mat(dd,jj-1) == 1 % started in good regime
%             if udraw > p11
%                 st_temp = 2; % switch from good to bad
%             else
%                 st_temp = 1; % don't switch and remain in good
%             end
%
%         else % start in bad regime
%             if udraw > p22
%                 st_temp = 1; % switch from bad to good
%             else
%                 st_temp = 2; % don't switch and remain in bad
%             end
%
%         end
%
%         eta2 = shocks_sim(jj,:); % GDP shock
%
%         if st_temp==2
%             y_draw = dgp.c_3_1_sync_2 + dgp.a0_3_1_sync_2*f_mat(draw,jj) + dgp.a0_3_2_sync_2*m_mat(draw,jj) + dgp.s_3_3_sync_2*eta2;
%         else
%             y_draw = dgp.c_3_1_sync_1 + dgp.a0_3_1_sync_1*f_mat(draw,jj) + dgp.a0_3_2_sync_1*m_mat(draw,jj) + dgp.s_3_3_sync_1*eta2;
%         end
%
% %         if dd==rr % randomly assign a state for time jj
% %             ind_st = round(rand(1)*(length(st_temp))+0.5,0);
% %             st_draw(jj) = st_temp(ind_st);
% %         end
%
% % Store objects:
%         st_mat(dd,jj) = st_temp;
%         y_mat2(dd,jj) = y_draw;
%         p12_draw(dd,jj) = p12;
%         p21_draw(dd,jj) = p21;
%     end
% end
%
%
% %
% % 12-months ahead
% y_mat_fut2 = NaN(size(y_mat2,1),size(y_mat2,2)-12);
% for hh=1:size(y_mat2,1)
%     for ww=1:size(y_mat2,2)-12
%         y_mat_fut2(hh,ww) = mean(y_mat2(hh,ww+1:ww+12));
%     end
% end


% %%
% dY_25_fut2 = prctile(y_mat_fut2,25)';
% dY_75_fut2 = prctile(y_mat_fut2,75)';
