%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all; clc;

rng(777); % Fixed seed for replication purposes

% Important paths
addpath('../RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
saveit    = 1; % 1 = save graphs
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
FF = db_full.FF.data;
MF = db_full.MF.data;
trend_in = 0; % 0 without, 1 in
if trend_in ==0
    GDP = db_full.GDPG.data;
elseif trend_in ==1
    GDP = db_full.GDPG.data+db_full.TRENDH.data; %TRENDH IS ACTUALLY THE MONTLHY AND NOT THE 1-YEAR-AHEAD (FIX THIS)
end

REC_US_temp = readtable('rec_us.xlsx','Sheet','Sheet1');
REC_US = table2array(REC_US_temp(:,2));
load('p_sim');

Reg_variable = 2; % 1, for 'p_reg2_filtered'
                  % 2, for 'REC_US'
                  % 3, for GDP lower than 1 s.d.
                  
Yraw = [MF,FF,GDP]; GDPraw = GDP; MFraw = MF; FFraw = FF;

reg_low = zeros(size(GDP));
reg_low(GDP<((std(GDP)*(-1))+mean(GDP)))=1;

REC_US = REC_US(nlags+1:end,1);
p_bad = false(size(REC_US));
p_bad(REC_US==1) = true; % If 1, recession

p_norm = false(size(p_bad));
p_norm(p_bad==true) = false;
p_norm(p_bad==false) = true; % If 1, normal/good regime


[Traw, M] = size(Yraw);
Ylag = mlag2(Yraw,nlags);

MFlag =  mlag2(MFraw,nlags);
FFlag =  mlag2(FFraw,nlags);

Y_s = GDPraw(nlags+1:Traw,:);
MF_s = MFraw(nlags+1:Traw,:);
FF_s = FFraw(nlags+1:Traw,:);

X = [Ylag(nlags+1:Traw,:) ones(Traw-nlags,1)];
X_s = [MFlag(nlags+1:Traw,:) FFlag(nlags+1:Traw,:) ones(Traw-nlags,1)];
X_MF_s = [MFlag(nlags+1:Traw,:) ones(Traw-nlags,1)];
X_FF_s = [FFlag(nlags+1:Traw,:) ones(Traw-nlags,1)];
[~,K] = size(X);
Y = Yraw(nlags+1:Traw,:);
T = Traw - nlags;

%% VAR - OLS COEFFICIENTS (BAD REGIME)

% Reduced form: Y_t = A_OLS*X_t + e_t
% Structural form: A0_OLS*Y_t = Al_OLS*X_t + inv(A0_OLS)*u_t
% OBS: X_t = [Y_{t-1} constant]
%      Al_OLS = inv(A0_OLS)*A_OLS(1:nvar,:)

% VAR

Y_bad = Y(p_bad,:); T_bad = size(Y_bad,1)+nlags;
X_bad = X(p_bad,:);
Y_s_bad = Y_s(p_bad,:); X_s_bad = X_s(p_bad,:);

A_OLS_bad = (X_bad'*X_bad)\(X_bad'*Y_bad); % This is the matrix of regression coefficients
a_OLS_bad = A_OLS_bad(:);         % This is the vector of coefficients
% that a_OLS = vec(A_OLS)
SSE_bad = (Y_bad - X_bad*A_OLS_bad)'*(Y_bad - X_bad*A_OLS_bad);
SIGMA_OLS_bad = SSE_bad./(T_bad-K);
RESID_OLS_bad = (Y_bad - X_bad*A_OLS_bad); % residuals
VCV_OLS_bad = cov(RESID_OLS_bad,1);
A0_OLS_bad = (chol(SIGMA_OLS_bad,'lower'))\eye(size(chol(SIGMA_OLS_bad,'lower'),2));
A1_OLS_bad =  ((A0_OLS_bad')\(A_OLS_bad'))';

% Simple regression

A_S_bad = (X_s_bad'*X_s_bad)\(X_s_bad'*Y_s_bad); % This is the matrix of regression coefficients
a_S_bad = A_S_bad(:);         % This is the vector of coefficients
% that a_S = vec(A_S)
SSE_bad = (Y_s_bad - X_s_bad*A_S_bad)'*(Y_s_bad - X_s_bad*A_S_bad);
SIGMA_S_bad = SSE_bad./(T_bad-K);
RESID_S_bad = (Y_s_bad - X_s_bad*A_S_bad); % residuals
STD_S_bad = std(RESID_S_bad);



%% VAR - OLS COEFFICIENTS (NORMAL REGIME)

Y_norm = Y(p_norm,:); T_norm = size(Y_norm,1)+nlags;
X_norm = X(p_norm,:);
Y_s_norm = Y_s(p_norm,:); X_s_norm = X_s(p_norm,:);

A_OLS_norm = (X_norm'*X_norm)\(X_norm'*Y_norm); % This is the matrix of regression coefficients
a_OLS_norm = A_OLS_norm(:);         % This is the vector of coefficients
% that a_OLS = vec(A_OLS)
SSE_norm = (Y_norm - X_norm*A_OLS_norm)'*(Y_norm - X_norm*A_OLS_norm);
SIGMA_OLS_norm = SSE_norm./(T_norm-K);
RESID_OLS_norm = (Y_norm - X_norm*A_OLS_norm); % residuals
VCV_OLS_norm = cov(RESID_OLS_norm,1);
A0_OLS_norm = (chol(SIGMA_OLS_norm,'lower'))\eye(size(chol(SIGMA_OLS_norm,'lower'),2));
A1_OLS_norm =  ((A0_OLS_norm')\(A_OLS_norm'))';

% Simple regression

A_S_norm = (X_s_norm'*X_s_norm)\(X_s_norm'*Y_s_norm); % This is the matrix of regression coefficients
a_S_norm = A_S_norm(:);         % This is the vector of coefficients
% that a_S = vec(A_S)
SSE_norm = (Y_s_norm - X_s_norm*A_S_norm)'*(Y_s_norm - X_s_norm*A_S_norm);
SIGMA_S_norm = SSE_norm./(T_norm-K);
RESID_S_norm = (Y_s_norm - X_s_norm*A_S_norm); % residuals
STD_S_norm = std(RESID_S_norm);

%% TOTAL

RESID_S_total = NaN(T,1);
RESID_S_total(p_bad,:) = RESID_S_bad;
RESID_S_total(p_norm,:) = RESID_S_norm;

%% AR(1) MF AND FF

A_MF = (X_MF_s'*X_MF_s)\(X_MF_s'*MF_s); % This is the matrix of regression coefficients
RESID_MF = (MF_s - X_MF_s*A_MF); % residuals
STD_MF = std(RESID_MF);
A_FF = (X_FF_s'*X_FF_s)\(X_FF_s'*FF_s); % This is the matrix of regression coefficients
RESID_FF = (FF_s - X_FF_s*A_FF); % residuals
STD_FF = std(RESID_FF);

%% ESTIMATE c_11, a_12 AND a_13

X_res_bad = [RESID_FF(p_bad,:) RESID_MF(p_bad,:)];
A_res_bad = (X_res_bad'*X_res_bad)\(X_res_bad'*RESID_S_bad); % This is the matrix of regression coefficients
X_res_norm = [RESID_FF(p_norm,:) RESID_MF(p_norm,:)];
A_res_norm = (X_res_norm'*X_res_norm)\(X_res_norm'*RESID_S_norm); % This is the matrix of regression coefficients

a12_2 = A_res_bad(1); a12_1 = A_res_norm(1); % This should've be negative
a13_2 = A_res_bad(2); a13_1 = A_res_norm(2); % This should've be negative
c11_2 = std(RESID_S_bad - X_res_bad*A_res_bad); c11_1 = std(RESID_S_norm-X_res_norm*A_res_norm);

%% SIM

mu_2 = A_S_bad(end); mu_1 = A_S_norm(end);
gamma_2 = A_S_bad(1); gamma_1 = A_S_norm(1);
beta_2 = A_S_bad(2); beta_1 = A_S_norm(2);
sig_2 = STD_S_bad; sig_1 = STD_S_norm;
% sig_2 = 20; sig_1 = 1;
b22 = A_FF(1); b33 = A_MF(1);


% 2 = BAD, 1 = NORMAL
burn = 500; % 500 burn, 500 keep
tt=T+burn; % time periods
rr = 600; % draws
f_draw = zeros(1,tt);m_draw = zeros(1,tt);y_draw = zeros(1,tt);
p12_draw = NaN(tt-burn,1); p21_draw = NaN(tt-burn,1);
y_mat = NaN(rr,tt-burn);
f_mat = NaN(1,tt-burn);
m_mat = NaN(1,tt-burn);
st_mat = NaN(rr,tt); st_mat(:,1) = ones(rr,1);
a12 = 2; b12 = 0.1; c12 = 0.1;
a21 = 2; b21 = 0.1; c21 = 0.1;
for jj=2:tt
    eta       = randn(3,1); % financial and macro shocks
    f_draw(jj) = b22*f_draw(jj-1)+eta(2,1);
    m_draw(jj) = b33*m_draw(jj-1)+eta(3,1);
    p12        = 1/(1+exp(a12-b12*(f_draw(jj))+c12*(m_draw(jj)))); % 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))' %USE THIS!
    p21        = 1/(1+exp(a21+b21*(f_draw(jj))-c21*(m_draw(jj)))); % 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))' %USE THIS!
    
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    st_temp = NaN(rr,1);
    for dd=1:rr
        
        udraw = rand(1);
        
        if st_mat(dd,jj-1) == 1 % started in good regime
            if udraw > p11
                st_temp(dd) = 2; % switch from good to bad
            else
                st_temp(dd) = 1; % don't switch and remain in good
            end
            
        else % start in bad regime
            if udraw > p22
                st_temp(dd) = 1; % switch from bad to good
            else
                st_temp(dd) = 2; % don't switch and remain in bad
            end
            
        end
        
        eta2 = randn(1,1); % GDP shock
       
        if st_temp(dd)==2
            y_draw = mu_2+beta_2*f_draw(jj)+gamma_2*m_draw(jj)+sig_2*eta2; % MAKE IT DEPENDENT ON THE F AND M SHOCKS
        else
            y_draw = mu_1+beta_1*f_draw(jj)+gamma_1*m_draw(jj)+sig_1*eta2;
        end
        
%         if dd==rr % randomly assign a state for time jj
%             ind_st = round(rand(1)*(length(st_temp))+0.5,0);
%             st_draw(jj) = st_temp(ind_st);
%         end
        st_mat(dd,jj) = st_temp(dd);
        if jj>burn
            y_mat(dd,jj-burn) = y_draw;
        end
    end
    if jj>burn
        f_mat(jj-burn) = f_draw(jj); m_mat(jj-burn) = m_draw(jj);
        p12_draw(jj-burn) = p12; p21_draw(jj-burn) = p21;
    end
end

%% Quantiles

% 12-months ahead
y_mat_fut = NaN(size(y_mat,1),size(y_mat,2)-11);
for hh=1:size(y_mat,1)
    for ww=1:size(y_mat,2)-11
        y_mat_fut(hh,ww) = mean(y_mat(hh,ww:ww+11));
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

figure; clf;
hold on
l1=plot(1:tt-burn, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:tt-burn, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:tt-burn, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:tt-burn, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution 1-month ahead','FontSize',16','Interpreter','Latex');
axis tight

figure; clf;
hold on
l1=plot(1:tt-burn-11, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:tt-burn-11, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:tt-burn-11, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:tt-burn-11, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');

legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
axis tight

% States
figure; clf;
hold on
plot(1:tt-burn, mean(st_mat(:,burn+1:end)),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','States');
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([1 2]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Average states','FontSize',16','Interpreter','Latex');


figure; clf;
hold on
l1=plot(1:tt-burn, p12_draw,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Good to bad');
l2=plot(1:tt-burn, p21_draw,'Color',colors(15,:),'LineWidth', 1.5,'DisplayName','Bad to good');
legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([0 1]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Switching probabilities','FontSize',16','Interpreter','Latex');

%% SAVE
Res_iterated.y_mat = y_mat; Res_iterated.FF = f_mat; Res_iterated.MF = m_mat;
Res_iterated.dY_25 = dY_25; Res_iterated.dY_75 = dY_75;
Res_iterated.dY_10 = dY_10; Res_iterated.dY_90 = dY_90;

Res_direct.y_mat = y_mat_fut; Res_direct.FF = f_mat(1:end-12); Res_direct.MF = m_mat(1:end-12);
Res_direct.dY_25 = dY_25_fut; Res_direct.dY_75 = dY_75_fut;
Res_direct.dY_10 = dY_10_fut; Res_direct.dY_90 = dY_90_fut;


% FF = c_1_1 + a0_1_1*FF + a0_1_2*MF + a0_1_3*GDP + a1_1_1*FF(-1) + a1_1_2*MF(-1) + a1_1_3*GDP(-1)
% FF s.d. = s_1_1
% MF = c_2_1 + a0_2_1*FF + a0_2_2*MF + a0_2_3*GDP + a1_2_1*FF(-1) + a1_2_2*MF(-1) + a1_2_3*GDP(-1)
% MF s.d. = s_2_2
% GDP, normal = c_3_1_sync_1 + a0_3_1_sync_1*FF + a0_3_2_sync_1*MF + a0_3_3_sync_1*GDP
% GDP, normal s.d. = s_3_3_sync_1
% GDP, bad = c_3_1_sync_2 + a0_3_1_sync_2*FF + a0_3_2_sync_2*MF + a0_3_3_sync_2*GDP
% GDP, bad s.d. = s_3_3_sync_2
% 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
% 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'

dgp.c_1_1 = A_FF(2); dgp.a0_1_1 = 0; dgp.a0_1_2 = 0; dgp.a0_1_3 = 0; dgp.a1_1_1 = A_FF(1); dgp.a1_1_2 = 0; dgp.a1_1_3=0; dgp.s_1_1 = STD_FF;
dgp.c_2_1 = A_MF(2); dgp.a0_2_1 = 0; dgp.a0_2_2 = 0; dgp.a0_2_3 = 0; dgp.a1_2_1 = 0; dgp.a1_2_2 = A_MF(1); dgp.a1_2_3=0; dgp.s_2_2 = STD_MF;
dgp.c_3_1_sync_1 = mu_1; dgp.a0_3_1_sync_1 = beta_1; dgp.a0_3_2_sync_1 = gamma_1; dgp.a0_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = sig_1;
dgp.c_3_1_sync_2 = mu_2; dgp.a0_3_1_sync_2 = beta_2; dgp.a0_3_2_sync_2 = gamma_2; dgp.a0_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = sig_2;
dgp.a12 = a12; dgp.a21 = a21; dgp.b12 = b12; dgp.b21 = b21; dgp.c12 = c12; dgp.c21 = c21;

save Res_sim Res_iterated Res_direct dgp
