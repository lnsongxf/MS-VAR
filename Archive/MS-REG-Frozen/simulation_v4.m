%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all; clc;

rng(777); % Fixed seed for replication purposes

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
saveit    = 1; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution


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


%% SIM


% 2 = BAD, 1 = NORMAL
burn = 500; % 500 burn, 500 keep
tt=1200+burn; % time periods
rr = 600; % draws
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
% Make y_mat of GDP(t+1):
% y_mat = y_mat(:,2:end);
% f_mat = f_mat(:,1:end-1);
% m_mat = m_mat(:,1:end-1);
% st_mat = st_mat(:,1:end-1);
% p12_draw = p12_draw(1:end-1,:);
% p21_draw = p21_draw(1:end-1,:);


%% Quantiles

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
title('Simulated GDP distribution in t','FontSize',16','Interpreter','Latex');
axis tight

figure; clf;
hold on
l1=plot(1:tt-burn-12, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:tt-burn-12, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:tt-burn-12, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:tt-burn-12, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');

legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
axis tight

% States
figure; clf;
hold on
plot(1:tt-burn, mean(st_mat),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','States');
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([1 2]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Average states','FontSize',16','Interpreter','Latex');


figure; clf;
hold on
l1=plot(1:tt-burn, mean(p12_draw),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Good to bad');
l2=plot(1:tt-burn, mean(p21_draw),'Color',colors(15,:),'LineWidth', 1.5,'DisplayName','Bad to good');
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

% dgp.c_1_1 = A_FF(2); dgp.a0_1_1 = 0; dgp.a0_1_2 = 0; dgp.a0_1_3 = 0; dgp.a1_1_1 = A_FF(1); dgp.a1_1_2 = 0; dgp.a1_1_3=0; dgp.s_1_1 = STD_FF;
% dgp.c_2_1 = A_MF(2); dgp.a0_2_1 = 0; dgp.a0_2_2 = 0; dgp.a0_2_3 = 0; dgp.a1_2_1 = 0; dgp.a1_2_2 = A_MF(1); dgp.a1_2_3=0; dgp.s_2_2 = STD_MF;
% dgp.c_3_1_sync_1 = mu_1; dgp.a0_3_1_sync_1 = beta_1; dgp.a0_3_2_sync_1 = gamma_1; dgp.a0_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = sig_1;
% dgp.c_3_1_sync_2 = mu_2; dgp.a0_3_1_sync_2 = beta_2; dgp.a0_3_2_sync_2 = gamma_2; dgp.a0_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = sig_2;
% dgp.a12 = a12; dgp.a21 = a21; dgp.b12 = b12; dgp.b21 = b21; dgp.c12 = c12; dgp.c21 = c21;

% save Res_sim Res_iterated Res_direct dgp
