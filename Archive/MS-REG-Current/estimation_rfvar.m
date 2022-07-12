%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t exogenous
% Author: Francesca Loria
% This Version: June 2020


%% housekeeping
clear
close all
% clc()

saveit=0;
save_mcmc=1;

addpath('/scratch/m1fxl02/RISE_toolbox');
addpath(genpath('scripts'));

%% Load RISE

rise_startup()

%% load the data

load data

%% plot the data

figure;
subplot(3,1,1)
plot(db.GDPGH,'linewidth',2)
title('Average Future GDP Growth')
subplot(3,1,2)
plot(db.FF,'linewidth',2)
title('Financial Factor')
subplot(3,1,3)
plot(db.MF,'linewidth',2)
title('Macroeconomic Factor')

varlist={'FF','MF','GDPGH'}';

%% set up Markov chains

chain=2;

if chain==1
    
    markov_chains=struct('name','sync',...
        'number_of_states',2,...
        'controlled_parameters',{{'c'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    model = 'Exogenous_ConstVol';

elseif chain==2
    
    markov_chains=struct('name','sync',...
        'number_of_states',2,...
        'controlled_parameters',{{'c','s'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    model = 'Exogenous';
    
elseif chain==3

    markov_chains=struct('name','vol',...
        'number_of_states',2,...
        'controlled_parameters',{{'s'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

    markov_chains(2)=struct('name','coef',...
        'number_of_states',2,...
        'controlled_parameters',{{'c'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
elseif chain==4

    markov_chains=struct('name','coef1',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,1)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

    markov_chains(2)=struct('name','coef2',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,2)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

    
    markov_chains(3)=struct('name','coef3',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,3)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

elseif chain==5

    markov_chains=struct('name','vol',...
        'number_of_states',2,...
        'controlled_parameters',{{'s'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

    markov_chains(2)=struct('name','coef1',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,1)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
    markov_chains(3)=struct('name','coef2',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,2)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
    markov_chains(4)=struct('name','coef3',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,3)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);


elseif chain==6

    markov_chains=struct('name','meanvol',...
        'number_of_states',2,...
        'controlled_parameters',{{'s','c(1,1)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
    markov_chains(2)=struct('name','coef2',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,2)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
    markov_chains(3)=struct('name','coef3',...
        'number_of_states',2,...
        'controlled_parameters',{{'c(1,3)'}},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);

end


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

%% set up VAR model


endog={'GDPGH'};
nlags=0;
exog={'FF','MF'};

const=true;

panel=[];

v=rfvar(endog,exog,nlags,const,panel,markov_chains);

%% set priors

var_prior=rfvar.prior_template();
% Sims-Zha Prior
% L1 : Overall tightness
% L2 : Cross-variable specific variance parameter
% L3 : Speed at which lags greater than 1 converge to zero
% L4 : tightness on deterministic/exogenous terms
% L5 : covariance dummies
% L6 : co-persistence
% L7 : Own-persistence
var_prior.type='sz';
var_prior.L1=1;
var_prior.L2=1;
var_prior.L3=1;
var_prior.L4=100;
var_prior.L5=1; % T-k
var_prior.L6=10; % cointegration
var_prior.L7=10; % unit root
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both
%Table A-1: Sims and Zha (1998) Normal-Wishart Prior Hyperparameters
% Parameter  Range      Chosen Value              Interpretation
% lambda1   [0, 1]          1            Overall tightness of the prior
% lambda2   [0, 1]          1            Tightness of the prior around own lags
% lambda3     = 1           1            Weight on own versus other lags
% lambda4     > 0           1            Lag decay (1 for harmonic)
% lambda5    >=0           0            Tightness of prior around the intercept
% lambda6    >=0 [0.01/0.02/0.04/infty] Sum-of-coefficients prior weight (unit root)
% lambda7    >=0         infty          Dummy-initial-observations prior weight (cointegration)

if chain==1 || chain==2
    % priors for the sync transition probabilities
    %------------------------------------------------
    switch_prior=struct();
    switch_prior.dirichlet_1={0.2,'sync_tp_1_2',0.2};
    switch_prior.dirichlet_2={0.2,'sync_tp_2_1',0.2};
elseif chain==3
    % priors for the coef transition probabilities
    %------------------------------------------------
    switch_prior=struct();
    switch_prior.coef_tp_1_2={0.5,0.1,0.3,'beta'};
    switch_prior.coef_tp_2_1={0.5,0.1,0.3,'beta'};
    % priors for the vol transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_1={0.2,'vol_tp_1_2',0.1};
    switch_prior.dirichlet_2={0.2,'vol_tp_2_1',0.1};
elseif chain==4
    % priors for the coef1 transition probabilities
    %------------------------------------------------
    switch_prior=struct();
    switch_prior.dirichlet_1={0.2,'coef1_tp_1_2',0.1};
    switch_prior.dirichlet_2={0.2,'coef1_tp_2_1',0.1};
    % priors for the coef2 transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_3={0.2,'coef2_tp_1_2',0.1};
    switch_prior.dirichlet_4={0.2,'coef2_tp_2_1',0.1};
    % priors for the coef3 transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_5={0.2,'coef3_tp_1_2',0.1};
    switch_prior.dirichlet_6={0.2,'coef3_tp_2_1',0.1};
elseif chain==5
    % priors for the coef1 transition probabilities
    %------------------------------------------------
    switch_prior=struct();
    switch_prior.dirichlet_1={0.2,'coef1_tp_1_2',0.1};
    switch_prior.dirichlet_2={0.2,'coef1_tp_2_1',0.1};
    % priors for the coef2 transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_3={0.2,'coef2_tp_1_2',0.1};
    switch_prior.dirichlet_4={0.2,'coef2_tp_2_1',0.1};
    % priors for the coef3 transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_5={0.2,'coef3_tp_1_2',0.1};
    switch_prior.dirichlet_6={0.2,'coef3_tp_2_1',0.1};
    % priors for the vol transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_7={0.2,'vol_tp_1_2',0.1};
    switch_prior.dirichlet_8={0.2,'vol_tp_2_1',0.1};
elseif chain==6
    % priors for the coef2 transition probabilities
    %------------------------------------------------
    switch_prior=struct();
    switch_prior.dirichlet_1={0.2,'coef2_tp_1_2',0.1};
    switch_prior.dirichlet_2={0.2,'coef2_tp_2_1',0.1};
    % priors for the coef3 transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_3={0.2,'coef3_tp_1_2',0.1};
    switch_prior.dirichlet_4={0.2,'coef3_tp_2_1',0.1};
    % priors for the meanvol transition probabilities
    %------------------------------------------------
    switch_prior.dirichlet_5={0.2,'meanvol_tp_1_2',0.1};
    switch_prior.dirichlet_6={0.2,'meanvol_tp_2_1',0.1};
end

prior=struct();
prior.var=var_prior;
prior.nonvar=switch_prior;


%% Find posterior mode

options = optimset('TolFun',1e-14,'MaxIter',100000);
sv = estimate(v,db,{'1973M1','2019M5'},prior,[],'fminunc');
% fminsearch, fmincon, fminunc, fminbnd


%% estimates

pmode=posterior_mode(sv)


%% Printing solution

diary ([log_folder 'reduced.txt'])
print_solution(sv)
diary off

%% plots probabilities against data

%{
%plot_probabilities(sv)
plot_data_against_probabilities(sv,'regime')
%}

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

mddobj=mdd(results,ff,lb,ub);

diary ([log_folder 'mdd.txt'])

fprintf('Importance sampling::%0.2f\n',is(mddobj,[],mdd.global_options))

fprintf('Reciprocal Importance sampling::%0.2f\n',ris(mddobj,[],mdd.global_options))

fprintf('Meng and Wong''s bridge::%0.2f\n',bridge(mddobj,true,mdd.global_options))

fprintf('Ulrich Mueller::%0.2f\n',mueller(mddobj,[],mdd.global_options))

fprintf('Laplace::%0.2f\n',laplace(mddobj))

fprintf('Sims, Waggoner and Zha::%0.2f\n',swz(mddobj,[],mdd.global_options))

fprintf('Laplace MCMC::%0.2f\n',laplace_mcmc(mddobj))

fprintf('Chib and Jeliazkov::%0.2f\n',cj(mddobj,[],mdd.global_options))

diary off

%% irfs distributions

%{
shock_names=[];

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

%%{
[Resids,Fits]=residuals(sv);

% Residuals
res_names=fieldnames(Resids);

titel='Residuals';

figure('name',titel);
    
thisname=res_names{ii};

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
%}


%% Solve

params=[results.pop.x];
myparams=params;
T = length(db.GDPGH.data)-nlags;
p_reg1 = zeros(T,options.N);
p_reg2 = zeros(T,options.N);
y_erg = zeros(T,options.N);

parfor i=1:options.N
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
    
    [Resids,Fits]=residuals(sv,myparams(:,i));
    [~,~,~,f]=filter(sv,myparams(:,i));
    
    fit1 = Fits.GDPGH.data(:,1)+db.TRENDH.data(nlags+1:end);
    fit2 = Fits.GDPGH.data(:,2)+db.TRENDH.data(nlags+1:end);
    resid1 = Resids.shock_1.data(:,1);
    resid2 = Resids.shock_1.data(:,2);
    y1 = fit1; 
    y2 = fit2; 
    
    p_reg1(:,i) = f.smoothed_regime_probabilities.regime1.data;
    p_reg2(:,i) = f.smoothed_regime_probabilities.regime2.data;    
    y_erg(:,i) = y1.*p_reg1(:,i) + y2.*p_reg2(:,i);
    
    
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

q10 = prctile(y_erg,10,2);
q50 = prctile(y_erg,50,2);
q90 = prctile(y_erg,90,2);


%% Time Series of Quantiles

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
    saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles']));
end

%% Time Series of Quantiles - Zoomed

numticks = 12;  

start_plot = '2006-Jan';
end_plot   = '2012-Dec';
sdz = find(datenum(start_plot,dataformat)==dates);
edz = find(datenum(end_plot,dataformat)==dates);

fig=figure;
hold on
plot(dates(sdz:edz), q10(sdz:edz),'Color',colors(1,:),'LineWidth', 2)
plot(dates(sdz:edz), q50(sdz:edz),'--','Color',colors(2,:),'LineWidth', 2)
plot(dates(sdz:edz), q90(sdz:edz),':','Color',colors(3,:),'LineWidth', 2)
hold off
leg = {'10th Quantile','Median','90th Quantile'};
ylabel('Percent','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sdz:numticks:edz));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sdz), dates(edz)])
rr=recessionplot;
legend(leg,'Orientation','Vertical','Location','Best','interpreter','Latex')
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
    print('-dpdf',fig,[fig_folder 'Quantiles_Zoom'],'-bestfit');
    saveas(fig,sprintf('%s.png',[fig_folder 'Quantiles_Zoom']));
end


%% Density for Specific Period

colors = [0    0.4470    0.7410;    
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;    
        0.4940    0.1840    0.5560;...
        0.4660 0.6740 0.1888];

%[~,mydate] = min(q10);
%md = datestr(dates(mydate),dataformat);

md = '2007-Oct';
mydate = find(datenum(md,dataformat)==dates);
od = '2015-May';
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
    saveas(fig,sprintf('%s.png',[fig_folder 'Histogram']));
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
    saveas(fig,sprintf('%s.png',[fig_folder 'KernelDensity']));
end


% Comparison with Quantile Regression Densities
if exist('ResMatch')==0
    load('/mcr/data_bfi/m1fxl02/3. Miscellanea/Coronavirus VIX/Text_Paper/Figures_Monthly/US/horizon_12/results___Country=US___GDP=GDP___SpecWith=ff_mf___Scenario=Baseline___Lags=0___Detrended=given___Sample=1973-Jan_to_2020-May.mat')
end

% Select Time Periods for Which to Plot PDFs
period_pdf = {md,od}; 

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
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
xmin = -6; 
xmax = 6; 
for t=1:length(period_pdf)
        % plot figure for in-sample pdfs    
    if t==1
        ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    end
    % PDF Plot
    yyaxis left
    plot(ax1,ResMatch.YY,PST(tind(t),:)','-','Color',dcolor(t,:),'LineWidth',2,'DisplayName',[char(period_pdf(t)) ' - QR']);
    hold(ax1,'on')
    yyaxis right
    [pdf,xi]=ksdensity(y_erg(tind(t),:));
    plot(ax1,xi,pdf,'-.','LineWidth', 2,'Color',dcolor(t,:),'DisplayName',[char(period_pdf(t)) ' - MS-SVAR']);
end
axis tight
yl = ylim;
realized = db.GDPGH.data + db.TRENDH.data;
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
    saveas(fig,sprintf('%s.png',[fig_folder 'DensityQR']));
end
    


%% Distribution of Regime Probabilities Against Data

numticks = 48;
left_color = [0 0 0];
right_color = colors(1,:);

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
hold on
plot(dates(sd:ed),db.FF.data(sd:ed),'k-.','LineWidth', 1.5)
plot(dates(sd:ed),db.MF.data(sd:ed),':','Color',colors(2,:),'LineWidth', 1.5)
hold off
yyaxis right
plotConfidenceBands(dates(sd:ed),[p_reg2_LB(sd:ed) p_reg2_M(sd:ed) p_reg2_UB(sd:ed)],2)
ylabel('Probability','interpreter','Latex','fontsize',10)
ax=gca;
ax.XTick = datenum(dates(sd:numticks:ed));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(ed)])
axis tight
rr=recessionplot;
legend('Financial','Macroeconomic','Location','Best','interpreter','Latex')
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
set(fig, 'Position', [figSize(1)/5 figSize(2) figSize(1) figSize(2)]);
tightfig;
if saveit==1
    print('-dpdf',fig,[fig_folder 'RegimeProbs'],'-bestfit');
    saveas(fig,sprintf('%s.png',[fig_folder 'RegimeProbs']));
end
