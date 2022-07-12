% Underatanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script produces posterior mode estimates and predictive
% densities at the posterior mode for any dataset.
%
% There are two estimation optionns
%  1. Direct: estimates the measurement equation at horizon hh.
%  2. Iterated: estimates the contemporanous measurement equation and
%     then iterates the predictive density through horizon hh.
%=====================================================================

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('../RISE_toolbox');         % RISE Toolbox
addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts
addpath(genpath('textools'));       % Tools to produce tables
addpath(genpath('auxtools'));       % Tools for plotting and other 
addpath(genpath('cbrewer'));        % Color mixing tools

% Date formats
inputformat   = 'yyyy-MMM';
dataformat    = 'yyyy-mmm';

%========================================
% USER INPUT
%========================================
dir_or_it    = 1;     % 1 = direct, 2 = iterated
hh           = 1;    % forecast horizon
const        = 1;     % 1 = have a constant in transition probability
normal       = 0;     % 1 = use normal distribution, 0 = gamma distribution
saveit       = 1;     % 1 = save figures
nParamDraws  = 1000;  % Parameter draws from posterior. has to be less thatn options.N;    
paramsUse    = 'mode';

% Data vintage, sample and country selection
datafilename = '11302020';
sheetuse     = 'US_DFM';
start_date   = '1973-Jan';
end_date     = '2020-Oct';

% Define target periods for density cuts
target_dates = [...%datenum('2020-Jan',dataformat),...
 %%datenum('2020-Feb',dataformat),...
 datenum('2020-Mar',dataformat),...
 datenum('2020-Apr',dataformat),...
 ];

% Define dates for plots
start_plot    = '1973-Feb';
end_plot      = '2019-Oct';
end_plot_full = '2020-Oct';

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

% MCMC options
options=struct();
options.alpha=0.234;
options.thin= 10;
options.burnin=10^3;
options.N = 2*10^4;
options.nchain=1;


%========================================
% MODEL SPECIFIC VARIABLES
%========================================

if dir_or_it ==1
    fcstype= 'direct';
    varlist={'FF','MF','GDPGH'}';
    % Folders
    paramsfolder = ['results/Params/' paramsUse '/Direct/'];
    texfolder   = ['results/Params/' paramsUse '/Direct/Tex/'];
    figfolder   = ['results/Figures/' paramsUse '/Direct/'];

elseif dir_or_it ==2
    fcstype= 'iterated';
    varlist={'FF','MF','GDPG'}';
    % Folders
    paramsfolder = ['results/Params/' paramsUse '/Iterated/'];
    texfolder   = ['results/Params/' paramsUse '/Iterated/Tex/'];
    figfolder   = ['results/Figures/' paramsUse '/Iterated/'];
end

%========================================
% CONFIGURE DATES AND PLOT OPTIONS
%========================================

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Vector of dates for plotting
dates    = dates_full(sd:ed);
tperiods = length(dates_full);  % Number of time-periods in simulation of GDP

% Figure options
FontSize  = 16;
numticks  = 48;
figSize   = [12 6];
linestyle = {'-','--',':'};
colors    = cbrewer('div', 'RdYlBu', 64);
colors2   = cbrewer('qual', 'Set1', 8);
left_color= [0 0 0];
right_color= colors(2,:);

% Create folders to store results
if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end
if exist(texfolder,'dir')==0;  mkdir(texfolder); end
if exist(figfolder,'dir')==0;  mkdir(figfolder); end


%==========================================================================
% LOAD DATA
%==========================================================================
% Load RISE
rise_startup()

[db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,hh);

% Collect MF, FF and trend (estimation sample)
FF    = db.FF.data(nlags+1:end);
MF    = db.MF.data(nlags+1:end);
TR    = db.TRENDH.data(nlags+1:end);
GDPG  = db.GDPG.data(nlags+1:end);
GDPGH = db.GDPGH.data(nlags+1:end);

% Collect MF, FF and trend (full sample)
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);
TR_full = db_full.TRENDH.data(nlags+1:end);

% Model name and title for posterior mode table
% Note: iterated model posterior mode depends on hh through the estimation sample

modelname = [datafilename '_' sheetuse '_' startdb enddb ...
    '_' 'C' num2str(const) 'N' num2str(normal) 'H' num2str(hh)];


%==========================================================================
% RISE ESTIMATION
%==========================================================================

if strcmp(paramsUse,'mode')
    
    % Load results or estimate posterior mode
    if exist([paramsfolder modelname '.mat'],'file')==0
        % Estimate posterior mode
        run scriptEstimateMode.m;
        
        % Store solution structure
        save([paramsfolder modelname '.mat'],'sPMode')
    else
        % Load stored results
        load([paramsfolder modelname '.mat']);
        
        % Map objects:
        sv = sPMode.sv;
    end
    
else
    
    % Load or estimate mcmc
    if exist([paramsfolder modelname '.mat'],'file')==0
        % Estimat posterior distribution
        run scriptEstimateMCMC.m

        % Store solution structure
        save([paramsfolder modelname '.mat'],'sPmcmc')

    else
        
        % Load stored results
        load([paramsfolder modelname '.mat']);
        
        % Map objects:
        sv = sPmcmc.sv;                
        pmcmc = sPmcmc.params;
    end
    
end

%==========================================================================
% Posterior Mode Analysis
%==========================================================================

% Posterior mode values
pmode=posterior_mode(sv);

% Smoothed probabilities
[~,~,~,f]=filter(sv);
p_reg1 = f.smoothed_regime_probabilities.regime1.data;
p_reg2 = f.smoothed_regime_probabilities.regime2.data;

% Filtered probabilities
p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;

% Adjusted filtered states
st_filtered = p_reg2_filtered(nlags+1:end);

% Map posterior mode coefficient estimates
fnames = fieldnames(pmode);
for jj=1:length(fnames)
    eval([fnames{jj} ' = pmode.' fnames{jj} ';']);
end
% Print Structural Form
print_structural_form(sv)

% Print Reduced Form
print_solution(sv)

return
%==========================================================================
% Predictive Density at Posterior Mode or MCMC
%==========================================================================

if strcmp(paramsUse,'mode')
    params_in   = pmode;
    nParamDraws = length(pmode.a12);
    nDraws      = 1000;  % Number of simulations for eacht time period
elseif strcmp(paramsUse,'mcmc')
    params_in   = pmcmc;    
    nDraws      = 10; 
end

% Extract the index for the desired density plots
date_index = find(ismember(dates_full, target_dates));

if dir_or_it==1
    % Direct forecast
    [quantiles,pdensity ] = fPredDensityDirect(TR_full,st_filtered,params_in,FF_full,MF_full,nParamDraws,nDraws,tperiods,hh,date_index);
elseif dir_or_it==2
    % Iterated forecast
    [quantiles, pdensity] = fPredDensityIterated(TR_full,st_filtered,params_in,FF_full,MF_full,nParamDraws,nDraws,tperiods,hh,date_index,'forecastst');
end

% % Predictive density at specified episodes
for ii=1:numel(date_index)
     [pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity(:,ii));
end
%%
%==========================================================================
% Quantile Plots 
%==========================================================================
figureTitleTag = [paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(hh) ', ' startdb ':' enddb];

% Add additional objects to structure
fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
subplot(211)
hold on
l5=plot(dates_full(sd:ed_full), quantiles.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th');
l6=plot(dates_full(sd:ed_full), quantiles.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th');
legend([l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
if dir_or_it==1
    title(['Quantiles of GDP Direct  Model: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
else
    title(['Quantiles of GDP Iterated Model: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
end
axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
%
subplot(212)
hold on
l3=plot(dates_full(sd:ed_full), 1-quantiles.st_t_mean,'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Forecast s(t)');
l4=plot(dates_full(sd:ed), st_filtered(sd:ed),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Draws s(t)');
legend([l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
%set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
ylim([0 1])
tightfig;



%==========================================================================
% Density plots
%=======================================================;===================
linetype= {'-','--','-.',':','.'};

fig2=figure;
hold on
for ii=1:numel(date_index)
plot(xi(ii,:),pdfi(ii,:),linetype{ii},'Color',colors(ii+(ii-1)*4,:),'LineWidth', 3,'DisplayName',datestr(dates_full(date_index(ii)),dataformat)); hold on;
end
if dir_or_it==1
    title(['Densities  Direct Model: ' figureTitleTag ],'Interpreter','Latex');
else
    title(['Densities Iterated Model: ' figureTitleTag ],'Interpreter','Latex');
end
ylimits = ylim; xlimtis = xlim;
vline(0,'k--');
legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel([num2str(hh) '-Months-Ahead GDP Growth'],'fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'XTick',-18:2:8)
set(gca,'TickLabelInterpreter','Latex')
axis tight
%xlim([-18 8])
set(fig2,'PaperOrientation','portrait');
set(fig2, 'PaperSize', figSize);
set(fig2, 'PaperUnits', 'inches');
set(fig2, 'Units','inches');
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

if saveit
print('-dpdf',fig,[figfolder 'Quantiles_' modelname],'-bestfit');
print('-dpdf',fig2,[figfolder 'DensitiesCOVID_' modelname],'-bestfit');
end
