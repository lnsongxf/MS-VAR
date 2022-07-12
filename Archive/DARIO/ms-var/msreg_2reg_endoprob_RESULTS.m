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
%===================================================================

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

horizons     = [3 12];     % forecast horizon
const        = 1;     % 1 = have a constant in transition probability
normal       = 0;     % 1 = use normal distribution, 0 = gamma distribution
saveit       = 1;     % 1 = save figures
num_rep      = 5000;  % Number of repetitions for the coverage bands of the DGP
nParamDraws  = 1000; % Number of parameter draws. Has to be less thatn sPmcmc.options.N; 
nDraws       = 10;    % Number of simulations for each parameter draw
paramsUse    = 'mcmc';% Select between mode and mcmc


% Data vintage, sample and country selection
datafilename = '11302020';
sheetuse     = 'US_DFM';
start_date   = '1973-Jan';
end_date     = '2020-Oct';

% Define target periods for density cuts
% yyyy-MMM
target_dates = [...
    datenum('2020-Mar',dataformat)...
    datenum('2020-Apr',dataformat)...
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


%========================================
% MODEL SPECIFIC VARIABLES
%========================================
% Path for stored parameters
paramsfolder_Direct   = ['results/Params/' paramsUse '/Direct/'];
paramsfolder_Iterated = ['results/Params/' paramsUse '/Iterated/'];

% Path for saving simulated quantiles
simulfolder_Direct   = ['results/Simulation/' paramsUse '/Direct/'];
simulfolder_Iterated = ['results/Simulation/' paramsUse '/Iterated/'];

% Path for stored densities
densityfolder_Direct   = ['results/Densities/' paramsUse '/Direct/'];
densityfolder_Iterated = ['results/Densities/' paramsUse '/Iterated/'];

% Path for Figures
fig_folder   = ['results/Figures/' paramsUse '/'];


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
if exist(fig_folder,'dir')==0; mkdir(fig_folder);  end

% Load RISE
rise_startup()

%==========================================================================
% LOOP OVER HORIZONS
%==========================================================================

for hh = horizons
    
    %==========================================================================
    % LOAD DATA
    %==========================================================================
    
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
    
    figureTitleTag = [strrep(sheetuse,'_',', ')  ', h=' num2str(hh) ', ' startdb ':' enddb];
    
    %==========================================================================
    % LOAD PREDICTIVE DENSITIES
    %==========================================================================
    % Load results or estimate posterior mode
    
    % Extract the index for the desired density plots
    date_index = find(ismember(dates_full, target_dates));
    
    % 1) Direct
    
    if exist([paramsfolder_Direct modelname '.mat'],'file')==0
        error(['First estimate Direct Predictive Density for: ' modelname] )
    else

        % Load stored results
        load([paramsfolder_Direct modelname '.mat']);

        if strcmp(paramsUse,'mode')
            % Map objects
            sv = sPMode;
            params_in=posterior_mode(sv);

        else
            % Map objects:
            sv = sPmcmc.sv;                
            params_in = sPmcmc.params;
            % Overwrite nDraws in mcmc
            nDraws    = 10;
        end


        % Filtered states at posterior mode
        [~,~,~,f]=filter(sv);
        st_filtered = f.filtered_regime_probabilities.regime2.data(nlags+1:end);
        
        % Compute quantiles and predictive density
        [quantDirect,pdenDirect] = fPredDensityDirect(TR_full,st_filtered,params_in,FF_full,MF_full,nParamDraws,nDraws,tperiods,hh,date_index);
        
        % Save simulated quantiles
        save([simulfolder_Direct modelname '.mat'],'quantDirect');
        save([densityfolder_Direct modelname '.mat'],'pdenDirect','target_dates');
        
    end
    
    % 2) Iterated
    
    if exist([paramsfolder_Iterated modelname '.mat'],'file')==0
        error(['First estimate Iterated Predictive Density for: ' modelname] )
    else
        
        % Load estimation results
        load([paramsfolder_Iterated modelname '.mat']);

        if strcmp(paramsUse,'mode')
            % Map objects
            sv = sPMode;
            params_in=posterior_mode(sv);

        else
            % Map objects:
            sv = sPmcmc.sv;                
            params_in = sPmcmc.params;
            % Overwrite nDraws in mcmc
            nDraws    = 10;

        end

        % Filtered states at posterior mode        
        [~,~,~,f]=filter(sv);
        st_filtered = f.filtered_regime_probabilities.regime2.data(nlags+1:end);
        
        % Compute quantiles and densities
        [quantIterated,pdenIterated] = fPredDensityIterated(TR_full,st_filtered,params_in,FF_full,MF_full,nParamDraws,nDraws,tperiods,hh,date_index,'forecastst');
        
        % Save simulated quantiles
        save([simulfolder_Iterated modelname '.mat'],'quantIterated');
        save([densityfolder_Iterated modelname '.mat'],'pdenIterated','target_dates');

        % Simulation drawing st
        % *** Uncomment this line only to compare drawing vs forecating st ***
        %[quantIteratedDrawst,pdenIteratedDrawst] = fPredDensityIterated(TR_full,st_filtered,params_in,FF_full,MF_full,nParamDraws,nDraws,tperiods,hh,date_index,'drawst');
    end
    
    
    %==========================================================================
    % Predictive density at specified episodes
    %==========================================================================
    for ii=1:numel(date_index)
        % Direct
        [pdfi_direct(ii,:),xi_direct(ii,:)]  = ksdensity(pdenDirect(:,ii));
        
        % Iterated
        [pdfi_iterated(ii,:),xi_iterated(ii,:)]  = ksdensity(pdenIterated(:,ii));
    end
    %%
    %==========================================================================
    % Density plots at specified dates
    %==========================================================================
    linetype= {'-','--','-.',':','.'};

        fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
        hold on
        for ii=1:numel(date_index)
            plot(xi_direct(ii,:),pdfi_direct(ii,:),'Color',colors(ii+(ii-1)*15,:),'LineWidth', 3,'DisplayName',['Direct: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
            plot(xi_iterated(ii,:),pdfi_iterated(ii,:),'--','Color',colors(ii+(ii-1)*15,:),'LineWidth', 3,'DisplayName',['Iterated: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
        end
        ylimits = ylim; xlimtis = xlim;
        vline(0,'k--');
        legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
        title(['Posterior Mode Predictive Density: ' figureTitleTag],'FontSize',16','Interpreter','Latex');
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(hh) '-Months-Ahead GDP Growth'],'fontsize',10,'Interpreter','Latex')
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', FontSize);
        set(gca,'Layer','top')
        set(gca,'XTick',-18:2:8)
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        xlim([-18 8])
        set(fig1,'PaperOrientation','portrait');
        set(fig1, 'PaperSize', figSize);
        set(fig1, 'PaperUnits', 'inches');
        set(fig1, 'Units','inches');
        set(fig1, 'PaperPositionMode', 'auto');
        set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
        tightfig;
        
        if saveit==1
            print('-dpdf',fig1,[fig_folder 'Density_' datestr(target_dates(1),dataformat) '_' modelname],'-bestfit');
        end
        %%
    
    %==========================================================================
    % Compare Time Series PLot of Quantiles
    %==========================================================================
    
    % Plot figure
    fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
    hold on
    l1=plot(dates_full, quantDirect.dYsim_10,'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','Direct 10th');
    l2=plot(dates_full, quantDirect.dYsim_90,'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','Direct 90th');
    l5=plot(dates_full, quantIterated.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','Iterated 10th');
    l6=plot(dates_full, quantIterated.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','Iterated 90th');
    legend([l1 l2 l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    ylabel('Percent','interpreter','Latex','fontsize',10)
    setfig(fig2,dates_full,sd,ed_full,numticks,FontSize,figSize,['Posterior Mode Quantiles: ' figureTitleTag] );
    hold off
    
    if saveit==1
        print('-dpdf',fig2,[fig_folder 'Quantiles_' modelname],'-bestfit');
    end

    %==========================================================================
    % Compare Simulated Probability of Bad Regime
    %==========================================================================

    
    fig3=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
    hold on
    l1=plot(dates_full, 1-quantDirect.st_t_mean,'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct');
    l2=plot(dates_full, 1-quantIterated.st_t_mean,'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Iterated');
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    %ylabel('Percent','interpreter','Latex','fontsize',10)
    setfig(fig3,dates_full,sd,ed_full,numticks,FontSize,figSize,['Simulated Probability of Bad Regime: ' figureTitleTag] );
    hold off
    
    if saveit==1
        print('-dpdf',fig3,[fig_folder 'St_' modelname],'-bestfit');
    end

end
