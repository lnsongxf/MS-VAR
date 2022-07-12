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
nParamDraws  = 10000; % samp.options.N;    % Number of parameter draws
nDraws       = 1500;  % Number of simulations for each parameter draw
paramsUse    = 'mcmc';


% Data vintage, sample and country selection
datafilename = '11302020';
sheetuse     = 'US_DFM';
start_date   = '1973-Jan';
end_date     = '2020-Oct';

startdb = '1973M1';
enddb   = '2019M10';



% Define dates for plots
start_plot    = '1973-Feb';
end_plot      = '2019-Oct';
end_plot_full = '2020-Oct';

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);


% Figure options
FontSize  = 16;
numticks  = 48;
figSize   = [12 6];
linestyle = {'-','--',':'};
colors    = cbrewer('div', 'RdYlBu', 64);
colors2   = cbrewer('qual', 'Set1', 8);
left_color= [0 0 0];
right_color= colors(2,:);


%========================================
% MODEL SPECIFIC VARIABLES
%========================================
% Path for stored parameters
paramsfolder_Direct   = ['results/Params/' paramsUse '/Direct/'];
paramsfolder_Iterated = ['results/Params/' paramsUse '/Iterated/'];

% Path for saving simulated quantiles
simulfolder_Direct   = ['results/Simulation/' paramsUse '/Direct/'];
simulfolder_Iterated = ['results/Simulation/' paramsUse '/Iterated/'];

% Path for saving densities
densityfolder_Direct   = ['results/Densities/' paramsUse '/Direct/'];
densityfolder_Iterated = ['results/Densities/' paramsUse '/Iterated/'];


% Path for Figures
fig_folder   = ['results/Figures/' paramsUse '/'];

%===============================
% LOAD
for hh=horizons
   modelname_dfm = [datafilename '_US_DFM_' startdb enddb ...
        '_' 'C' num2str(const) 'N' num2str(normal) 'H' num2str(hh)];
   modelname_ipebp = [datafilename '_US_IPEBP_' startdb enddb ...
        '_' 'C' num2str(const) 'N' num2str(normal) 'H' num2str(hh)];

    figureTitleTag = ['h=' num2str(hh) ', ' startdb ':' enddb];

sim1 = load([simulfolder_Direct modelname_dfm '.mat']);
sim2 = load([simulfolder_Direct modelname_ipebp '.mat']);


sim3 = load([simulfolder_Iterated modelname_dfm '.mat']);
sim4 = load([simulfolder_Iterated modelname_ipebp '.mat']);

% Plot figure
fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on;
l1=plot(dates_full, sim1.quantDirect.dYsim_10,'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','DFM: Direct 10th');
l2=plot(dates_full, sim1.quantDirect.dYsim_90,'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','DFM: Direct 90th');
l3=plot(dates_full, sim2.quantDirect.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','IPEBP: Direct 10th');
l4=plot(dates_full, sim2.quantDirect.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','IPEBP: Direct 90th');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
ylabel('Percent','interpreter','Latex','fontsize',10)
hold off
title(['Quantiles Comparison Direct Model: ' figureTitleTag],'Interpreter','Latex');
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(1), dates_full(end)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(1), dates_full(end)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig1,'PaperOrientation','portrait');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

fig1b=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on;
l1=plot(dates_full, sim3.quantIterated.dYsim_10,'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','DFM: Iterated 10th');
l2=plot(dates_full, sim3.quantIterated.dYsim_90,'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','DFM: Iterated 90th');
l3=plot(dates_full, sim4.quantIterated.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','IPEBP: Iterated 10th');
l4=plot(dates_full, sim4.quantIterated.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','IPEBP: Iterated 90th');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
ylabel('Percent','interpreter','Latex','fontsize',10)
hold off
title(['Quantiles Comparison Iterated Model: ' figureTitleTag],'Interpreter','Latex');
ax=gca;
ax.XTick = datenum(dates_full(sd:numticks:ed_full));
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(1), dates_full(end)])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'TickLabelInterpreter','Latex')
set(fig1b,'PaperOrientation','portrait');
set(fig1b, 'PaperSize', figSize);
set(fig1b, 'PaperUnits', 'inches');
set(fig1b, 'Units','inches');
set(fig1b, 'PaperPositionMode', 'auto');
set(fig1b, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;


if saveit==1
    print('-dpdf',fig1,[fig_folder 'CompareQuantilesDirect_' strrep(modelname_dfm,'US_DFM','US_DFMIPEBP')],'-bestfit');
    print('-dpdf',fig1b,[fig_folder 'CompareQuantilesIterated_' strrep(modelname_dfm,'US_DFM','US_DFMIPEBP')],'-bestfit');

end


%=========================================================
% COMPARE DENSITIES FOR SELECTED DATES
%=========================================================

den1 = load([densityfolder_Direct modelname_dfm '.mat']);
den2 = load([densityfolder_Direct modelname_ipebp '.mat']);
den3 = load([densityfolder_Iterated modelname_dfm '.mat']);
den4 = load([densityfolder_Iterated modelname_ipebp '.mat']);

target_dates = den1.target_dates;
date_index = find(ismember(dates_full, den1.target_dates));

   %==========================================================================
    % Predictive density at specified episodes
    %==========================================================================
    for ii=1:numel(date_index)
        % Direct
        [pdfi_direct1(ii,:),xi_direct1(ii,:)]  = ksdensity(den1.pdenDirect(:,ii));
        [pdfi_direct2(ii,:),xi_direct2(ii,:)]  = ksdensity(den2.pdenDirect(:,ii));
        
        % Iterated
        [pdfi_iterated1(ii,:),xi_iterated1(ii,:)]  = ksdensity(den3.pdenIterated(:,ii));
        [pdfi_iterated2(ii,:),xi_iterated2(ii,:)]  = ksdensity(den4.pdenIterated(:,ii));
    end
    
    %==========================================================================
    % Density plots at specified dates
    %==========================================================================
    
    
        fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;        
        hold on
        for ii=1:numel(date_index)
        plot(xi_direct1(ii,:),pdfi_direct1(ii,:),'Color',colors(ii+(ii-1)*15,:),'LineWidth', 3,'DisplayName',['Direct DFM: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
        plot(xi_direct2(ii,:),pdfi_direct2(ii,:),'--','Color',colors(ii+(ii-1)*15,:),'LineWidth', 3,'DisplayName',['Direct IPEBP: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
        end
        ylimits = ylim; xlimtis = xlim;
        vline(0,'k--');
        legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
        title(['Posterior Mode Predictive Density, Direct: ' figureTitleTag],'FontSize',16','Interpreter','Latex');
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(hh) '-Months-Ahead GDP Growth' ],'fontsize',10,'Interpreter','Latex')
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', FontSize);
        set(gca,'Layer','top')
        set(gca,'XTick',-18:2:8)
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        xlim([-18 8])
        set(fig2,'PaperOrientation','portrait');
        set(fig2, 'PaperSize', figSize);
        set(fig2, 'PaperUnits', 'inches');
        set(fig2, 'Units','inches');
        set(fig2, 'PaperPositionMode', 'auto');
        set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
        tightfig;
        
        fig3=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
        hold on
        for ii=1:numel(date_index)
        plot(xi_iterated1(ii,:),pdfi_iterated1(ii,:),'Color',colors(ii+(ii-1)*4,:),'LineWidth', 3,'DisplayName',['Iterated DFM: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
        plot(xi_iterated2(ii,:),pdfi_iterated2(ii,:),'--','Color',colors(ii+(ii-1)*4,:),'LineWidth', 3,'DisplayName',['Iterated IPEBP: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
        end
        ylimits = ylim; xlimtis = xlim;
        vline(0,'k--');
        legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
        title(['Posterior Mode Predictive Density, Iterated: ' figureTitleTag],'FontSize',16','Interpreter','Latex');
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(hh) '-Months-Ahead GDP Growth'],'fontsize',10,'Interpreter','Latex')
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', FontSize);
        set(gca,'Layer','top')
        set(gca,'XTick',-18:2:8)
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        xlim([-18 8])        
        set(fig3,'PaperOrientation','portrait');
        set(fig3, 'PaperSize', figSize);
        set(fig3, 'PaperUnits', 'inches');
        set(fig3, 'Units','inches');
        set(fig3, 'PaperPositionMode', 'auto');
        set(fig3, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
        tightfig;
        
        if saveit==1
            print('-dpdf',fig2,[fig_folder 'CompareDensityDirect_' datestr(target_dates(1),dataformat) strrep(modelname_dfm,'US_DFM','US_DFMIPEBP')],'-bestfit');
            print('-dpdf',fig3,[fig_folder 'CompareDensityIterated_' datestr(target_dates(1),dataformat) strrep(modelname_dfm,'US_DFM','US_DFMIPEBP')],'-bestfit');
        end
        
    


end
