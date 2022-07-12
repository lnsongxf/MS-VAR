clear; clc;
% addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         % RISE Toolbox
addpath('/if/gms_research_gar/GAR/MS-VAR/MS-REG-Current/auxtools/');         % RISE Toolbox

% rise_startup()


% Default Settings
set(0,'defaulttextinterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

savefigures = 1;

% LOAD DATA
load inputirf_new.mat
load params_irf.mat
load filter_results.mat

% Load Estimated Parameter Structure
% load([modelname '.mat']);
           
% sv = sPMode.sv;
% [~,~,~,f]=filter(sv);

% Filtered probabilities
% p_reg2_updated  = f.updated_regime_probabilities.regime2.data;
 
s2data = ones(1,length(p_reg2_updated));
s2data(p_reg2_updated>0.9)=2;

% USER OPTIONS
hh = opt.hh;
nlags = 1;     %[Fixed at 1]

% number of draws of regime shocks
opt.nDrawsU  = 2000; 
opt.nDrawsEPS= 1000; 

%% IRFs
%Baseline shocks shock matrices
shocks_baseline = randn(3,hh,opt.nDrawsEPS);
if opt.nDrawsEPS ==1; shocks_baseline = shocks_baseline.*0; end

%Generate cell array of target dates and convert full date array to cell
%array
% target_dates_cell = {'01-Jan-1999','01-Jan-2000','01-Jan-2007','01-Jan-2008','01-Aug-2008','01-Jan-2009'};

target_dates_cell = {'01-Jan-1999','01-Jan-2008'};
% target_dates_cell = {'01-Jan-1999'};
% target_dates_cell = {'01-Jan-2008'};

dates_full_cell = cellstr(datestr(dates_full));

%Initialize init structure
init_irf.s0 = [];
sd_base = 0; % keep at zero

% Impulses
sd_shock.f = +3;
sd_shock.m = -3;
sd_shock.y = -3;
opt.scenario=0;

[irf, pctl] = fAllDatesAlgorithm9_exante(dates_full_cell,target_dates_cell,FF_full,MF_full,GDPG_full,s2data,p_reg2_updated,init_irf,shocks_baseline,param,opt,sd_base,sd_shock);
%% PDFs
close all

cut_tails = 0.005; % Cut of the tails
band_pdf = [0.1 0.1]; % Bandwidth adjustment (fin. and macro)
sum_cut = [cut_tails (1-cut_tails)];

for d = 1:length(dates_full_cell)
    if any(strcmp(dates_full_cell{d},target_dates_cell)) %check if the current date is one of the target dates
        
        % Get current date
        current_date = datetime(dates_full_cell{d}, 'InputFormat', 'dd-MMM-yyyy');
        
        % Transform current date to YYYYMM format
        yearmonth = strcat('Y', num2str(year(datetime(current_date))), 'M', num2str(month(datetime(current_date))));
        
        irf_temp_f = irf.(yearmonth).f_shock.y; irf_temp_m = irf.(yearmonth).m_shock.y;
        
        irf_dens_f =mean(irf_temp_f,2);
        irf_dens_f_12 =irf_temp_f(:,end);
        [pdfi_f,xi_f]  = ksdensity(irf_dens_f,'Bandwidth',band_pdf(1));
        [pdfi_f_12,xi_f_12]  = ksdensity(irf_dens_f_12,'Bandwidth',band_pdf(1));
        sum_pdfi_f = cumsum(pdfi_f);
        sum_pdfi_f = sum_pdfi_f/sum(pdfi_f);
        sum_pdfi_f_12 = cumsum(pdfi_f_12);
        sum_pdfi_f_12 = sum_pdfi_f_12/sum(pdfi_f_12);
        
        xlim_inf_temp = find(sum_pdfi_f>sum_cut(1),1,'first');
        xlim_sup_temp = find(sum_pdfi_f>sum_cut(2),1,'first');
        xlim_inf_temp_12 = find(sum_pdfi_f_12>sum_cut(1),1,'first');
        xlim_sup_temp_12 = find(sum_pdfi_f_12>sum_cut(2),1,'first');
        
        xlim_inf = min([xi_f(xlim_inf_temp) xi_f_12(xlim_inf_temp_12)]);
        xlim_sup = max([xi_f(xlim_sup_temp) xi_f_12(xlim_sup_temp_12)]);
        
        fig2=figure('Name',['PDF_financial_U_' num2str(opt.nDrawsU) '_EPS_' num2str(opt.nDrawsEPS) '_' yearmonth '']);
        hold on
        plot(xi_f,pdfi_f,'Color',[0 0 0 ],'LineWidth', 3,'DisplayName','Fin., avg. 12-month');
        plot(xi_f_12,pdfi_f_12,'Color',[0.5 0.5 0.5],'LineWidth', 3,'DisplayName','Fin., 12-month-ahead');
        axis tight
        vline(0,'k--');
        legend('Orientation','Vertical','Location','best','interpreter','Latex');legend boxoff;
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', 12);
        set(gca,'Layer','top')
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        xlim([xlim_inf-0.5 xlim_sup+0.5])
        set(fig2,'PaperOrientation','portrait');
        set(fig2, 'PaperUnits', 'inches');
        set(fig2, 'Units','inches');
        set(fig2, 'PaperPositionMode', 'auto');
        tightfig;
        
        irf_dens_m =mean(irf_temp_m,2);
        irf_dens_m_12 =irf_temp_m(:,end);
        [pdfi_m,xi_m]  = ksdensity(irf_dens_m,'Bandwidth',band_pdf(1));
        [pdfi_m_12,xi_m_12]  = ksdensity(irf_dens_m_12,'Bandwidth',band_pdf(1));
        
        sum_pdfi_m = cumsum(pdfi_m);
        sum_pdfi_m = sum_pdfi_m/sum(pdfi_m);
        sum_pdfi_m_12 = cumsum(pdfi_m_12);
        sum_pdfi_m_12 = sum_pdfi_m_12/sum(pdfi_m_12);
        
        xlim_inf_temp = find(sum_pdfi_m>sum_cut(1),1,'first');
        xlim_sup_temp = find(sum_pdfi_m>sum_cut(2),1,'first');
        xlim_inf_temp_12 = find(sum_pdfi_m_12>sum_cut(1),1,'first');
        xlim_sup_temp_12 = find(sum_pdfi_m_12>sum_cut(2),1,'first');
        
        xlim_inf = min([xi_m(xlim_inf_temp) xi_m_12(xlim_inf_temp_12)]);
        xlim_sup = max([xi_m(xlim_sup_temp) xi_m_12(xlim_sup_temp_12)]);
        
        fig3=figure('Name',['PDF_macro_U_' num2str(opt.nDrawsU) '_EPS_' num2str(opt.nDrawsEPS) '_' yearmonth '']);
        hold on
        plot(xi_m,pdfi_m,'Color',[0 0 0 ],'LineWidth', 3,'DisplayName','Macro, avg. 12-month');
        plot(xi_m_12,pdfi_m_12,'Color',[0.5 0.5 0.5],'LineWidth', 3,'DisplayName','Macro, 12-month-ahead');
        axis tight
        vline(0,'k--');
        legend('Orientation','Vertical','Location','best','interpreter','Latex');legend boxoff;
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
        
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', 12);
        set(gca,'Layer','top')
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        xlim([xlim_inf-0.5 xlim_sup+0.5])
        set(fig3,'PaperOrientation','portrait');
        set(fig3, 'PaperUnits', 'inches');
        set(fig3, 'Units','inches');
        set(fig3, 'PaperPositionMode', 'auto');
        tightfig;
    else
        continue
    end
end




%%
fields = fieldnames(pctl);
for i=1:length(fields)
    plot_fan_gdp_alg9exante(pctl.(fields{i}),hh,strcat('irfAlg9Exante_U_',num2str(opt.nDrawsU),'_EPS_',num2str(opt.nDrawsEPS),'_',fields{i},''),strcat('A9: Ex-Ante IRF, sd=',num2str(sd_shock.f),', f0=',num2str(pctl.(fields{i}).f0),', m0=',num2str(pctl.(fields{i}).m0),', ',fields{i}));
end

if savefigures ==1
    h = get(0,'children');
    h = sort(h);
    for i=1:length(h)
        set(h(i),'Units','Inches');
        pos = get(h(i),'Position');
        set(h(i),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h(i),['results/' get(h(i),'Name')],'-dpdf','-r0')
%         print(h(i),get(h(i),'Name'),'-depsc','-r0')
    end 
end
return

%% Scenario

% Impulses
sd_shock.f = +2;
sd_shock.m = -2;
sd_shock.y = 0;
opt.scenario=1;

[irf, pctl] = fAllDatesAlgorithm9_exante(dates_full_cell,target_dates_cell,FF_full,MF_full,GDPG_full,s2data,p_reg2_updated,init_irf,shocks_baseline,param,opt,sd_base,sd_shock);
fields = fieldnames(pctl);
for i=1:length(fields)
    plot_fan_gdp_alg9exante(pctl.(fields{i}),hh,strcat('Algorithm9/scenarioAlg9Exante_',fields{i},'.pdf'),strcat('A9: Ex-Ante Scenario, sd=',num2str([sd_shock.f sd_shock.m sd_shock.y]),', f0=',num2str(pctl.(fields{i}).f0),', m0=',num2str(pctl.(fields{i}).m0),', ',fields{i}));
end