%======================================================================
% COVID-19 DENSITIES
% called in msreg_2reg_endoprob_SVAR_FULL.m
%======================================================================

DATA_MAR13 = readtable('data/MacroRisk_March13_US_only.xlsx','Sheet','DFM_73_Monthly');
DATA_APR2 = readtable('data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
DATA_CURRENT = readtable('data/temp.xlsx','Sheet','US_DFM');

[~, db_mar13,~,~,~] = fLoadData('MacroRisk_March13_US_only','DFM_73_Monthly',start_date,end_date,opt.hh);
[~, db_apr2,~,~,~] = fLoadData('MacroRisk_April2_US_only','DFM_73_Monthly',start_date,end_date,opt.hh);


%% Dates for density

%     covid_dates=[
%         datenum('2020-Apr',dataformat)
%         ];
covid_dates=[
    datenum('2020-Mar',dataformat)
    ];
opt_covid = opt;
opt_covid.date_index = find(ismember(dates_full, covid_dates));
opt_covid.nDraws = 5000;
%     opt_covid.tperiods = opt_covid.date_index;
% March-13-2020
FF_COVID = db_mar13.FF.data;
MF_COVID = db_mar13.MF.data;
GDPGH_COVID = db_mar13.GDPGH.data;
if opt_covid.dir_or_it==1
    % Direct forecast
    [quantiles_mar,pdensity_covid_mar] = fPredDensityDirectFull(sv,params_in,...
        FF_COVID,MF_COVID,GDPGH_COVID,TR_full,opt_covid,simtype);
elseif opt_covid.dir_or_it==2
    % Iterated forecast
    [quantiles_mar, pdensity_covid_mar] = fPredDensityIteratedFull(sv,params_in,...
        FF_COVID,MF_COVID,GDPGH_COVID,TR_full,opt_covid,simtype);
end
%%
% April-2-2020
FF_COVID = db_apr2.FF.data;
MF_COVID = db_apr2.MF.data;
GDPGH_COVID = db_apr2.GDPGH.data;

if opt_covid.dir_or_it==1
    % Direct forecast
    [quantiles_apr,pdensity_covid_apr] = fPredDensityDirectFull(sv,params_in,...
        FF_COVID,MF_COVID,GDPGH_COVID,TR_full,opt_covid,simtype);
elseif opt_covid.dir_or_it==2
    % Iterated forecast
    [quantiles_apr, pdensity_covid_apr] = fPredDensityIteratedFull(sv,params_in,...
        FF_COVID,MF_COVID,GDPGH_COVID,TR_full,opt_covid,simtype);
end

%%
% Last Obs.


%     if opt_covid.dir_or_it==1
%         % Direct forecast
%         [quantiles_last,pdensity_covid_last] = fPredDensityDirectFull_2per(sv,params_in,...
%         FF_full,MF_full,GDPGH_full,TR_full,opt_covid,modelspec);
%     elseif opt_covid.dir_or_it==2
%         % Iterated forecast
%         [quantiles_last, pdensity_covid_last] = fPredDensityIteratedFull(sv,params_in,...
%             FF_COVID,MF_COVID,GDPGH_COVID,TR_COVID,opt_covid,'forecastst',modelspec);
%     end

%% Compute densities for MS model
% results correspond to the posterior parameter draws

[pdfi_mar,xi_mar]  = ksdensity(pdensity_covid_mar.realized);
[pdfi_mar_good,xi_mar_good]  = ksdensity(pdensity_covid_mar.good);
[pdfi_mar_bad,xi_mar_bad]  = ksdensity(pdensity_covid_mar.bad);
preg_mar= 1-quantiles_mar.st_t_mean(opt_covid.date_index);
preg_apr= 1-quantiles_apr.st_t_mean(opt_covid.date_index);

[pdfi_apr,xi_apr]  = ksdensity(pdensity_covid_apr.realized);
[pdfi_apr_good,xi_apr_good]  = ksdensity(pdensity_covid_apr.good);
[pdfi_apr_bad,xi_apr_bad]  = ksdensity(pdensity_covid_apr.bad);

% [pdfi_last,xi_last]  = ksdensity(pdensity_covid_last.realized);
% [pdfi_last_good,xi_last_good]  = ksdensity(pdensity_covid_last.good);
% [pdfi_last_bad,xi_last_bad]  = ksdensity(pdensity_covid_last.bad);

%% Figures

linetype= {'--','-','-.',':','-.'};

if opt.showfig==1

    %% PDFs: Ergodic and Regime-Specific
    
    fig20=figure;
    hold on
    % march
    plot(xi_mar,pdfi_mar,linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName','Mar 13, 2020');
    plot(xi_mar_bad,pdfi_mar_bad,linetype{1},'Color',[1 0 0 ],'LineWidth', 3,'DisplayName','Mar 13, 2020 (bad)');
    plot(xi_mar_good,pdfi_mar_good,linetype{1},'Color',[0 0 1 ],'LineWidth', 3,'DisplayName','Mar-2020 (good)');
    % april
    plot(xi_apr,pdfi_apr,linetype{2},'Color',colors(50,:),'LineWidth', 3,'DisplayName','Apr 2, 2020');
    plot(xi_apr_bad,pdfi_apr_bad,linetype{4},'Color',[1 0 0 ],'LineWidth', 3,'DisplayName','Apr 2, 2020 (bad)');
    plot(xi_apr_good,pdfi_apr_good,linetype{4},'Color',[0 0 1 ],'LineWidth', 3,'DisplayName','Apr 2, 2020 (good)');
    % plot(xi_last,pdfi_last,linetype{2},'Color',colors(20,:),'LineWidth', 3,'DisplayName','Apr-2020');
    % plot(xi_last_bad,pdfi_last_bad,linetype{4},'Color',[1 0 0 ],'LineWidth', 3,'DisplayName','Apr-2020 (bad)');
    % plot(xi_last_good,pdfi_last_good,linetype{4},'Color',[0 0 1 ],'LineWidth', 3,'DisplayName','Apr-2020 (good)');
    axis tight
    ylimits = ylim; xlimits = xlim;
    %plot(zeros(1,length(xi_mar)),zeros(1,length(pdfi_mar)),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(567))*100,1)) '\%']);
    vline(0,'k--');
    %vline(GDPGH_full_wt(567),'r:','Realized Value');
    legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
    ylabel('PDF','fontsize',10,'interpreter','Latex')
    xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    %set(gca,'XTick',unique(sort([0 round(GDPGH_full(opt.date_index(ii)),2) xlimits(1):2:xlimits(2)])))
    set(gca,'TickLabelInterpreter','Latex')
    axis tight
    xlim([-10 8])
    set(fig20,'PaperOrientation','portrait');
    set(fig20, 'PaperSize', figSize);
    set(fig20, 'PaperUnits', 'inches');
    set(fig20, 'Units','inches');
    set(fig20, 'PaperPositionMode', 'auto');
    set(fig20, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    if opt.saveit
        print('-dpdf',fig20,[figfolder 'CovidRegimePDF_' modelname ],'-bestfit');
    end
    
    %% Ergodic PDFs
    
    fig21=figure;
    hold on
    plot(xi_mar,pdfi_mar,linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName','Mar 13, 2020');
    plot(xi_apr,pdfi_apr,linetype{2},'Color',colors(50,:),'LineWidth', 3,'DisplayName','Apr 2, 2020');
    % plot(xi_last,pdfi_last,linetype{2},'Color',colors(20,:),'LineWidth', 3,'DisplayName','Apr-2020');
    axis tight
    ylimits = ylim; xlimits = xlim;
    %plot(zeros(1,length(xi_mar)),zeros(1,length(pdfi_mar)),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(567))*100,1)) '\%']);
    vline(0,'k--');
    %vline(GDPGH_full_wt(567),'r:','Realized Value');
    legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
    ylabel('PDF','fontsize',10,'interpreter','Latex')
    xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    %set(gca,'XTick',unique(sort([0 round(GDPGH_full(opt.date_index(ii)),2) xlimits(1):2:xlimits(2)])))
    set(gca,'TickLabelInterpreter','Latex')
    axis tight
    xlim([-10 8])
    set(fig21,'PaperOrientation','portrait');
    set(fig21, 'PaperSize', figSize);
    set(fig21, 'PaperUnits', 'inches');
    set(fig21, 'Units','inches');
    set(fig21, 'PaperPositionMode', 'auto');
    set(fig21, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    if opt.saveit
        print('-dpdf',fig21,[figfolder 'CovidPDF_' modelname ],'-bestfit');
    end
    
    %% Quantiles
    fig22=figure;
    l1=plot(dates_full,1-quantiles_mar.st_t_mean,'LineWidth',2,'DisplayName','March Vintage'); 
    hold on
    l2=plot(dates_full,1-quantiles_apr.st_t_mean,'--','LineWidth',2,'DisplayName','April Vintage');
    axis tight
    ylabel('Probability of Bad Regime','fontsize',10,'interpreter','Latex')
    xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
    datetick('x','yyyy','keepticks')
    %set(gca, 'XLim', [dates_full(sd), dates_full(end)])
    ylim([0 1])
    hBands=recessionplot;
    set(hBands,'FaceColor','r','FaceAlpha',0.3)
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    % Set figure appereance
    set(fig22,'PaperOrientation','portrait');
    set(fig22, 'PaperSize', figSize);
    set(fig22, 'PaperUnits', 'inches');
    set(fig22, 'Units','inches');
    set(fig22, 'PaperPositionMode', 'auto');
    set(fig22, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    if opt.saveit
        print('-dpdf',fig22,[figfolder 'CovidQuantiles_' modelname ],'-bestfit');
    end

end