%==========================================================================
% This script compares the simulated quantiles and the predicitive density
% for a specified period in the iterated model under two alternatives ways
% to obtain the current state s(t): a) drawing from filtered probabilities
% b) by simulating conditional on s(t-1). 
%
%
% To run this script you need to :
% 1. Set tperiods = ed; %length(dates_full);  % Number of time-periods in simulation of GDP
% 2. Set an individual horizon
% 3. Select a target date within sample
%==========================================================================

% Figure to check difference between drawsing st or forecasting st  in the
% iterated model

fig3=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(1:tperiods), quantIterated.dYsim_10,'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','Iterated 10th, forecast st');
l2=plot(dates_full(1:tperiods), quantIterated.dYsim_90,'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','Iterated 90th, forecast st');
l5=plot(dates_full(1:tperiods), quantIteratedDrawst.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','Iterated 10th, draw st');
l6=plot(dates_full(1:tperiods), quantIteratedDrawst.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','Iterated 90th, draw st');
legend([l1 l2 l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
ylabel('Percent','interpreter','Latex','fontsize',10)
setfig(fig3,dates_full,sd,tperiods,numticks,FontSize,figSize,['Posterior Mode Quantiles: ' figureTitleTag] );
hold off


for ii=1:numel(date_index)
    % Iterated
    [pdfi_iterated(ii,:),xi_iterated(ii,:)]  = ksdensity(pdenIterated(:,ii));
    
    % Iterated drawst
    [pdfi_iterated_drawst(ii,:),xi_iterated_drawst(ii,:)]  = ksdensity(pdenIteratedDrawst(:,ii));
end

fig4=figure;
hold on
for ii=1:numel(date_index)
    plot(xi_iterated(ii,:),pdfi_iterated(ii,:),'Color',colors(ii+(ii-1)*4,:),'LineWidth', 3,'DisplayName',['Iterated Forecast St: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
    plot(xi_iterated_drawst(ii,:),pdfi_iterated_drawst(ii,:),'--','Color',colors(ii+(ii-1)*4,:),'LineWidth', 3,'DisplayName',['Iterated Draw St: ' datestr(dates_full(date_index(ii)),dataformat)] ); hold on;
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
set(gca,'XTick',-16:2:6)
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-12 8])
set(fig4,'PaperOrientation','portrait');
set(fig4, 'PaperSize', figSize);
set(fig4, 'PaperUnits', 'inches');
set(fig4, 'Units','inches');
set(fig4, 'PaperPositionMode', 'auto');
set(fig4, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

if saveit==1
    print('-dpdf',fig3,[fig_folder 'QuantilesIteratedCompareSt_' modelname],'-bestfit');
    print('-dpdf',fig4,[fig_folder 'DensityIteratedCompareSt_' datestr(target_dates,dataformat) '_' modelname],'-bestfit');
end