load([quantfolder modelname 'quantiles_iterated_alg2.mat']);
load([quantfolder modelname 'quantiles_iterated_dgp_alg2.mat']);
load([quantfolder modelname 'quantiles_direct_alg2.mat']);

fig4=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
hold on
l1=plot(dates_full(sd:ed), quantiles_iterated_alg2.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Iterated)');
l2=plot(dates_full(sd:ed), quantiles_iterated_alg2.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Iterated)');
l3=plot(dates_full(sd:ed), quantiles_direct_alg2.dYsim_10(sd:ed),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Direct)');
l4=plot(dates_full(sd:ed), quantiles_direct_alg2.dYsim_90(sd:ed),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Direct)');
l7=plot(dates_full(sd:ed), quantiles_iterated_dgp_alg2.dYsim_10(sd:ed),'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (d.g.p.)');
l8=plot(dates_full(sd:ed), quantiles_iterated_dgp_alg2.dYsim_90(sd:ed),'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (d.g.p.)');
l10=plot(dates_full(sd:ed), quantiles_iterated_dgp_alg2.GDPGH(sd:ed),'k-','LineWidth', 3,'DisplayName','d.g.p.');
legend([l1 l2 l3 l4 l7 l8 l10],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Quantiles comparison','Interpreter','Latex','FontSize',16);
axis tight
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
tightfig;