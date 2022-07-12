% Figure options
figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
FontSize  = 16;
numticks  = 48;
figSize   = [12 6];
linestyle = {'-','--',':'};
colors    = cbrewer('div', 'RdYlBu', 64);
colors2   = cbrewer('qual', 'Set1', 8);
left_color= [0 0 0];
right_color= colors(2,:);
%%
fig1=figure;clf;
subplot(2,1,1)
l1=plot(dates(sd:end),benchmark.dYsim_10(sd:end),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
l2=plot(dates(sd:end),benchmark.dYsim_90(sd:end),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
l3=plot(dates(sd:end),cf{1}.dYsim_10(sd:end),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 1');
l4=plot(dates(sd:end),cf{1}.dYsim_90(sd:end),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 1');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
if opt.dir_or_it==1
    title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
else
    title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
end
axis tight
numticks = 36;
set(gca,'XTick',datenum(dates(sd:numticks:(end))))
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates(end)])

subplot(2,1,2)
l1=plot(dates(sd:end),benchmark.dYsim_10(sd:end),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
l2=plot(dates(sd:end),benchmark.dYsim_90(sd:end),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
l3=plot(dates(sd:end),cf{2}.dYsim_10(sd:end),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 2');
l4=plot(dates(sd:end),cf{2}.dYsim_90(sd:end),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 2');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
if opt.dir_or_it==1
    title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
else
    title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
end
axis tight
numticks = 36;
set(gca,'XTick',datenum(dates(sd:numticks:(end))))
datetick('x','yyyy','keepticks')
set(gca, 'XLim', [dates(sd), dates((end))])

% Set figure appereance
set(fig1,'PaperOrientation','portrait');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;
%%
print('-dpdf',fig1,[figfolder 'QuantilesCF_' modelname ],'-bestfit');
%%
rtable.dates = datestr(dates(sd:end));
rtable.dYsim_10 = benchmark.dYsim_10;
rtable.dYsim_25 = benchmark.dYsim_25;
rtable.dYsim_75 = benchmark.dYsim_75;
rtable.dYsim_90 = benchmark.dYsim_90;    
rtable.dYsim_10_cf1 = cf{1}.dYsim_10;
rtable.dYsim_25_cf1 = cf{1}.dYsim_25;
rtable.dYsim_75_cf1 = cf{1}.dYsim_75;
rtable.dYsim_90_cf1 = cf{1}.dYsim_90;
rtable.dYsim_10_cf2 = cf{2}.dYsim_10;
rtable.dYsim_25_cf2 = cf{2}.dYsim_25;
rtable.dYsim_75_cf2 = cf{2}.dYsim_75;
rtable.dYsim_90_cf2 = cf{2}.dYsim_90;
rtable.data     = GDPGH_full_wt(sd_full:end);
tabout = struct2table(rtable);
writetable(tabout,'R/dataquantilesDirect.csv')
%% FIGURE WITH DENSITY PLOTS
count = 1;
fig2 = figure;clf; 
for ix = [92 99 242 244]
    subplot(2,2,count)
    l1=plot(benchmark.xi{ix},benchmark.pdfi{ix},'k','LineWidth',2,'DisplayName','Estimated'); hold on;
    l2=plot(cf{1}.xi{ix},cf{1}.pdfi{ix},'--','LineWidth',2,'Color',colors(10*1,:),'DisplayName',['No Macro Factor']);
    l3=plot(cf{2}.xi{ix},cf{2}.pdfi{ix},'--','LineWidth',2,'Color',colors(10*2,:),'DisplayName',['No Financial Factor']);
    title(datestr(datenum(datetime(efirstok,'InputFormat',inputformat)+calmonths(ix)),dataformat),'Interpreter','Latex')
    legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children'))); xlim([-12 6]);
    count = count + 1;
end
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     if opt.dir_or_it==1        
%         text( 0.5, 0.1, ['OOS Direct: ' figureTitleTag ], 'FontSize', 14', 'FontWeight', 'Bold', ...
%             'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
%     else
%         text( 0.5, 0.2, ['OOS Iterated: ' figureTitleTag ], 'FontSize', 14', 'FontWeight', 'Bold', ...
%             'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','Latex' ) ;        
%     end
% Set figure appereance
set(fig2,'PaperOrientation','portrait');
set(fig2, 'PaperSize', figSize);
set(fig2, 'PaperUnits', 'inches');
set(fig2, 'Units','inches');
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

print('-dpdf',fig2,[figfolder 'DensityCF_' modelname ],'-bestfit');


%% SAVE RESULTS
save([counterfolder modelname '.mat'],'benchmark','cf');


%%
return
