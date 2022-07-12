function plot_fan_gdp_alg9exante(pctl_in,periods,filename,supt,pctl2_in)
time = linspace(1,periods,periods)';
patchx = [1:periods periods:-1:1]';
lb99 = [0.7137, 0.8667, 0.8902];
lb95 = [0.4275, 0.6706, 0.7608];
lb90 = [0.3882, 0.5647, 0.7216];
db = [38, 49, 109]/256;
lr1 = [0.8314, 0.6549, 0.6549];
lr5 = [0.8000, 0.5922, 0.5922];
lr10 = [0.6392, 0.3725, 0.3725];
dr = [174,0,0]/256;
dg = [0 0.5 0];

myalpha25 = 0.75;
myalpha10 = 0.65;
myalpha5 = 0.55;
myalpha1 = 0.4;

cmedian = [dr; dg];

% Set y-axis bounds
y_ub = 0.5; y_lb = -2;
f_ub = 0.5; f_lb = -2;
m_ub = 0.5; m_lb = -2;


% Y lim for GDP response
min_y = min([min(min(pctl_in.f_shock.y)); min(min(pctl_in.m_shock.y)); min(min(pctl_in.y_shock.y))]);
max_y = max([max(max(pctl_in.f_shock.y)); max(max(pctl_in.m_shock.y)); max(max(pctl_in.y_shock.y))]);

if min_y < y_lb;  y_lb = floor(min_y); end
if max_y > y_ub;  y_ub = ceil(max_y); end


% Y lim for FF response
min_f = min([min(min(pctl_in.f_shock.f)); min(min(pctl_in.m_shock.f)); min(min(pctl_in.y_shock.f))]);
max_f = max([max(max(pctl_in.f_shock.f)); max(max(pctl_in.m_shock.f)); max(max(pctl_in.y_shock.f))]);

if min_f < f_lb;  f_lb = floor(min_f); end
if max_f > f_ub;  f_ub = ceil(max_f); end


% Y lim for MM response
min_m = min([min(min(pctl_in.f_shock.m)); min(min(pctl_in.m_shock.m)); min(min(pctl_in.y_shock.m))]);
max_m = max([max(max(pctl_in.f_shock.m)); max(max(pctl_in.m_shock.m)); max(max(pctl_in.y_shock.m))]);

if min_m < m_lb;  m_lb = floor(min_m); end
if max_m > m_ub;  m_ub = ceil(max_m); end


% Label names for fixed regime
fnames2 = fields(pctl2_in);
irf_fixed_good = pctl2_in.(fnames2{1});
irf_fixed_bad = pctl2_in.(fnames2{2});

figure('Position',[0, 0, 1200, 800]);
subplot(2,3,1); hold on;
patch(patchx,[pctl_in.f_shock.y(1,:)'; flipud(pctl_in.f_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.f_shock.y(2,:)'; flipud(pctl_in.f_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.f_shock.y(3,:)'; flipud(pctl_in.f_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.f_shock.y(4,:)'; flipud(pctl_in.f_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.f_shock.y(5,:)'; flipud(pctl_in.f_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.f_shock.y(6,:)'; flipud(pctl_in.f_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.f_shock.y(7,:)'; flipud(pctl_in.f_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.f_shock.y(8,:)'; flipud(pctl_in.f_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.f_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.f_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Financial Shock');
title('GDP response')

subplot(2,3,2); hold on;
patch(patchx,[pctl_in.m_shock.y(1,:)'; flipud(pctl_in.m_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.m_shock.y(2,:)'; flipud(pctl_in.m_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.m_shock.y(3,:)'; flipud(pctl_in.m_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.m_shock.y(4,:)'; flipud(pctl_in.m_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.m_shock.y(5,:)'; flipud(pctl_in.m_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.m_shock.y(6,:)'; flipud(pctl_in.m_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.m_shock.y(7,:)'; flipud(pctl_in.m_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.m_shock.y(8,:)'; flipud(pctl_in.m_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.m_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.m_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Macro Shock');
title('GDP response')

subplot(2,3,3); hold on;
patch(patchx,[pctl_in.y_shock.y(1,:)'; flipud(pctl_in.y_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.y_shock.y(2,:)'; flipud(pctl_in.y_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.y_shock.y(3,:)'; flipud(pctl_in.y_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.y_shock.y(4,:)'; flipud(pctl_in.y_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.y_shock.y(5,:)'; flipud(pctl_in.y_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.y_shock.y(6,:)'; flipud(pctl_in.y_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.y_shock.y(7,:)'; flipud(pctl_in.y_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.y_shock.y(8,:)'; flipud(pctl_in.y_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.y_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.y_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
lgd = legend('Location','eastoutside'); legend boxoff
ylabel('- GDP Shock');
title('GDP response')

subplot(2,3,4)
plot(pctl_in.p_bad_path_shock.f,'LineWidth',2.5,'Color',dr,'DisplayName','Shocked'); hold on;
plot(pctl_in.p_bad_path_shockbase.f,'--','LineWidth',2,'Color',db,'DisplayName','Baseline'); hold on;
legend('Location','N','Orientation','Horizontal');title('p(s=bad)'); legend boxoff
xticks([3 6 9 12]); box off;
set(gca,'FontSize',14); xlim([1 12]);
ylim([0,1.1]);

subplot(2,3,5)
plot(pctl_in.p_bad_path_shock.m,'LineWidth',2.5,'Color',dr,'DisplayName','Shocked'); hold on;
plot(pctl_in.p_bad_path_shockbase.m,'--','LineWidth',2,'Color',db,'DisplayName','Baseline'); hold on;
legend('Location','N','Orientation','Horizontal');title('p(s=bad)'); legend boxoff
xticks([3 6 9 12]); box off;
set(gca,'FontSize',14); xlim([1 12]);
ylim([0,1.1]);

subplot(2,3,6)
plot(pctl_in.p_bad_path_shock.y,'LineWidth',2.5,'Color',dr,'DisplayName','Shocked'); hold on;
plot(pctl_in.p_bad_path_shockbase.y,'--','LineWidth',2,'Color',db,'DisplayName','Baseline'); hold on;
legend('Location','N','Orientation','Horizontal');title('p(s=bad)'); legend boxoff
xticks([3 6 9 12]); box off;
ylim([0,1.1]); xlim([1 12]);
set(gca,'FontSize',14);

set(lgd,'Position',[0.5 0.06 0.0 0.0],'Orientation','Horizontal'); 

saveas(gcf,[filename]);

%%

figure('Position',[0, 0, 1200, 800]);
subplot(3,3,1); hold on;
patch(patchx,[pctl_in.f_shock.y(1,:)'; flipud(pctl_in.f_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.f_shock.y(2,:)'; flipud(pctl_in.f_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.f_shock.y(3,:)'; flipud(pctl_in.f_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.f_shock.y(4,:)'; flipud(pctl_in.f_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.f_shock.y(5,:)'; flipud(pctl_in.f_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.f_shock.y(6,:)'; flipud(pctl_in.f_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.f_shock.y(7,:)'; flipud(pctl_in.f_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.f_shock.y(8,:)'; flipud(pctl_in.f_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.f_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.f_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Financial Shock');
title('GDP response')

subplot(3,3,2); hold on;
patch(patchx,[pctl_in.m_shock.y(1,:)'; flipud(pctl_in.m_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.m_shock.y(2,:)'; flipud(pctl_in.m_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.m_shock.y(3,:)'; flipud(pctl_in.m_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.m_shock.y(4,:)'; flipud(pctl_in.m_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.m_shock.y(5,:)'; flipud(pctl_in.m_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.m_shock.y(6,:)'; flipud(pctl_in.m_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.m_shock.y(7,:)'; flipud(pctl_in.m_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.m_shock.y(8,:)'; flipud(pctl_in.m_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.m_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.m_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Macro Shock');
title('GDP response')

subplot(3,3,3); hold on;
patch(patchx,[pctl_in.y_shock.y(1,:)'; flipud(pctl_in.y_shock.y(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.y_shock.y(2,:)'; flipud(pctl_in.y_shock.y(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.y_shock.y(3,:)'; flipud(pctl_in.y_shock.y(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.y_shock.y(4,:)'; flipud(pctl_in.y_shock.y(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.y_shock.y(5,:)'; flipud(pctl_in.y_shock.y(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.y_shock.y(6,:)'; flipud(pctl_in.y_shock.y(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.y_shock.y(7,:)'; flipud(pctl_in.y_shock.y(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.y_shock.y(8,:)'; flipud(pctl_in.y_shock.y(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,irf_fixed_good.y_shock.y(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:),'DisplayName','Good Regime');
plot(time,irf_fixed_bad.y_shock.y(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:),'DisplayName','Bad Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
lgd = legend('Location','eastoutside'); legend boxoff
ylabel('- GDP Shock');
title('GDP response')

subplot(3,3,4); hold on;
patch(patchx,[pctl_in.f_shock.f(1,:)'; flipud(pctl_in.f_shock.f(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.f_shock.f(2,:)'; flipud(pctl_in.f_shock.f(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.f_shock.f(3,:)'; flipud(pctl_in.f_shock.f(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.f_shock.f(4,:)'; flipud(pctl_in.f_shock.f(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.f_shock.f(5,:)'; flipud(pctl_in.f_shock.f(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.f_shock.f(6,:)'; flipud(pctl_in.f_shock.f(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.f_shock.f(7,:)'; flipud(pctl_in.f_shock.f(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.f_shock.f(8,:)'; flipud(pctl_in.f_shock.f(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.f_shock.f(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.f_shock.f(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([f_lb f_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Financial Shock');
title('FF response')

subplot(3,3,5); hold on;
patch(patchx,[pctl_in.m_shock.f(1,:)'; flipud(pctl_in.m_shock.f(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.m_shock.f(2,:)'; flipud(pctl_in.m_shock.f(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.m_shock.f(3,:)'; flipud(pctl_in.m_shock.f(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.m_shock.f(4,:)'; flipud(pctl_in.m_shock.f(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.m_shock.f(5,:)'; flipud(pctl_in.m_shock.f(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.m_shock.f(6,:)'; flipud(pctl_in.m_shock.f(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.m_shock.f(7,:)'; flipud(pctl_in.m_shock.f(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.m_shock.f(8,:)'; flipud(pctl_in.m_shock.f(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.m_shock.f(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.m_shock.f(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([f_lb f_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Macro Shock');
title('FF response')

subplot(3,3,6); hold on;
patch(patchx,[pctl_in.y_shock.f(1,:)'; flipud(pctl_in.y_shock.f(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.y_shock.f(2,:)'; flipud(pctl_in.y_shock.f(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.y_shock.f(3,:)'; flipud(pctl_in.y_shock.f(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.y_shock.f(4,:)'; flipud(pctl_in.y_shock.f(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.y_shock.f(5,:)'; flipud(pctl_in.y_shock.f(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.y_shock.f(6,:)'; flipud(pctl_in.y_shock.f(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.y_shock.f(7,:)'; flipud(pctl_in.y_shock.f(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.y_shock.f(8,:)'; flipud(pctl_in.y_shock.f(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.y_shock.f(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.y_shock.f(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([f_lb f_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- GDP Shock');
title('FF response')


subplot(3,3,7); hold on;
patch(patchx,[pctl_in.f_shock.m(1,:)'; flipud(pctl_in.f_shock.m(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.f_shock.m(2,:)'; flipud(pctl_in.f_shock.m(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.f_shock.m(3,:)'; flipud(pctl_in.f_shock.m(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.f_shock.m(4,:)'; flipud(pctl_in.f_shock.m(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.f_shock.m(5,:)'; flipud(pctl_in.f_shock.m(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.f_shock.m(6,:)'; flipud(pctl_in.f_shock.m(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.f_shock.m(7,:)'; flipud(pctl_in.f_shock.m(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.f_shock.m(8,:)'; flipud(pctl_in.f_shock.m(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.f_shock.m(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.f_shock.m(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([m_lb m_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Financial Shock');
title('MM response')

subplot(3,3,8); hold on;
patch(patchx,[pctl_in.m_shock.m(1,:)'; flipud(pctl_in.m_shock.m(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.m_shock.m(2,:)'; flipud(pctl_in.m_shock.m(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.m_shock.m(3,:)'; flipud(pctl_in.m_shock.m(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.m_shock.m(4,:)'; flipud(pctl_in.m_shock.m(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.m_shock.m(5,:)'; flipud(pctl_in.m_shock.m(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.m_shock.m(6,:)'; flipud(pctl_in.m_shock.m(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.m_shock.m(7,:)'; flipud(pctl_in.m_shock.m(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.m_shock.m(8,:)'; flipud(pctl_in.m_shock.m(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.m_shock.m(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.m_shock.m(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([m_lb m_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- Macro Shock');
title('MM response')

subplot(3,3,9); hold on;
patch(patchx,[pctl_in.y_shock.m(1,:)'; flipud(pctl_in.y_shock.m(2,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.y_shock.m(2,:)'; flipud(pctl_in.y_shock.m(3,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.y_shock.m(3,:)'; flipud(pctl_in.y_shock.m(4,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.y_shock.m(4,:)'; flipud(pctl_in.y_shock.m(5,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.y_shock.m(5,:)'; flipud(pctl_in.y_shock.m(6,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.y_shock.m(6,:)'; flipud(pctl_in.y_shock.m(7,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.y_shock.m(7,:)'; flipud(pctl_in.y_shock.m(8,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.y_shock.m(8,:)'; flipud(pctl_in.y_shock.m(9,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,irf_fixed_good.y_shock.m(5,:)','-','LineWidth',2.5,'Color',cmedian(2,:));
plot(time,irf_fixed_bad.y_shock.m(5,:)','--','LineWidth',2.5,'Color',cmedian(1,:));
xticks([3 6 9 12]);  xlim([1 12]);ylim([m_lb m_ub]);
set(gca,'FontSize',14);
hold off;
ylabel('- GDP Shock');
title('MM response')



set(lgd,'Position',[0.5 0.06 0.0 0.0],'Orientation','Horizontal'); 



%%
sgtitle(supt);
orient(gcf,'landscape');
saveas(gcf,['All_' filename]);

% saveas(gcf,[filename]);

end