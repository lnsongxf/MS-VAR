function plot_fan_gdp_alg10exante_mcmc(pctl_in,periods,filename,supt,pctl2_in)
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

myalpha25 = 0.75;
myalpha10 = 0.65;
myalpha5 = 0.55;
myalpha1 = 0.4;


% Set y-axis bounds
y_ub = 0.5;
y_lb = -2;
   
% Get percentiles 
min_y = min([min(min(pctl_in.f_shock.y)); min(min(pctl_in.m_shock.y)); min(min(pctl_in.y_shock.y))]);
max_y = max([max(max(pctl_in.f_shock.y)); max(max(pctl_in.m_shock.y)); max(max(pctl_in.y_shock.y))]);

if min_y < y_lb
    y_lb = floor(min_y);
end

if max_y > y_ub
    y_ub = ceil(max_y);
end



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
plot(time,pctl2_in.f_shock.y(5,:)','-','LineWidth',2.5,'Color','k','DisplayName','Fixed Regime');
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
plot(time,pctl2_in.m_shock.y(5,:)','-','LineWidth',2.5,'Color','k','DisplayName','Fixed Regime');
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
plot(time,pctl2_in.y_shock.y(5,:)','-','LineWidth',2.5,'Color','k','DisplayName','Fixed Regime');
xticks([3 6 9 12]);  xlim([1 12]);ylim([y_lb y_ub]);
set(gca,'FontSize',14);
hold off;
lgd = legend('Location','eastoutside'); legend boxoff
ylabel('- GDP Shock');
title('GDP response')

subplot(2,3,4); hold on;
patch(patchx,[pctl_in.p_bad_path_shock.f(1,:)'; flipud(pctl_in.p_bad_path_shock.f(2,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1); 
patch(patchx,[pctl_in.p_bad_path_shock.f(2,:)'; flipud(pctl_in.p_bad_path_shock.f(3,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.p_bad_path_shock.f(3,:)'; flipud(pctl_in.p_bad_path_shock.f(4,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.p_bad_path_shock.f(4,:)'; flipud(pctl_in.p_bad_path_shock.f(5,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.p_bad_path_shock.f(5,:)'; flipud(pctl_in.p_bad_path_shock.f(6,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25)
patch(patchx,[pctl_in.p_bad_path_shock.f(6,:)'; flipud(pctl_in.p_bad_path_shock.f(7,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10)
patch(patchx,[pctl_in.p_bad_path_shock.f(7,:)'; flipud(pctl_in.p_bad_path_shock.f(8,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5)
patch(patchx,[pctl_in.p_bad_path_shock.f(8,:)'; flipud(pctl_in.p_bad_path_shock.f(9,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1)
plot(time,pctl2_in.p_bad_path_shock.f(5,:)','-','LineWidth',2.5,'Color','k');
b=plot(time,pctl2_in.p_bad_path_shockbase.f(5,:)','--','LineWidth',2.5,'Color',db,'DisplayName','Baseline Mean');
legend(b,'Location','north');
xticks([3 6 9 12]); box off;
ylim([0,1.1]);
set(gca,'FontSize',14); xlim([1 12]);
hold off;
title('p(s=bad)');


subplot(2,3,5); hold on;
patch(patchx,[pctl_in.p_bad_path_shock.m(1,:)'; flipud(pctl_in.p_bad_path_shock.m(2,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.p_bad_path_shock.m(2,:)'; flipud(pctl_in.p_bad_path_shock.m(3,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.p_bad_path_shock.m(3,:)'; flipud(pctl_in.p_bad_path_shock.m(4,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.p_bad_path_shock.m(4,:)'; flipud(pctl_in.p_bad_path_shock.m(5,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.p_bad_path_shock.m(5,:)'; flipud(pctl_in.p_bad_path_shock.m(6,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.p_bad_path_shock.m(6,:)'; flipud(pctl_in.p_bad_path_shock.m(7,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.p_bad_path_shock.m(7,:)'; flipud(pctl_in.p_bad_path_shock.m(8,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.p_bad_path_shock.m(8,:)'; flipud(pctl_in.p_bad_path_shock.m(9,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,pctl2_in.p_bad_path_shock.m(5,:)','-','LineWidth',2.5,'Color','k','DisplayName','Fixed Regime');
b=plot(time,pctl2_in.p_bad_path_shockbase.m(5,:)','--','LineWidth',2.5,'Color',db,'DisplayName','Baseline Mean');
legend(b,'Location','north');
xticks([3 6 9 12]); box off;
ylim([0,1.1]);
set(gca,'FontSize',14); xlim([1 12]);
hold off;
title('p(s=bad)');

subplot(2,3,6); hold on;
patch(patchx,[pctl_in.p_bad_path_shock.y(1,:)'; flipud(pctl_in.p_bad_path_shock.y(2,:)')],lb99,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','1st'); 
patch(patchx,[pctl_in.p_bad_path_shock.y(2,:)'; flipud(pctl_in.p_bad_path_shock.y(3,:)')],lb95,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','5th')
patch(patchx,[pctl_in.p_bad_path_shock.y(3,:)'; flipud(pctl_in.p_bad_path_shock.y(4,:)')],lb90,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','10th')
patch(patchx,[pctl_in.p_bad_path_shock.y(4,:)'; flipud(pctl_in.p_bad_path_shock.y(5,:)')],db,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','25th')
patch(patchx,[pctl_in.p_bad_path_shock.y(5,:)'; flipud(pctl_in.p_bad_path_shock.y(6,:)')],dr,'EdgeColor','none','FaceAlpha',myalpha25,'DisplayName','75th')
patch(patchx,[pctl_in.p_bad_path_shock.y(6,:)'; flipud(pctl_in.p_bad_path_shock.y(7,:)')],lr10,'EdgeColor','none','FaceAlpha',myalpha10,'DisplayName','90th')
patch(patchx,[pctl_in.p_bad_path_shock.y(7,:)'; flipud(pctl_in.p_bad_path_shock.y(8,:)')],lr5,'EdgeColor','none','FaceAlpha',myalpha5,'DisplayName','95th')
patch(patchx,[pctl_in.p_bad_path_shock.y(8,:)'; flipud(pctl_in.p_bad_path_shock.y(9,:)')],lr1,'EdgeColor','none','FaceAlpha',myalpha1,'DisplayName','99th')
plot(time,pctl2_in.p_bad_path_shock.y(5,:)','-','LineWidth',2.5,'Color','k','DisplayName','Fixed Regime');
b=plot(time,pctl2_in.p_bad_path_shockbase.y(5,:)','--','LineWidth',2.5,'Color',db,'DisplayName','Baseline Mean');
legend(b,'Location','north');
xticks([3 6 9 12]); box off;
ylim([0,1.1]); xlim([1 12]);
set(gca,'FontSize',14);
hold off;
title('p(s=bad)');

set(lgd,'Position',[0.5 0.06 0.0 0.0],'Orientation','Horizontal'); 
%%
sgtitle(supt);
orient(gcf,'landscape');
saveas(gcf,[filename]);

end