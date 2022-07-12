set(0,'defaulttextinterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');


fig=figure('Name','Histograms'); clf;
subplot(3,2,1)
hold on
l1=histogram(s_3_3_sync_1(correct,:),50,'EdgeColor',colors(60,:),'FaceColor',colors(60,:),'DisplayName','S.d. regime 1');
l2=histogram(s_3_3_sync_2(correct,:),50,'EdgeColor',colors(25,:),'FaceColor',colors(25,:),'DisplayName','S.d. regime 2');
if dir_or_it ==2
    l3=xline(dgp.s_3_3_sync_1,'Color','k','LineWidth', 1.5,'DisplayName','DGPs');
    xline(dgp.s_3_3_sync_2,'Color','k','LineWidth', 1.5);
    legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
elseif dir_or_it ==1
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
end
set(gca,'children',flipud(get(gca,'children')))
hold off
axis tight
subplot(3,2,2)
hold on
l1=histogram(a0_3_1_sync_1(correct,:)*(-1),50,'EdgeColor',colors(60,:),'FaceColor',colors(60,:),'DisplayName','$f_t$ regime 1');
l2=histogram(a0_3_1_sync_2(correct,:)*(-1),50,'EdgeColor',colors(25,:),'FaceColor',colors(25,:),'DisplayName','$f_t$ regime 2');
if dir_or_it ==2
    l3=xline(dgp.a0_3_1_sync_1,'Color','k','LineWidth', 1.5,'DisplayName','DGPs');
    xline(dgp.a0_3_1_sync_2,'Color','k','LineWidth', 1.5);
    legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
elseif dir_or_it ==1
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
end
set(gca,'children',flipud(get(gca,'children')))
hold off
axis tight
subplot(3,2,3)
hold on
l1=histogram(c_3_1_sync_1(correct,:),50,'EdgeColor',colors(60,:),'FaceColor',colors(60,:),'DisplayName','Constant regime 1');
l2=histogram(c_3_1_sync_2(correct,:),50,'EdgeColor',colors(25,:),'FaceColor',colors(25,:),'DisplayName','Constant regime 2');
if dir_or_it ==2
    l3=xline(dgp.c_3_1_sync_1,'Color','k','LineWidth', 1.5,'DisplayName','DGPs');
    xline(dgp.c_3_1_sync_2,'Color','k','LineWidth', 1.5);
    legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
elseif dir_or_it ==1
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
end
set(gca,'children',flipud(get(gca,'children')))
hold off
axis tight
subplot(3,2,4)
hold on
l1=histogram(a0_3_2_sync_1(correct,:)*(-1),50,'EdgeColor',colors(60,:),'FaceColor',colors(60,:),'DisplayName','$m_t$ regime 1');
l2=histogram(a0_3_2_sync_2(correct,:)*(-1),50,'EdgeColor',colors(25,:),'FaceColor',colors(25,:),'DisplayName','$m_t$ regime 2');
if dir_or_it ==2
    l3=xline(dgp.a0_3_2_sync_1,'Color','k','LineWidth', 1.5,'DisplayName','DGPs');
    xline(dgp.a0_3_2_sync_2,'Color','k','LineWidth', 1.5);
    legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
elseif dir_or_it ==1
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
end
set(gca,'children',flipud(get(gca,'children')))
hold off
axis tight
tightfig;
if saveit==1
    print('-dpdf',fig,[sim_folder '/' model 'Histograms'],'-bestfit');
    saveas(fig,sprintf('%s.png',[sim_folder '/' model 'Histograms']));
end


if dir_or_it ==2
       
    % Contemporaneous
    fig=figure; clf;
    hold on
    l1=plot(1:size_y, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
    l2=plot(1:size_y, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
    l3=plot(1:size_y, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
    l4=plot(1:size_y, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
    legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Estimated GDP Distribution in $t$','FontSize',16','Interpreter','Latex');
    axis tight
    if saveit==1
        print('-dpdf',fig,[sim_folder '/' model 'GDP_estimated'],'-bestfit');
        saveas(fig,sprintf('%s.png',[sim_folder '/' model 'GDP_estimated']));
    end
    
    % 12-months ahead
    fig=figure; clf;
    hold on
    l1=plot(1:size_y-12, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
    l2=plot(1:size_y-12, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
    l3=plot(1:size_y-12, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
    l4=plot(1:size_y-12, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
    
    legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Estimated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
    axis tight
    if saveit==1
        print('-dpdf',fig,[sim_folder '/' model 'GDP_estimated_one_year'],'-bestfit');
        saveas(fig,sprintf('%s.png',[sim_folder '/' model 'GDP_estimated_one_year']));
    end
    
    % Comparison 12-months ahead
    fig=figure; clf;
    hold on
    l1=plot(1:size_y-12, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
    % l2=plot(1:size_y-12, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
    % l3=plot(1:size_y-12, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
    l4=plot(1:size_y-12, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
    l11=plot(1:size_y-12, Res.Res_direct.dY_10(2:end-12),'-.','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','DGP 10th');
    % l22=plot(1:size_y-12, Res.Res_direct.dY_25(2:end-12),'-.','Color',colors(15,:),'LineWidth', 2,'DisplayName','DGP 25th');
    % l33=plot(1:size_y-12, Res.Res_direct.dY_75(2:end-12),'-.','Color',colors(45,:),'LineWidth', 2,'DisplayName','DGP 75th');
    l44=plot(1:size_y-12, Res.Res_direct.dY_90(2:end-12),'-.','Color',colors(55,:),'LineWidth', 1.5,'DisplayName','DGP 90th');
    legend([l1 l4 l11 l44],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Estimated and DGP of GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
    axis tight
    if saveit==1
        print('-dpdf',fig,[sim_folder '/' model 'GDP_estimated_and_dgp'],'-bestfit');
        saveas(fig,sprintf('%s.png',[sim_folder '/' model 'GDP_estimated_and_dgp']));
    end
    
elseif dir_or_it ==1
    % 12-months ahead
    fig=figure; clf;
    hold on
    l1=plot(1:size_y, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
    l2=plot(1:size_y, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
    l3=plot(1:size_y, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
    l4=plot(1:size_y, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
    
    legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Simulated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
    axis tight
    
    % Comparison 12-months ahead
    fig=figure; clf;
    hold on
    l1=plot(1:size_y, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
    % l2=plot(1:size_y, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
    % l3=plot(1:size_y, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
    l4=plot(1:size_y, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
    l11=plot(1:size_y, Res.Res_direct.dY_10(2:end),'-.','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','DGP 10th');
    % l22=plot(1:size_y, Res.Res_direct.dY_25(2:end),'-.','Color',colors(15,:),'LineWidth', 2,'DisplayName','DGP 25th');
    % l33=plot(1:size_y, Res.Res_direct.dY_75(2:end),'-.','Color',colors(45,:),'LineWidth', 2,'DisplayName','DGP 75th');
    l44=plot(1:size_y, Res.Res_direct.dY_90(2:end),'-.','Color',colors(55,:),'LineWidth', 1.5,'DisplayName','DGP 90th');
    legend([l1 l4 l11 l44],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Estimated and DGP of GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
    axis tight
    if saveit==1
        print('-dpdf',fig,[sim_folder '/' model 'GDP_estimated_and_dgp'],'-bestfit');
        saveas(fig,sprintf('%s.png',[sim_folder '/' model 'GDP_estimated_and_dgp']));
    end
end








