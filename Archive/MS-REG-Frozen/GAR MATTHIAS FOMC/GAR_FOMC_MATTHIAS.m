%% MS-VAR GAR DENISITIES
% March 2020 vs April 2021
% Dario Caldara, Danilo Cascaldi-Garcia, Pablo Cuba-Borda, Francesca Loria

clear; clc; close all

FontSize = 16;
figSize = [12 6];


staff = 2.14; % Q3/Q3
s_wave = -2.45; % Q3/Q3

load('GAR_FOMC_MATTHIAS.mat')

fig=figure;
hold on
l1=plot(xi_mar,pdf_mar,'Color',[0.6471 0 0.1490],'LineWidth', 3,'DisplayName',['March-2020 ']); hold on;
l2=plot(xi_apr21,pdf_apr21,'Color',[0.2392 0.3843 0.6667],'LineWidth', 3,'DisplayName',['April-2021 ']);
ylimits = ylim; xlimtis = xlim;
plot(zeros(10,1),linspace(0,ylimits(2),10),'k--');
%l3=plot([staff staff],[0 ylimits(2)],'b','LineWidth',1,'DisplayName',['Staff forecast = ' num2str(round(staff,1)) '\%']);
%l4=plot([s_wave s_wave],[0 ylimits(2)],'r','LineWidth',1,'DisplayName',['Second Waves = ' num2str(round(s_wave,1)) '\%']);
hleg = legend([l1 l2 ],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
title('MS-VAR: Densitites for GDP growth over the next 12 months (U.S.)','FontSize',16','Interpreter','Latex');
ylabel('PDF','fontsize',10,'interpreter','Latex')
xlabel('1-Year-Ahead GDP Growth','fontsize',10,'Interpreter','Latex')
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', FontSize);
set(gca,'Layer','top')
set(gca,'XTick',-16:2:6)
set(gca,'TickLabelInterpreter','Latex')
axis tight
xlim([-12 8])
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);

% save('GAR_FOMC_MATTHIAS.mat','xi_mar','pdf_mar','xi_apr21','pdf_apr21')