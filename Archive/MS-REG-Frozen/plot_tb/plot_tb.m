%% OCTOBER DENSITY with QR
%-------------------------------------------------------------------
% Comparison with Quantile Regression Densities
%-------------------------------------------------------------------

clear
close all

model_selec = 2; % 1, for Foreign
                 % 2, for U.S.

if model_selec ==1
    load('data_plot_foreign.mat');
    xmin = -2; xmax = 8;
elseif model_selec ==2
    load('data_plot_us.mat');
    xmin = -6; xmax = 6;
end

colors2 = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;...
    0.4660 0.6740 0.1888];

dataformat  = 'yyyy-mmm';
md = '2020-Oct';
mydate = find(datenum(md,dataformat)==dates_full);

FontSize = 16;
numticks = 48;
figSize = [12 6];

% Select Time Periods for Which to Plot PDFs
period_pdf = {md};

tind = find(datenum(char(period_pdf),dataformat)==dates_full)+1;

% rescale  pdf so that the range is between 0 and 1
PST = ResMatch.PST;
CDF = ResMatch.CST;

[~,i25_oct] = min(abs(CDF(tind,:)'-0.25));
[~,i25_ms_oct] = min(abs(cdf_oct-0.25));

ym = max(PST(tind,:)');
ymax = max([(max(ym)*1.05);max(pdf_oct)*1.05]);

fig=figure;
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
hold on;
l1=plot(ax1,ResMatch.YY,PST(tind,:)','-.','Color',colors2(5,:),'LineWidth',2,'DisplayName',['2020-Oct (QR)']);
l11=plot(ax1,xi_oct,pdf_oct,'-','LineWidth', 2,'Color',colors2(5,:),'DisplayName',['2020-Oct (MS-VAR)']);
l1.Color = [l1.Color 0.5];
plot(ax1,[0 0],[0 ymax],'k--','HandleVisibility','off');
plot(ax1,ResMatch.YY(:,i25_oct),PST(tind,i25_oct)','d','MarkerFaceColor',colors2(5,:),'MarkerEdgeColor','black','Markersize',10);
plot(ax1,xi_oct(i25_ms_oct),pdf_oct(i25_ms_oct),'s','MarkerFaceColor',colors2(5,:),'MarkerEdgeColor','black','Markersize',10);
axis tight
yl = [0 ymax];
hL=legend([l1 l11],'Location','NorthEast');
set(hL,'interpreter','Latex')
legend boxoff
title('QR and MS-VAR: Density for GDP growth over the next 12 months','FontSize',16','Interpreter','Latex');
set(ax1, 'XLim', [xmin, xmax])
set(ax1,'XTick',xmin:1:xmax)
set(ax1,'TickLabelInterpreter','Latex')
set(ax1, 'FontName', 'Times New Roman');
set(ax1, 'FontSize', FontSize);
set(ax1,'Layer','top')
set(fig,'PaperOrientation','portrait');
set(fig, 'PaperSize', figSize);
set(fig, 'PaperUnits', 'inches');
set(fig, 'Units','inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
