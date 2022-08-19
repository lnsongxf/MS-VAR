%%***************************************************************
%                         HOUSEKEEPING                          *
%%***************************************************************
clear variables;
close all;
clc;

load 'results.mat';

Horizon=60;

Ltildeq50=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 50th
Ltildeq16=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 16th
Ltildeq84=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 84th

Ltildeq025=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 025th
Ltildeq975=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 975th


for ii=1:size(Ltilde,1)
    for jj=1:size(Ltilde,2)
        for kk=1:size(Ltilde,3)
        Ltildeq50(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.5);
        Ltildeq16(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.16);
        Ltildeq84(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.84);
        
        Ltildeq025(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.025);
        Ltildeq975(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.975);
        end
    end
end

ftsizeaxis=11;
ftsizexlabel=11;
ftsizetitle=11;
ftlinewidth = 1.0;
medianwidth=1.0;

H=Horizon;



close all

hFig = figure(1);
set(hFig, 'Position', [20 20 700 350])

subplot(2,3,1)
plot(0:1:H,squeeze(Ltildeq50(:,1,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on

a=(squeeze(Ltildeq16(:,1,1)))';
b=(squeeze(Ltildeq84(:,1,1)))';

a95=(squeeze(Ltildeq025(:,1,1)))';
b95=(squeeze(Ltildeq975(:,1,1)))';

x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth )
set(gca,'YTick',[-1.2 -0.8 -0.4 0 0.4 0.8])
axis([0 60 -1.2 0.8])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Industrial Production','FontSize',ftsizetitle)
box on


subplot(2,3,2)
plot(0:1:H,squeeze(Ltildeq50(:,2,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on
a=(squeeze(Ltildeq16(:,2,1)))';
b=(squeeze(Ltildeq84(:,2,1)))';

a95=(squeeze(Ltildeq025(:,2,1)))';
b95=(squeeze(Ltildeq975(:,2,1)))';

x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
set(gca,'YTick',[-0.5 -0.25 0 0.25 0.5])
axis([0 60 -0.5 0.5])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' CPI','FontSize',ftsizetitle)
box on

subplot(2,3,3)
plot(0:1:H,squeeze(Ltildeq50(:,3,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on
a=(squeeze(Ltildeq16(:,3,1)))';
b=(squeeze(Ltildeq84(:,3,1)))';
a95=(squeeze(Ltildeq025(:,3,1)))';
b95=(squeeze(Ltildeq975(:,3,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
set(gca,'YTick',[-4 -2 0 2 4])
axis([0 60 -4 4])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Commodity Price Index','FontSize',ftsizetitle)
box on


subplot(2,3,4)
plot(0:1:H,squeeze(Ltildeq50(:,4,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on
a=(squeeze(Ltildeq16(:,4,1)))';
b=(squeeze(Ltildeq84(:,4,1)))';
a95=(squeeze(Ltildeq025(:,4,1)))';
b95=(squeeze(Ltildeq975(:,4,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
set(gca,'YTick',[-4 -2 0 2 4])
axis([0 60 -4 4])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Total Reserves','FontSize',ftsizetitle)
box on


subplot(2,3,5)
plot(0:1:H,squeeze(Ltildeq50(:,5,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on
a=(squeeze(Ltildeq16(:,5,1)))';
b=(squeeze(Ltildeq84(:,5,1)))';
a95=(squeeze(Ltildeq025(:,5,1)))';
b95=(squeeze(Ltildeq975(:,5,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth )
set(gca,'YTick',[-4 -2 0 2 4])
axis([0 60 -4 4])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Nonborrowed Reserves','FontSize',ftsizetitle)
box on



subplot(2,3,6)
plot(0:1:H,squeeze(Ltildeq50(:,6,1)),'-k','LineWidth',medianwidth)
hline(0,'-r')
hold on
a=(squeeze(Ltildeq16(:,6,1)))';
b=(squeeze(Ltildeq84(:,6,1)))';
a95=(squeeze(Ltildeq025(:,6,1)))';
b95=(squeeze(Ltildeq975(:,6,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb('mediumblue'),rgb('mediumblue'),0,0.5);
hold on
[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
set(gca,'XTick',[0;12;24;36;48;60])
set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth )
set(gca,'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4])
axis([0 60 -0.3 0.4])
xlabel('Years','FontSize',ftsizexlabel)
ylabel('Percentage points','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Federal Funds Rate','FontSize',ftsizetitle)
box on

set(gcf, 'PaperPositionMode', 'auto');



print -dpng  'restrictions_1_and_2_in_first_differences_1983m1_2007m6_IP_CPI.png'

