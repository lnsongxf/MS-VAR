%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
dir_or_it = 2; % 1 = direct, 2 = iterated
compare   = 1; % 1 = compare the interated and direct results
const     = 1; % 1 = have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution
saveit    = 0; % 1 = save figures

% Data vintage, sample and country selection
start_date   = '2100-Jan';
end_date     = '2198-Dec';
end_date_db  = '2199-Dec';

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

if dir_or_it ==1
    model = 'results_direct';
elseif dir_or_it ==2
    model = 'results_iterated';
end

sim_folder = 'Results/Simulation/';
if exist(sim_folder,'dir')==0
    mkdir(sim_folder)
end


% Create date formats for plotting
inputformat = 'yyyy-MMM';
dataformat  = 'yyyy-mmm';
start_plot    = '2100-Feb';
end_plot      = end_date_db;
end_plot_full = end_date;

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Vector of dates for plotting
dates   = dates_full(sd:ed);

% Figure options
FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);
left_color = [0 0 0];
right_color = colors(2,:);

%% Load RISE
rise_startup()

%% Prepare Data
if exist([sim_folder '' model '.mat'],'file')==2
    compare_temp = compare;
    load([sim_folder '' model '.mat'])
    compare = compare_temp;
else
    run estimation_loop.m;
    save([sim_folder '' model '.mat'])
end

%% Check if reg =1 is high growth, low vol
check_good = zeros(n_draws,1);
check_flip = zeros(n_draws,1);
check_wrong = zeros(n_draws,1);

for yy=1:n_draws
    if c_3_1_sync_1(yy)>c_3_1_sync_2(yy)
        if s_3_3_sync_1(yy)<s_3_3_sync_2(yy)
            check_good(yy)=1; % if mean is larger, and s.d. is smaller
        else
            check_wrong(yy)=1; % if mean is larger, and s.d. is larger
        end
    else
        if s_3_3_sync_1(yy)<s_3_3_sync_2(yy)
            check_wrong(yy)=1; % if mean is smaller, and s.d. is smaller
        else
            check_flip(yy)=1; % if mean is smaller, and s.d. is larger
        end
    end
end
%% Quantiles

y_mat = y_fit';
size_y = size(y_mat,2);

% Check for convergence errors
correct = s_3_3_sync_2<dgp.s_3_3_sync_2*10;


% Compute percentiles
dY_25 = prctile(y_mat(correct,:),25)'; dY_75 = prctile(y_mat(correct,:),75)';
dY_10 = prctile(y_mat(correct,:),10)'; dY_90 = prctile(y_mat(correct,:),90)';

run plot_histograms.m

if dir_or_it ==2
    
    % 12-months ahead from the fitted
    y_mat_fut = NaN(size(y_mat,1),size(y_mat,2)-12);
    for hh=1:size(y_mat,1)
        for ww=1:size(y_mat,2)-12
            y_mat_fut(hh,ww) = mean(y_mat(hh,ww+1:ww+12));
        end
    end
    dY_25_fut = prctile(y_mat_fut(correct,:),25)'; dY_75_fut = prctile(y_mat_fut(correct,:),75)';
    dY_10_fut = prctile(y_mat_fut(correct,:),10)'; dY_90_fut = prctile(y_mat_fut(correct,:),90)';
    
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

%% Compare both estimations
if compare ==1
    if exist([sim_folder 'results_iterated.mat'],'file')==2
        load([sim_folder 'results_iterated.mat'])
        compare =1;
        st_mat = st_fit';
        y_mat = y_fit';
        % 12-months ahead iterated
        y_mat_fut2 = NaN(size(y_mat,1),size(y_mat,2));
        % Set a 0.9 threshold to switch from good to bad regime
        % Check for convergence errors
        correct = s_3_3_sync_2<dgp.s_3_3_sync_2*10;
        
        for jj=1:size(y_mat,2)
            for dd=1:n_draws
                y_temp = [y_mat(dd,jj);NaN(12,1)];
                f_temp = [Res_iterated.FF(dd,jj);NaN(12,1)];
                m_temp = [Res_iterated.MF(dd,jj);NaN(12,1)];
                s_temp = [st_mat(dd,jj);NaN(12,1)];
                
                for ee=2:13
                    
                    eta1 = randn(2,1); % financial and macro shocks
                    
                    m_temp(ee) = c_2_1(dd) +                              a1_2_1(dd)*f_temp(ee-1) + a1_2_2(dd)*m_temp(ee-1) + s_2_2(dd)*eta1(2,1);
                    f_temp(ee) = c_1_1(dd) + a0_1_2(dd)*m_temp(ee)*(-1) + a1_1_1(dd)*f_temp(ee-1) + a1_1_2(dd)*m_temp(ee-1) + s_1_1(dd)*eta1(1,1);
                    
                    
                    
                    p12 = 1/(1+exp(a12(dd)-b12(dd)*(f_temp(ee))+c12(dd)*(m_temp(ee))));
                    p21 = 1/(1+exp(a21(dd)+b21(dd)*(f_temp(ee))-c21(dd)*(m_temp(ee))));
                    
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                    
                    udraw = rand(1);
                    
                    if s_temp(ee-1) == 1 % started in good regime
                        if udraw > p11
                            s_temp(ee) = 2; % switch from good to bad
                        else
                            s_temp(ee) = 1; % don't switch and remain in good
                        end
                    else % start in bad regime
                        if udraw > p22
                            s_temp(ee) = 1; % switch from bad to good
                        else
                            s_temp(ee) = 2; % don't switch and remain in bad
                        end
                        
                    end
                    
                    
                    eta2 = randn(1,1); % GDP shock
                    
                    if s_temp(ee)==2
                        y_temp(ee) = c_3_1_sync_2(dd)+a0_3_1_sync_2(dd)*f_temp(ee)*(-1) +a0_3_2_sync_2(dd)*m_temp(ee)*(-1) +s_3_3_sync_2(dd)*eta2;
                    else
                        y_temp(ee) = c_3_1_sync_1(dd)+a0_3_1_sync_1(dd)*f_temp(ee)*(-1) +a0_3_2_sync_1(dd)*m_temp(ee)*(-1) +s_3_3_sync_1(dd)*eta2;
                    end
                    
                    
                    if ee==13
                        y_mat_fut2(dd,jj) = mean(y_temp(2:end,1));
                    end
                end
            end
        end
        
        dY_25_fut2 = prctile(y_mat_fut2(correct,:),25)'; dY_75_fut2 = prctile(y_mat_fut2(correct,:),75)';
        dY_10_fut2 = prctile(y_mat_fut2(correct,:),10)'; dY_90_fut2 = prctile(y_mat_fut2(correct,:),90)';
        if exist([sim_folder 'results_direct.mat'],'file')==2
            load([sim_folder 'results_direct.mat'])
            y_mat = y_fit';
            size_y = size(y_mat,2);
            
            % Check for convergence errors
            correct = s_3_3_sync_2<dgp.s_3_3_sync_2*10;
            
            % Compute percentiles
            dY_25 = prctile(y_mat(correct,:),25)'; dY_75 = prctile(y_mat(correct,:),75)';
            dY_10 = prctile(y_mat(correct,:),10)'; dY_90 = prctile(y_mat(correct,:),90)';
            
            fig=figure; clf;
            hold on
            l1=plot(1:size_y, dY_10_fut2,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Iterated 10th');
            % l2=plot(1:size_y, dY_25_fut2,'Color',colors(15,:),'LineWidth', 2,'DisplayName','Iterated 25th');
            % l3=plot(1:size_y, dY_75_fut2,'Color',colors(45,:),'LineWidth', 2,'DisplayName','Iterated 75th');
            l4=plot(1:size_y, dY_90_fut2,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','Iterated 90th');
            l11=plot(1:size_y, dY_10,'-.','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Direct 10th');
            % l22=plot(1:size_y, dY_25,'-.','Color',colors(15,:),'LineWidth', 2,'DisplayName','Direct 25th');
            % l33=plot(1:size_y, dY_75,'-.','Color',colors(45,:),'LineWidth', 2,'DisplayName','Direct 75th');
            l44=plot(1:size_y, dY_90,'-.','Color',colors(55,:),'LineWidth', 1.5,'DisplayName','Direct 90th');
            legend([l1 l4 l11 l44],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Comparison GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
            axis tight
            if saveit==1
                print('-dpdf',fig,[sim_folder 'CompareGDP_dist_direct_iterated'],'-bestfit');
                saveas(fig,sprintf('%s.png',[sim_folder 'CompareGDP_dist_direct_iterated']));
            end
            
            example_draw = 10;
            
            fig=figure; clf;
            hold on
            l11=plot(1:size_y, y_mat(example_draw,:),'Color',colors(60,:),'LineWidth', 1.5,'DisplayName','Direct draw');
            l1=plot(1:size_y, y_mat_fut2(example_draw,:),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Iterated draw');
            legend([l1 l11],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Comparison GDP distribution 1-year ahead (single draw)','FontSize',16','Interpreter','Latex');
            axis tight
            if saveit==1
                print('-dpdf',fig,[sim_folder 'CompareGDP_draw_direct_iterated'],'-bestfit');
                saveas(fig,sprintf('%s.png',[sim_folder 'CompareGDP_draw_direct_iterated']));
            end
        end
    end
end


toc;
