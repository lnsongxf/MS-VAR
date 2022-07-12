%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
dir_or_it    = 2; % 1 = direct, 2 = iterated
const        = 1; % 1 = have a constant in transition probability
normal       = 0; % 1 = use normal distribution, 0 = gamma distribution
saveit       = 0; % 1 = save figures
example_draw = 10; % single draw for comparison
num_rep      = 5000; % Number of repetitions for the coverage bands of the DGP
tperiods     = 312;  % Number of time-periods in simulation of GDP
nParamDraws  = 10000; %samp.options.N;    % Number of parameter draws
nDraws       = 10;   % Number of simulations for each parameter draw

compare      = 2;    % 1 = compare the iterated and direct results (no parameter uncertainty)
% 2 = compare the iterated and direct results (parameter uncertainty)


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

%% Load results or estimate
if exist([sim_folder '' model '.mat'],'file')==2
    compare_temp = compare;
    load([sim_folder '' model '.mat'])
    compare = compare_temp;
else
    run estimation_loop.m; % This only estimates and stores the posterior mode
    save([sim_folder '' model '.mat'])
end

if compare ==2
    if exist([sim_folder '' model '_samp.mat'],'file')==2
        load([sim_folder '' model '_samp.mat'])
    else
        samp=posterior_sampling_example(example_draw,dir_or_it,nlags,start_date,const,normal);
        save([sim_folder '' model '_samp.mat'],'samp')
    end
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
end

%run plot_histograms_and_figs.m


%% PARAMETER UNCERTAINTY
close all;

if compare ==2
    if exist([sim_folder 'results_iterated.mat'],'file')==2
        load([sim_folder 'results_iterated.mat'])
        if exist([sim_folder 'results_iterated_samp.mat'],'file')==2;load([sim_folder 'results_iterated_samp.mat']);end
        
        
        %------------------------------------------------------------------
        % DGP
        %------------------------------------------------------------------
        
        
        % Initialize matrices
        y_mat_fut_ex_dgp = NaN(tperiods,1); % This stores all realizations for the ITERATED
        st_dgp           = NaN(tperiods,1); % This stores all realizations for the ITERATED
        st_dgp_t         = NaN(tperiods,1); % This stores all realizations for the ITERATED
        wb               = waitbar(0,'Calculating paths for the d.g.p.');
        
        for jj=1:tperiods
            waitbar(jj/tperiods,wb,'Calculating paths for the d.g.p.')
            
            dd=example_draw;
            
            for uu=1:num_rep
                
                % Macro and financial factor in period t
                f_t = Res_iterated.FF(dd,jj);
                m_t = Res_iterated.MF(dd,jj);
                
                % Transition probabilities
                p12 = 1/(1+exp(dgp.a12-dgp.b12*(f_t)+dgp.c12*(m_t)));
                p21 = 1/(1+exp(dgp.a21+dgp.b21*(f_t)-dgp.c21*(m_t)));
                
                % State probabilities t
                p11 = 1 - p12; % probability of remaining in normal
                p22 = 1 - p21; % probability of remaining in bad
                
                % Simulate s(t) conditional on m(t) and f(t)
                
                if jj==1
                    st_lag = 1; % Initialize in good regime
                    
                    st_sim = simulate_st(st_lag,p11,p22,1);
                    
                else
                    
                    st_lag = st_dgp_t(uu,jj-1);
                    
                    st_sim = simulate_st(st_lag,p11,p22,1);
                    
                end
                
                % This is s(t)
                st_temp_par = st_sim;
                
                
                % Simulate factors and transtiion probabilities
                % conditional on f_t and m_t
                [f_temp,m_temp,p11,p22] = simulate_factors(f_t, m_t, dgp, 12,1,'dgp');
                
                % Simulate states: s(t+1)...s(t+h)
                s_temp = simulate_st(st_temp_par,p11,p22,12);
                
                % Simulate GDP
                y_temp = simulate_gdp(s_temp,f_temp,m_temp,dgp,1,'dgp');
                
                % Average future GDP t+1:t+h
                y_mat_fut_ex_dgp(uu,jj) = mean(y_temp(1:end,1));
                
                % Collect t+h+1 state
                st_dgp(uu,jj)   = s_temp(12);
                
                % Collect t state
                st_dgp_t(uu,jj) = st_sim;
            end
            
        end
        close(wb);
        
        % Change st=2 (bad) to st = 0 for plotting
        % Mean(st) represents probability of good regime
        st_dgp(st_dgp==2) = 0;
        
        st_dgp_t(st_dgp_t==2) = 0;
        
        %%
        % Compute percentiles
        st_mean.dgp   = mean(st_dgp)';   % This is in t+12
        st_t_mean.dgp = mean(st_dgp_t)'; % This is in t
        dYsim_25.dgp = prctile(y_mat_fut_ex_dgp,25)';
        dYsim_75.dgp = prctile(y_mat_fut_ex_dgp,75)';
        dYsim_10.dgp = prctile(y_mat_fut_ex_dgp,10)';
        dYsim_90.dgp = prctile(y_mat_fut_ex_dgp,90)';
        
        if saveit
            % Plot figure
            fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(1:tperiods, dYsim_25.dgp,'--','Color',colors(5,:),'LineWidth', 1.5,'DisplayName', '25th');
            l2=plot(1:tperiods, dYsim_75.dgp,'--','Color',colors(50,:),'LineWidth', 1.5,'DisplayName','75th');
            l3=plot(1:tperiods, dYsim_10.dgp,'-','Color',colors(7,:),'LineWidth', 1.5,'DisplayName','10th');
            l4=plot(1:tperiods, dYsim_90.dgp,'-','Color',colors(45,:),'LineWidth', 1.5,'DisplayName','90th');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('GDP distribution 1-year ahead (DGP)','FontSize',16','Interpreter','Latex');
            axis tight
            
            subplot(212)
            hold on
            l1=plot(1:tperiods, 1-st_t_mean.dgp(1:tperiods),'--','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','DGP (mean)');
            legend([l1],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime (single draw)','FontSize',16','Interpreter','Latex');
            axis tight
        end
        %%
        
        
        %------------------------------------------------------------------
        % Iterated Draw: Estimated Parameter
        %------------------------------------------------------------------
        
        
        % Map stored objects
        st_mat = st_fit'; % Now use the fitted values
        y_mat = y_fit'; % Now use the fitted values
        
        % Pre-allocate
        y_mat_fut_ex_temp = NaN(nDraws*nParamDraws,tperiods); % This stores all realizations for the ITERATED
        st_it_temp        = NaN(nDraws*nParamDraws,tperiods);
        st_it_temp_t      = NaN(nDraws*nParamDraws,tperiods);
        
        
        % Waitbar and counters
        wb = waitbar(0,'Calculating paths for the iterated version');
        rcount = 1;
        
        % Loop over parameter draws
        for dd=1:nParamDraws
            
            % Counter for rows storing parameter*repetition
            waitbar(dd/nParamDraws,wb,'Calculating paths for the iterated version')
            
            % Loop over repetitions
            for uu=1:nDraws
                
                % Loop over periods
                for jj=1:tperiods
                    
                    % Macro and financial factor in period t
                    f_t = Res_iterated.FF(example_draw,jj);
                    m_t = Res_iterated.MF(example_draw,jj);
                    
                    % Transition probabilities
                    p12 = 1/(1+exp(samp.a12(dd)-samp.b12(dd)*(f_t)+samp.c12(dd)*(m_t)));
                    p21 = 1/(1+exp(samp.a21(dd)+samp.b21(dd)*(f_t)-samp.c21(dd)*(m_t)));
                    
                    
                    % State probabilities t
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                    % Simulate s(t) conditional on m(t) and f(t)
                    
                    if jj==1
                        st_lag = 1; % Initialize in good regime
                        
                        st_sim = simulate_st(st_lag,p11,p22,1);
                        
                    else
                        
                        st_lag = st_it_temp_t(rcount,jj-1);
                        
                        st_sim = simulate_st(st_lag,p11,p22,1);
                        
                    end
                    
                    % This is s(t)
                    st_temp_par = st_sim;
                    
                    % Simulate factors and transtion probabilities
                    % conditional on f_t and m_t
                    [f_temp,m_temp,p11_temp,p22_temp] = simulate_factors(f_t, m_t, samp, 12,dd,'iter');
                    
                    % Simulate states: s(t+1)...s(t+h)
                    s_temp = simulate_st(st_temp_par,p11_temp,p22_temp,12);
                    
                    % Simulate GDP
                    y_temp = simulate_gdp(s_temp,f_temp,m_temp,samp,dd,'iter');
                    
                    % Average future GDP t+1:t+h
                    y_mat_fut_ex_temp(rcount,jj) = mean(y_temp(1:end,1));
                    
                    % Collect t+h state
                    st_it_temp(rcount,jj)        = s_temp(12);
                    
                    % Collect t state
                    st_it_temp_t(rcount,jj)      = st_sim;
                    
                    
                end
                
                
                rcount = rcount + 1;
            end
            
            
        end
        
        close(wb);
        
        y_mat_fut_ex = y_mat_fut_ex_temp;   %NaN(nDraws*nParamDraws,tperiods); % This stores all realizations for the ITERATED
        st_it        = st_it_temp;          %NaN(nDraws*nParamDraws,tperiods);
        st_it_t      = st_it_temp_t;        %NaN(nDraws*nParamDraws,tperiods);
        
        % Change st=2 (bad) to st = 0 for plotting
        % Mean(st) represents probability of good regime
        st_it(st_it_temp==2) = 0;
        st_it_t(st_it_temp_t==2) = 0;
        
        
        st_mean.it  = mean(st_it)';   % This is in t+12
        st_t_mean.it= mean(st_it_t)'; % This is in t
        dYsim_25.it = prctile(y_mat_fut_ex,25)';
        dYsim_75.it = prctile(y_mat_fut_ex,75)';
        dYsim_10.it = prctile(y_mat_fut_ex,10)';
        dYsim_90.it = prctile(y_mat_fut_ex,90)';
        
        if saveit
            % Plot figure
            fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(1:tperiods, dYsim_10.dgp,'--','Color',colors(5,:),'LineWidth', 1.5,'DisplayName', 'DGP 10th');
            l2=plot(1:tperiods, dYsim_90.dgp,'--','Color',colors(50,:),'LineWidth', 1.5,'DisplayName','DGP 90th');
            l3=plot(1:tperiods, dYsim_10.it,'-','Color',colors(7,:),'LineWidth', 1.5,'DisplayName','Iterated 10th');
            l4=plot(1:tperiods, dYsim_90.it,'-','Color',colors(45,:),'LineWidth', 1.5,'DisplayName','Iterated 90th');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('GDP distribution 1-year ahead (DGP)','FontSize',16','Interpreter','Latex');
            axis tight
            
            subplot(212)
            hold on
            l1=plot(1:tperiods, 1-st_t_mean.dgp(1:tperiods),'--','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','DGP (mean)');
            l2=plot(1:tperiods, 1-st_t_mean.it(1:tperiods),'-','Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Iterated (mean)');
            legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime (single draw)','FontSize',16','Interpreter','Latex');
            axis tight
        end
        
        
        
        
        %------------------------------------------------------------------
        % Direct Draw
        %------------------------------------------------------------------
        
        if exist([sim_folder 'results_direct.mat'],'file')==2
            load([sim_folder 'results_direct.mat'])
            if exist([sim_folder 'results_direct_samp.mat'],'file')==2;load([sim_folder 'results_direct_samp.mat']);end
            
            % Map stored objects
            y_mat = y_fit';
            size_y = size(y_mat,2);
            
            %             ff_temp_par1 = Res_direct.FF(example_draw,:);
            %             mf_temp_par1 = Res_direct.MF(example_draw,:);
            %             st_temp_par1 = st_mat(example_draw,:);
            
            % Pre-allocate output
            y_mat_fut_direct_temp = NaN(nDraws*nParamDraws,tperiods); % This stores all realizations for the direct
            st_direct_temp_t      = NaN(nDraws*nParamDraws,tperiods);
            
            % Waitbar and counter
            wb = waitbar(0,'Calculating paths for the direct version');
            
            % ARe these the direct draws????????
            
            
            
            
            rcount = 1;
            
            % Loop over parameter draws
            for dd=1:nParamDraws
                
                % Counter for rows storing parameter*repetition
                waitbar(dd/nParamDraws,wb,'Calculating paths for the direct version')
                
                % Loop over repetitions
                for uu=1:nDraws
                    
                    % Loop over periods
                    for jj=1:tperiods
                        
                        
                        % Simulate s(t) conditional on realized path of m and f
                        % from t-11 to t.
                        
                        if jj<13
                            st_lag = 1; % Initialize in good regime
                            
                            st_sim = 1; % Remain in good regime
                            
                            % This is s(t)
                            st_temp = st_sim;
                            
                        else
                            
                            % Realized path of f and m from t-11,....,t
                            f_temp = Res_direct.FF(example_draw,jj-11:jj);     %ff_temp_par1(jj-11:jj);
                            m_temp = Res_direct.MF(example_draw,jj-11:jj);     % mf_temp_par1(jj-11:jj);
                            
                            % Compute transition probabilities for s(t-11),...,s(t)
                            p12 = NaN(12,1); p21 = NaN(12,1);
                            
                            for ww=1:12
                                p12(ww) = 1/(1+exp(samp.a12(dd)-samp.b12(dd)*(f_temp(ww))+samp.c12(dd)*(m_temp(ww))));
                                p21(ww) = 1/(1+exp(samp.a21(dd)+samp.b21(dd)*(f_temp(ww))-samp.c21(dd)*(m_temp(ww))));
                            end
                            
                            p11 = ones(12,1) - p12; % probability of remaining in normal
                            p22 = ones(12,1) - p21; % probability of remaining in bad
                            
                            % s(t-11)
                            st_lag = st_direct_temp_t(rcount,jj-12);
                            
                            % This returns s(t-11)...s(t))
                            st_sim = simulate_st(st_lag,p11,p22,12);
                            
                            % This is s(t)
                            st_temp = st_sim(12);
                            
                        end
                        
                        
                        
                        % Macro and financial factor in period t
                        f_t = Res_direct.FF(example_draw,jj);
                        m_t = Res_direct.MF(example_draw,jj);
                        
                        
                        % Simulate GDP
                        y_temp = simulate_gdp(st_temp,f_t,m_t,samp,dd,'direct');
                        
                        % Average future GDP t+1:t+h
                        y_mat_fut_direct_temp(rcount,jj) = y_temp(1);
                        
                        % Collect t state
                        st_direct_temp_t(rcount,jj)      = st_temp;
                        
                        
                    end
                    
                    
                    rcount = rcount + 1;
                end
                
                
            end
            
            close(wb);
            
            y_mat_fut_direct = y_mat_fut_direct_temp;   %NaN(nDraws*nParamDraws,tperiods); % This stores all realizations for the ITERATED
            st_direct_t          = st_direct_temp_t;        %NaN(nDraws*nParamDraws,tperiods);
            
            % Change st=2 (bad) to st = 0 for plotting
            % Mean(st) represents probability of good regime
            st_direct_t(st_direct_temp_t==2) = 0;
            
            
            
            st_t_mean.direct= mean(st_direct_t)'; % This is in t
            dYsim_25.direct = prctile(y_mat_fut_direct,25)';
            dYsim_75.direct = prctile(y_mat_fut_direct,75)';
            dYsim_10.direct = prctile(y_mat_fut_direct,10)';
            dYsim_90.direct = prctile(y_mat_fut_direct,90)';
            
            tperiods_graph = 312;
            %             tperiods_graph = tperiods;
            
            % Plot figure
            fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(13:tperiods_graph, dYsim_10.dgp(13:tperiods_graph),'--','Color',colors(5,:),'LineWidth', 2,'DisplayName', 'DGP 10th');
            l2=plot(13:tperiods_graph, dYsim_90.dgp(13:tperiods_graph),'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','DGP 90th');
            l3=plot(13:tperiods_graph, dYsim_10.it(13:tperiods_graph),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated 10th');
            l4=plot(13:tperiods_graph, dYsim_90.it(13:tperiods_graph),'-','Color',colors(45,:),'LineWidth', 2,'DisplayName','Iterated 90th');
            l5=plot(13:tperiods_graph, dYsim_10.direct(13:tperiods_graph),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct 10th');
            l6=plot(13:tperiods_graph, dYsim_90.direct(13:tperiods_graph),':','Color',colors(50,:),'LineWidth', 2.5,'DisplayName','Direct 90th');
            legend([l1 l2 l3 l4 l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
            axis tight
            
            subplot(212)
            hold on
            l1=plot(13:tperiods_graph, 1-st_t_mean.dgp(13:tperiods_graph),'--','Color',colors(5,:),'LineWidth', 2,'DisplayName','DGP (mean)');
            l2=plot(13:tperiods_graph, 1-st_t_mean.it(13:tperiods_graph),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated (mean)');
            l3=plot(13:tperiods_graph, 1-st_t_mean.direct(13:tperiods_graph),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct (mean)');
            legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
            axis tight
            ylim([0 1])
            
            tightfig
            
            if saveit ==1
                print('-dpdf',fig,[sim_folder 'pdens_par_unc_new'],'-bestfit');
                saveas(fig,sprintf('%s.png',[sim_folder 'pdens_par_unc_new']));
            end
            
             % Plot figure
            fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            plot(13:tperiods_graph, zeros(size(13:tperiods_graph,2),1),':k','LineWidth', 1.5);
%             l3=plot(13:tperiods_graph, dYsim_10.it(13:tperiods_graph)-dYsim_10.dgp(13:tperiods_graph),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated 10th');
%             l4=plot(13:tperiods_graph, dYsim_90.it(13:tperiods_graph)-dYsim_90.dgp(13:tperiods_graph),'-','Color',colors(45,:),'LineWidth', 2,'DisplayName','Iterated 90th');
            l5=plot(13:tperiods_graph, dYsim_10.direct(13:tperiods_graph)-dYsim_10.dgp(13:tperiods_graph),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct 10th');
            l6=plot(13:tperiods_graph, dYsim_90.direct(13:tperiods_graph)-dYsim_90.dgp(13:tperiods_graph),':','Color',colors(50,:),'LineWidth', 2.5,'DisplayName','Direct 90th');
            legend([l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Discrepancy of GDP distribution 1-year ahead (Direct - d.g.p.)','FontSize',16','Interpreter','Latex');
            axis tight
            
            subplot(212)
            hold on
            plot(13:tperiods_graph, zeros(size(13:tperiods_graph,2),1),':k','LineWidth', 1.5);
%             l2=plot(13:tperiods_graph, 1-st_t_mean.it(13:tperiods_graph)-(1-st_t_mean.dgp(13:tperiods_graph)),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated (mean)');
            l3=plot(13:tperiods_graph, 1-st_t_mean.direct(13:tperiods_graph)-(1-st_t_mean.dgp(13:tperiods_graph)),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct (mean)');
            legend(l3,'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Discrepancy of the Probability of Bad Regime (Direct - d.g.p.)','FontSize',16','Interpreter','Latex');
            axis tight
            
            if saveit ==1
                print('-dpdf',fig,[sim_folder 'pdens_discrepancy'],'-bestfit');
                saveas(fig,sprintf('%s.png',[sim_folder 'pdens_discrepancy']));
            end
        end
    end
end




toc;
