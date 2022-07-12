%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous

%% housekeeping
clear; close all; clc; tic;

cd('../')

% Important paths
addpath('../RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
dir_or_it    = 2; % 1 = direct, 2 = iterated
hh           = 12; % forecast horizon
const        = 1; % 1 = have a constant in transition probability
normal       = 0; % 1 = use normal distribution, 0 = gamma distribution
saveit       = 1; % 1 = save figures
num_rep      = 5000; % Number of repetitions for the coverage bands of the DGP
nParamDraws  = 10000; %samp.options.N;    % Number of parameter draws
nDraws       = 10;   % Number of simulations for each parameter draw


data_type    = 1; % 1, for LHS our GDP, RHS IP and EBP
% 2, for LHS Markit GDP, RHS MF and FF


% Data vintage, sample and country selection
datafilename = 'Data_exercise.xlsx';
sheetuse     = 'DFM_73_Monthly';
start_date   = '1973-Jan';
end_date     = '2020-Oct';


% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];

% Create date formats for plotting
inputformat = 'yyyy-MMM';
dataformat  = 'yyyy-mmm';
start_plot    = '1973-Feb';
end_plot      = '2019-Oct';
end_plot_full = '2020-Oct';

if data_type ==1
    ex_data = 'RHS';
elseif data_type ==2
    ex_data = 'LHS';
end

if dir_or_it ==1
    model = 'results_direct';
elseif dir_or_it ==2
    model = 'results_iterated';
end

sim_folder = 'Results/Data_exercise/';
if exist(sim_folder,'dir')==0
    mkdir(sim_folder)
end


% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Vector of dates for plotting
dates   = dates_full(sd:ed);
tperiods     = length(dates);  % Number of time-periods in simulation of GDP

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
if exist([sim_folder '' ex_data '_' model '.mat'],'file')==0
    run estimation_data.m; % This only estimates and stores the posterior mode
    save([sim_folder '' ex_data '_' model '.mat'])
end

if exist([sim_folder '' ex_data '_' model '_samp.mat'],'file')==0
    samp=posterior_sampling_data(dir_or_it,nlags,const,normal,db);
    save([sim_folder '' ex_data '_' model '_samp.mat'],'samp')
end



%% PARAMETER UNCERTAINTY
close all;

if exist([sim_folder '' ex_data '_results_iterated.mat'],'file')==2
    load([sim_folder '' ex_data '_results_iterated.mat'])
    if exist([sim_folder '' ex_data '_results_iterated_samp.mat'],'file')==2
        load([sim_folder '' ex_data '_results_iterated_samp.mat']);
    end
    
    trend_fit = db.TRENDH.data(nlags+1:end);
    
    %%
    %------------------------------------------------------------------
    % Iterated Draw: Estimated Parameter
    %------------------------------------------------------------------
    
    % Map stored objects
    st_mat = st_fit'; % Now use the fitted values
    y_mat = y_fit'; % Now use the fitted values
    
    % Pre-allocate
    y_mat_fut_ex_temp = NaN(nParamDraws,nDraws,tperiods); % This stores all realizations for the ITERATED
    st_it_temp        = NaN(nParamDraws,nDraws,tperiods);
    st_it_temp_t      = NaN(nParamDraws,nDraws,tperiods);
    
    
    % Loop over parameter draws
    parfor dd=1:nParamDraws
        
        % Loop over repetitions
        for uu=1:nDraws
            
            % Loop over periods
            for jj=1:tperiods
                
                % Macro and financial factor in period t
                f_t = FF(jj);
                m_t = MF(jj);
                
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
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
                
                % This is s(t)
                st_temp_par = st_sim;
                
                % Simulate factors and transtion probabilities
                % conditional on f_t and m_t from t+1 to t+h
                [f_temp,m_temp,p11_temp,p22_temp] = simulate_factors(f_t, m_t, samp, hh,dd,'iter');
                
                % Simulate states: s(t+1)...s(t+h)
                s_temp = simulate_st(st_temp_par,p11_temp,p22_temp,hh);
                
                % Simulate GDP
                y_temp = simulate_gdp(s_temp,f_temp,m_temp,samp,dd,'iter');
                
                % Average future GDP t+1:t+h
                y_mat_fut_ex_temp(dd,uu,jj) = mean(y_temp(1:end,1))+trend_fit(jj);
                
                % Collect t+h state
                st_it_temp(dd,uu,jj)        = s_temp(hh);
                
                % Collect t state
                st_it_temp_t(dd,uu,jj)      = st_sim;
                
                % This will be the st_lag for the next jj period iteration
                st_lag = st_sim;
                
            end
        end
    end
        
    % Reshape matrices
    y_mat_fut_ex = reshape(y_mat_fut_ex_temp,[nDraws*nParamDraws,tperiods]);
    st_it = reshape(st_it_temp,[nDraws*nParamDraws,tperiods]);
    st_it_t = reshape(st_it_temp_t,[nDraws*nParamDraws,tperiods]);
    
    % Change st=2 (bad) to st = 0 for plotting
    % Mean(st) represents probability of good regime
    st_it(st_it==2) = 0;
    st_it_t(st_it_t==2) = 0;
    
    st_mean.it  = mean(st_it)';   % This is in t+hh
    st_t_mean.it= mean(st_it_t)'; % This is in t
    dYsim_25.it = prctile(y_mat_fut_ex,25)';
    dYsim_75.it = prctile(y_mat_fut_ex,75)';
    dYsim_10.it = prctile(y_mat_fut_ex,10)';
    dYsim_90.it = prctile(y_mat_fut_ex,90)';
    
    
    %%
    
    %------------------------------------------------------------------
    % Direct Draw
    %------------------------------------------------------------------
    
    if exist([sim_folder '' ex_data '_results_direct.mat'],'file')==2
        load([sim_folder '' ex_data '_results_direct.mat'])
        if exist([sim_folder '' ex_data '_results_direct_samp.mat'],'file')==2
            load([sim_folder '' ex_data '_results_direct_samp.mat']);
        end
        
        % Map stored objects
        y_mat = y_fit';
        size_y = size(y_mat,2);
        
        % Pre-allocate output
        y_mat_fut_direct_temp = NaN(nParamDraws,nDraws,tperiods); % This stores all realizations for the direct
        st_direct_temp_t      = NaN(nParamDraws,nDraws,tperiods);
        st_lag_temp           = NaN(nParamDraws,nDraws,tperiods);
        
        % Loop over parameter draws
        parfor dd=1:nParamDraws
            
            % Counter for rows storing parameter*repetition
            %             waitbar(dd/nParamDraws,wb,'Calculating paths for the direct version')
            
            % Loop over repetitions
            for uu=1:nDraws
                st_sim = 1; % initialization of the simulated states
                
                % Loop over periods
                for jj=1:tperiods
                    
                    % Simulate s(t) conditional on realized path of m and f
                    % from t-11 to t.
                    if jj<hh+1
                        
                        % This is s(t)
                        st_temp = 1;
                        
                    else
                        
                        % Realized path of f and m from t-11,....,t
                        f_temp = FF(jj-11:jj);     %ff_temp_par1(jj-11:jj);
                        m_temp = MF(jj-11:jj);     % mf_temp_par1(jj-11:jj);
                        
                        % Compute transition probabilities for s(t-11),...,s(t)
                        p12 = NaN(hh,1); p21 = NaN(hh,1);
                        
                        for ww=1:hh
                            p12(ww) = 1/(1+exp(samp.a12(dd)-samp.b12(dd)*(f_temp(ww))+samp.c12(dd)*(m_temp(ww))));
                            p21(ww) = 1/(1+exp(samp.a21(dd)+samp.b21(dd)*(f_temp(ww))-samp.c21(dd)*(m_temp(ww))));
                        end
                        
                        p11 = ones(hh,1) - p12; % probability of remaining in normal
                        p22 = ones(hh,1) - p21; % probability of remaining in bad
                        
                        % s(t-11)
                        st_lag = st_sim(1);
                        
                        % This returns s(t-11)...s(t))
                        st_sim = simulate_st(st_lag,p11,p22,hh);
                        
                        % This is s(t)
                        st_temp = st_sim(hh);
                        
                    end
                    
                    % Macro and financial factor in period t
                    f_t = FF(jj);
                    m_t = MF(jj);
                    
                    % Simulate GDP
                    y_temp = simulate_gdp(st_temp,f_t,m_t,samp,dd,'direct')+trend_fit(jj);
                    
                    % Average future GDP t+1:t+h
                    y_mat_fut_direct_temp(dd,uu,jj) = y_temp(1);
                    
                    % Collect t state
                    st_direct_temp_t(dd,uu,jj)      = st_temp;
                    
                end
            end
        end
                
        % Reshape matrices
        y_mat_fut_direct = reshape(y_mat_fut_direct_temp,[nDraws*nParamDraws,tperiods]);
        st_direct_t = reshape(st_direct_temp_t,[nDraws*nParamDraws,tperiods]);
        
        % Change st=2 (bad) to st = 0 for plotting
        % Mean(st) represents probability of good regime
        st_direct_t(st_direct_t==2) = 0;
        
        st_t_mean.direct= mean(st_direct_t)'; % This is in t
        dYsim_25.direct = prctile(y_mat_fut_direct,25)';
        dYsim_75.direct = prctile(y_mat_fut_direct,75)';
        dYsim_10.direct = prctile(y_mat_fut_direct,10)';
        dYsim_90.direct = prctile(y_mat_fut_direct,90)';
        
        tperiods_graph = tperiods;
        
        % Plot figure
        fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
        subplot(211)
        hold on
        l3=plot(dates_full(sd:ed), dYsim_10.it,'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated 10th');
        l4=plot(dates_full(sd:ed), dYsim_90.it,'-','Color',colors(45,:),'LineWidth', 2,'DisplayName','Iterated 90th');
        l5=plot(dates_full(sd:ed), dYsim_10.direct,':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct 10th');
        l6=plot(dates_full(sd:ed), dYsim_90.direct,':','Color',colors(50,:),'LineWidth', 2.5,'DisplayName','Direct 90th');
        legend([l3 l4 l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
        set(gca,'children',flipud(get(gca,'children')))
        hold off
        ylabel('Percent','interpreter','Latex','fontsize',10)
        title('GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
        axis tight
        datetick('x','yyyy','keepticks')
        set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
        %
        subplot(212)
        hold on
        l2=plot(dates_full(sd:ed), 1-st_t_mean.it,'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated (mean)');
        l3=plot(dates_full(sd:ed), 1-st_t_mean.direct,':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct (mean)');
        legend([l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
        set(gca,'children',flipud(get(gca,'children')))
        hold off
        ylabel('Percent','interpreter','Latex','fontsize',10)
        title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
        axis tight
        datetick('x','yyyy','keepticks')
        set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
        ylim([0 1])
        tightfig;
        if saveit ==1
            print('-dpdf',fig,[sim_folder 'pdens_par_unc_new_data'],'-bestfit');
            saveas(fig,sprintf('%s.png',[sim_folder 'pdens_par_unc_new_data']));
        end
        
        %{
        % Plot figure
        fig=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
        subplot(211)
        hold on
        plot((hh+1):tperiods_graph, zeros(size((hh+1):tperiods_graph,2),1),':k','LineWidth', 1.5);
%             l3=plot((hh+1):tperiods_graph, dYsim_10.it((hh+1):tperiods_graph)-dYsim_10.dgp((hh+1):tperiods_graph),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated 10th');
%             l4=plot((hh+1):tperiods_graph, dYsim_90.it((hh+1):tperiods_graph)-dYsim_90.dgp((hh+1):tperiods_graph),'-','Color',colors(45,:),'LineWidth', 2,'DisplayName','Iterated 90th');
        l5=plot((hh+1):tperiods_graph, dYsim_10.direct((hh+1):tperiods_graph)-dYsim_10.dgp((hh+1):tperiods_graph),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct 10th');
        l6=plot((hh+1):tperiods_graph, dYsim_90.direct((hh+1):tperiods_graph)-dYsim_90.dgp((hh+1):tperiods_graph),':','Color',colors(50,:),'LineWidth', 2.5,'DisplayName','Direct 90th');
        legend([l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
        set(gca,'children',flipud(get(gca,'children')))
        hold off
        ylabel('Percent','interpreter','Latex','fontsize',10)
        title('Discrepancy of GDP distribution 1-year ahead (Direct - d.g.p.)','FontSize',16','Interpreter','Latex');
        axis tight

        subplot(212)
        hold on
        plot((hh+1):tperiods_graph, zeros(size((hh+1):tperiods_graph,2),1),':k','LineWidth', 1.5);
%             l2=plot((hh+1):tperiods_graph, 1-st_t_mean.it((hh+1):tperiods_graph)-(1-st_t_mean.dgp((hh+1):tperiods_graph)),'-','Color',colors(7,:),'LineWidth', 2,'DisplayName','Iterated (mean)');
        l3=plot((hh+1):tperiods_graph, 1-st_t_mean.direct((hh+1):tperiods_graph)-(1-st_t_mean.dgp((hh+1):tperiods_graph)),':','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Direct (mean)');
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
        %}
    end
end


toc;
