%------------------------------------------------------------------
% Predictive Density - Direct simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out,pden_mat] = fPredDensityDirectFull_RISE(sv,param_in,FF,MF,GDP,trend_fit,opt,temp,dates_full,db_full)

tperiods = opt.tperiods;
nlags = sv.nlags;

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(opt.nParamDraws,opt.nDraws,tperiods); 
y_mat_fut_ex_temp_good = NaN(opt.nParamDraws,opt.nDraws,tperiods); 
y_mat_fut_ex_temp_bad = NaN(opt.nParamDraws,opt.nDraws,tperiods); 

% Pre-allocate state matrices
stemp_t   = NaN(opt.nParamDraws,opt.nDraws,tperiods);

% Collect filtered probabilities
[~,~,~,f]=filter(sv);
p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;

% Loop over parameter draws
for dd=1:opt.nParamDraws
    
    param = scriptParams_FULL_companion(sv,param_in,dd); % works for any lag and model
    
    pbad = p_reg2_filtered;
    
    
    % Loop over repetitions
    parfor uu=1:opt.nDraws
%     for uu=1:opt.nDraws
        % Loop over periods
        if mod(uu,1000)==0
        fprintf('Current draw: %i of %i\n',uu,opt.nDraws);
        end
        
        st_sim = NaN(tperiods,1);
        
        % Loop over periods
        for tt=1:tperiods
            
            if opt.simul==0

                % Simulate s(t) conditional on realized path of m and f
                % from t-opt.hh+1 to t.
                if tt<=opt.hh
                    % In the first hh we can construct s(t) from the filtered probabilities

                    % Use filtered probability to determine current state
                    udraw = rand(1); % draw a coin that determines the regime
                    if udraw < (1-pbad(tt))
                        st_temp = 1;
                    else
                        st_temp = 2;
                    end
                else

                    % Get current state
                    % Here we have two options:
                    % drawst = draw from filtered states every period
                    % forecastst = recursively simulate st

                    if strcmp(temp,'drawst')
                        % Use filtered probability to determine current state
                        udraw = rand(1); % draw a coin that determines the regime
                        if udraw < (1-pbad(tt))
                            st_temp = 1;
                        else
                            st_temp = 2; % bad regime
                        end

                    elseif strcmp(temp,'forecastst')
                        % s(t-opt.hh)
                        st_lag = st_sim(tt-opt.hh);

                        % Realized path of f and m from t-opt.hh+1 to t
                        f_temp = FF(tt-opt.hh+1:tt);
                        m_temp = MF(tt-opt.hh+1:tt);

                        % Compute transition probabilities for s(t-opt.hh+1),...,s(t)
                        [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);

                        p11 = ones(opt.hh,1) - p12; % probability of remaining in normal
                        p22 = ones(opt.hh,1) - p21; % probability of remaining in bad

                        % This returns s(t-opt.hh+1)...s(t)
                        st_sim_temp = simulate_st(st_lag,p11,p22,opt.hh);

                        % This is s(t)
                        st_temp = st_sim_temp(opt.hh);
                    end
                end

                st_sim(tt) = st_temp;

                % Macro and financial factor in period t
                f_t = FF(tt);
                m_t = MF(tt);
                gdp_t = GDP(tt);

                % Simulate GDP
                [y_temp,y_temp1,y_temp2] = simulate_DirectFull(st_temp,f_t,m_t,gdp_t,param);

                % Average future GDP t+1:t+h
                y_mat_fut_ex_temp(dd,uu,tt) = y_temp+trend_fit(tt);

                y_mat_fut_ex_temp_good(dd,uu,tt) = y_temp1 + trend_fit(tt);
                y_mat_fut_ex_temp_bad(dd,uu,tt) = y_temp2+ trend_fit(tt);

                % Collect t state
                stemp_t(dd,uu,tt) = st_temp;
                
            else
                
                %___________________________________________________________________%
                %-------------------------------------------------------------------%
                % RISE OPTION %

                % forecast function in RISE does not work for first period
                % since it cuts the data
                if tt==nlags
                    
                    udraw = rand(1); % draw a coin that determines the regime

                    if udraw < (1-pbad(tt))
                        st_temp = 1;
                    else
                        st_temp = 2; % bad regime
                    end
                    
                    % Macro and financial factor in period t
                    f_t = FF(tt);
                    m_t = MF(tt);
                    gdp_t = GDP(tt);

                    % Simulate GDP
                    [y_temp,y_temp1,y_temp2] = simulate_DirectFull(st_temp,f_t,m_t,gdp_t,param);
                    
                else
                
                    % Construct Date of first period in forecast    
                    yeart = year(dates_full(tt));
                    montht= month(dates_full(tt));            
                    start_period = [num2str(yeart) 'm' num2str(montht)];

                    % Forecast     
                    shock_uncertainty = 'true';
                    mycast=forecast(sv,db_full,start_period,[],opt.hh,shock_uncertainty);

                    % Collect GDPG simulation 
                    % Discard first observation that corresponds to initial conditions
                    y_temp = mycast.GDPGH.data;

                    % Average future GDP t+1:t+h
                    y_mat_fut_ex_temp(dd,uu,tt) = mean(y_temp(1:opt.hh,1)) + trend_fit(tt);
                    y_mat_fut_ex_temp_good(dd,uu,tt) = y_mat_fut_ex_temp(dd,uu,tt);
                    y_mat_fut_ex_temp_bad(dd,uu,tt) = y_mat_fut_ex_temp(dd,uu,tt);

                    % Collect state in period t
                    stemp_t(dd,uu,tt) = 1;

                end
                %___________________________________________________________________%
                %-------------------------------------------------------------------%
                
            end
            
        end
    end
end

% Reshape matrices
y_mat_fut_ex  = reshape(y_mat_fut_ex_temp,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);
y_mat_fut_ex_bad  = reshape(y_mat_fut_ex_temp_bad,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);


stmat_t        = reshape(stemp_t,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);


% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;


out.st_t_mean= mean(stmat_t)'; % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean      = mean(y_mat_fut_ex)';

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xselected_periods
if ~isempty(opt.date_index)
    pden_mat.full     = y_mat_fut_ex;
    pden_mat.realized = y_mat_fut_ex(:,opt.date_index);
    pden_mat.good     = y_mat_fut_ex_good(:,opt.date_index);
    pden_mat.bad      = y_mat_fut_ex_bad(:,opt.date_index);
else
    pden_mat = [];
end
