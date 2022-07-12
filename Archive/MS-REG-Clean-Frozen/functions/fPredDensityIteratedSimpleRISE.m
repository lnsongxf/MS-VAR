%------------------------------------------------------------------
% Predictive Density - Iterated simulation using RISE forecast function
%
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIteratedSimpleRISE(trend_fit,p_reg2_filtered,param_in,FF, MF, nParamDraws,nDraws,tperiods,hh,dates_index,temp,dates_full,sv)


% Pre-allocate
y_mat_fut_ex_temp = NaN(nParamDraws,nDraws,tperiods-12); 
stemp_hh          = NaN(nParamDraws,nDraws,tperiods-12);
stemp_t           = NaN(nParamDraws,nDraws,tperiods-12);


% Set forecast horizon
% Need to forecast at least 2 periods otherwise RISE crashes
nsteps=hh+1;
shock_uncertainty = true;

% Loop over parameter draws
for dd=1:nParamDraws
    
    % Loop over repetitions
    for uu=1:nDraws
        % Loop over periods
        % *** Only works over estimation sample that ends 
        % *** 12 months before last observation
        fprintf('Current draw: %i \n',uu);
        
        parfor tt=1:(tperiods-12)

            % Construct Date of first period in forecast                
            yeart = year(dates_full(tt+1));
            montht= month(dates_full(tt+1));            
            start_period = [num2str(yeart) 'm' num2str(montht)];

            % Forecast     
            mycast=forecast(sv,[],start_period,[],nsteps,shock_uncertainty);
              
            % Collect GDPG simulation 
            % Discard first observation that corresponds to initial conditions
            y_temp = mycast.GDPG.data(2:end);

            % Average future GDP t+1:t+h
            y_mat_fut_ex_temp(dd,uu,tt) = mean(y_temp(1:hh,1)) + trend_fit(tt);

            % % Collect state in period t + hh
            % stemp_hh(dd,uu,tt) = s_temp(hh);
            
            % % Collect state in period t
            % stemp_t(dd,uu,tt) = st_sim;
                      
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex  = reshape(y_mat_fut_ex_temp,[nDraws*nParamDraws,tperiods-12]);


% stmat_hh     = reshape(stemp_hh,[nDraws  *nParamDraws,tperiods]);
% stmat_t      = reshape(stemp_t,[nDraws*  nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
% stmat_hh(stmat_hh==2) = 0;
% stmat_t(stmat_t==2) = 0;

% out.st_h_mean = mean(stmat_hh)';    % This is in t+hh
% out.st_t_mean = mean(stmat_t)';     % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean     = mean(y_mat_fut_ex)';

% This is the predictive density matrix (nDraws*nParamDraws)xselected_periods
if ~isempty(dates_index)
    pden_mat = y_mat_fut_ex(:,dates_index);
else
    pden_mat = [];
end
