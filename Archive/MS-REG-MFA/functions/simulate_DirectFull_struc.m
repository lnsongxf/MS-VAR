% simulate_DirectFull.m 
%
% Simulates the measurement equation for the GDP equation conditional on a
% path for m(t), f(t), s(t) and shocks eps(t). 
%
% y(t) = c(s(t)) + alpha(s(t))*m(t) + beta*(s(t))*f(t) + gamma(s(t))*eps(t)
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2] = simulate_DirectFull_struc(st,f_temp,m_temp,gdp_temp,f_temp_lag,m_temp_lag,gdp_temp_lag,param)

% Number of periods
nsim = length(st);

% GDP shock
eta2 = randn(3,nsim); 

% Store GDP simulations
Ynext   = NaN(3,nsim);

%% Collect reduced-form matrices
C_sync_1 = param.C_sync_1;
C_sync_2 = param.C_sync_2;
A0_sync_1 = -param.A0_sync_1;
A0_sync_2 = -param.A0_sync_2;
A1_sync_1 = param.A1_sync_1;
A1_sync_2 = param.A1_sync_2;
SIG_sync_1 = param.SIG_sync_1;
SIG_sync_2 = param.SIG_sync_2;

Ymat = [f_temp;m_temp;gdp_temp];
Ymat_lag = [f_temp_lag;m_temp_lag;gdp_temp_lag];

% this only creates valid forecasts for
y_out_1 = C_sync_1(end,:) + A0_sync_1(end,1:2)*Ymat(1:2) +A1_sync_1(end,:)*Ymat_lag + SIG_sync_1(end,:)*eta2;
y_out_2 = C_sync_2(end,:) + A0_sync_2(end,1:2)*Ymat(1:2) +A1_sync_2(end,:)*Ymat_lag + SIG_sync_2(end,:)*eta2;

% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)    
if st==2
    % Bad regime
    Ynext = y_out_2;
else
    % Good regime
    Ynext= y_out_1;
end

y_out = Ynext;
