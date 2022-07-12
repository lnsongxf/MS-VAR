% simulate_DirectFull.m 
%
% Simulates the measurement equation for the GDP equation conditional on a
% path for m(t), f(t), s(t) and shocks eps(t). 
%
% y(t) = c(s(t)) + alpha(s(t))*m(t) + beta*(s(t))*f(t) + gamma(s(t))*eps(t)
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2] = simulate_DirectFull(st,f_temp,m_temp,gdp_temp,param)

% Number of periods
nsim = length(st);

% GDP shock
eta2 = randn(3,nsim); 

% Store GDP simulations
Ynext   = NaN(3,nsim);

%% Collect reduced-form matrices
D_sync_1 = param.D_sync_1;
D_sync_2 = param.D_sync_2;
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;

Ymat = [f_temp;m_temp;gdp_temp];

Ynext_1 = D_sync_1 + B_sync_1*Ymat + O_sync_1*eta2;
Ynext_2 = D_sync_2 + B_sync_2*Ymat + O_sync_2*eta2;

y_out_1 = Ynext_1(end,:);
y_out_2 = Ynext_2(end,:);

% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
    
if st==2
    % Bad regime
    Ynext = Ynext_2;
else
    % Good regime
    Ynext = Ynext_1;
end
        

y_out = Ynext(end,:);

