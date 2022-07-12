% simulate_gdpSimple.m 
%
% Simulates the measurement equation for the GDP equation conditional on a
% path for m(t), f(t), s(t) and shocks eps(t). 
%
% y(t) = c(s(t)) + alpha(s(t))*m(t) + beta*(s(t))*f(t) + gamma(s(t))*eps(t)
%--------------------------------------------------------------------------

function y_out = simulate_gdpSimple(st,f_temp,m_temp,param,dd,sign_a0)

% Number of periods
nsim = length(st);

% GDP shock
eta2 = randn(1,nsim); 

% Store GDP simulations
y_out = NaN(nsim,1);


% Map parameters
c_2_1_sync_1 = param.c_2_1_sync_1(dd);
c_2_1_sync_2 = param.c_2_1_sync_2(dd);
s_2_2_sync_1 = param.s_2_2_sync_1(dd);
s_2_2_sync_2 = param.s_2_2_sync_2(dd);

% Impact coefficients
a0_2_1_sync_1 = -param.a0_2_1_sync_1(dd);
a0_2_1_sync_2 = -param.a0_2_1_sync_2(dd);
a0_2_2_sync_1 = -param.a0_2_2_sync_1(dd);
a0_2_2_sync_2 = -param.a0_2_2_sync_2(dd);


% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
for tt=1:nsim 
    
    if st(tt)==2
        % Bad regime
        y_out(tt) = c_2_1_sync_2 + a0_2_1_sync_2*f_temp(tt) + a0_2_2_sync_2*m_temp(tt) + s_2_2_sync_2*eta2(1,tt);
    else
        % Good regime
        y_out(tt) = c_2_1_sync_1 + a0_2_1_sync_1*f_temp(tt) + a0_2_2_sync_1*m_temp(tt) + s_2_2_sync_1*eta2(1,tt);
    end
        
end
