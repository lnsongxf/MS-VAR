% simulate_gdp.m 
%
% Simulates the measurement equation for the GDP equation conditional on a
% path for m(t), f(t), s(t) and shocks eps(t). 
%
% y(t) = c(s(t)) + alpha(s(t))*m(t) + beta*(s(t))*f(t) + gamma(s(t))*eps(t)
%--------------------------------------------------------------------------

function y_out = simulate_gdp(st,f_temp,m_temp,param,dd,sign_a0)

% Number of periods
nsim = length(st);

% GDP shock
eta2 = randn(1,nsim); 

% Store GDP simulations
y_out = NaN(nsim,1);


% Map parameters
c_3_1_sync_1 = param.c_3_1_sync_1(dd);
c_3_1_sync_2 = param.c_3_1_sync_2(dd);
s_3_3_sync_1 = param.s_3_3_sync_1(dd);
s_3_3_sync_2 = param.s_3_3_sync_2(dd);

if strcmp(sign_a0,'dgp')
    a0_3_1_sync_1 = param.a0_3_1_sync_1(dd);
    a0_3_1_sync_2 = param.a0_3_1_sync_2(dd);
    a0_3_2_sync_1 = param.a0_3_2_sync_1(dd);
    a0_3_2_sync_2 = param.a0_3_2_sync_2(dd);
else
    a0_3_1_sync_1 = -param.a0_3_1_sync_1(dd);
    a0_3_1_sync_2 = -param.a0_3_1_sync_2(dd);
    a0_3_2_sync_1 = -param.a0_3_2_sync_1(dd);
    a0_3_2_sync_2 = -param.a0_3_2_sync_2(dd);
end


% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
for tt=1:nsim 
    
    if st(tt)==2
        % Bad regime
        y_out(tt) = c_3_1_sync_2+a0_3_1_sync_2*f_temp(tt) +a0_3_2_sync_2*m_temp(tt) + s_3_3_sync_2*eta2(1,tt);
    else
        % Good regime
        y_out(tt) = c_3_1_sync_1+a0_3_1_sync_1*f_temp(tt) + a0_3_2_sync_1*m_temp(tt) + s_3_3_sync_1*eta2(1,tt);
    end
        
end
