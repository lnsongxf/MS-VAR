% simulate_factors.m
%
% Simulate factors for nsim periods starting from m_lag and f_lag
%
%--------------------------------------------------------------------------

function [f_out, m_out,p11, p22] = simulate_factors(f_lag, m_lag, param,nsim,dd,sign_a0)

% Allocate objects
m_out = NaN(nsim,1);
f_out = NaN(nsim,1);

% Map parameters
c_2_1 = param.c_2_1(dd);
c_1_1 = param.c_1_1(dd);
a1_2_1= param.a1_2_1(dd);
a1_1_1= param.a1_1_1(dd);
a1_2_2= param.a1_2_2(dd);
a1_1_2= param.a1_1_2(dd);
s_2_2 = param.s_2_2(dd);
s_1_1 = param.s_1_1(dd);
a12   = param.a12(dd);
b12   = param.b12(dd);
c12   = param.c12(dd);
a21   = param.a21(dd);
b21   = param.b21(dd);
c21   = param.c21(dd);

if strcmp(sign_a0,'dgp')
    a0_1_2= param.a0_1_2(dd);
else
    a0_1_2= -param.a0_1_2(dd);
end

% financial and macro shocks
eta1 = randn(2,nsim); 

for tt=1:nsim
    
    % Simulate factors
    m_out(tt) = c_2_1 +                     a1_2_1*f_lag + a1_2_2*m_lag + s_2_2*eta1(2,tt);
    f_out(tt) = c_1_1 + a0_1_2*m_out(tt) + a1_1_1*f_lag + a1_1_2*m_lag + s_1_1*eta1(1,tt);
    
    % update state
    m_lag = m_out(tt);
    f_lag = f_out(tt);
end

% Transition probabilities
p12 = 1./(1+exp(a12-b12*(f_out)+c12*(m_out)));
p21 = 1./(1+exp(a21+b21*(f_out)-c21*(m_out)));

% State probabilities t+1,...,t+h
p11 = ones(12,1) - p12; % probability of remaining in normal
p22 = ones(12,1) - p21; % probability of remaining in bad
