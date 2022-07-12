% simulate_DirectFull.m 
%
% Simulates the measurement equation for the GDP equation conditional on a
% path for m(t), f(t), s(t) and shocks eps(t). 
%
% y(t) = c(s(t)) + alpha(s(t))*m(t) + beta*(s(t))*f(t) + gamma(s(t))*eps(t)
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2] = simulate_DirectFull(st,Ymat,param,opt)

% Number of periods
nsim = length(st);

% Matrix of shocks
ETA = randn(opt.nvars,nsim); 

%% Collect reduced-form matrices
D_sync_1 = param.D_sync_1;
D_sync_2 = param.D_sync_2;
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;

%% Simulate variables based on simulated path for st and factors m(tt) and f(tt)
Ynext_1 = D_sync_1 + B_sync_1*Ymat + O_sync_1*ETA;
Ynext_2 = D_sync_2 + B_sync_2*Ymat + O_sync_2*ETA;

if st==2
    % Bad regime
    Ynext = Ynext_2;
else
    % Good regime
    Ynext = Ynext_1;
end
        
%% Collect the simulation
y_out = Ynext;
y_out_1 = Ynext_1;
y_out_2 = Ynext_2;

