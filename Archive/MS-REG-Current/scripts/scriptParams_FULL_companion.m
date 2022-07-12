function [params] = scriptParams_FULL_companion(sv,params_in,dd)

% Y_t = B Y_t-1 + B*u_t
F = fieldnames(params_in);
C = num2cell(struct2cell(params_in),3);
C = cellfun(@(c)[c{:}],C,'uni',0);
params_aux = vertcat(C{:});
myparams = params_aux(:,dd);
% sv.estim_.estim_param=params_in(:,dd); % -cannot set this read-only
% property, but can push it within the "solve" command below
sol=solve(sv,myparams);
nx = sol.nx;
nvars = sol.nvars;

% VAR autoregressive matrix
solB = sol(1).B;
B=solB(:,nx+1:end,:);
B_sync_1 = B(:,:,1);
B_sync_2 = B(:,:,2);

% VCV matrix of structural coefficients
SIG_sync_1 = diag(sol(1).S0(:,1));
SIG_sync_2 = diag(sol(1).S0(:,2));

C_sync_1  = sol.A(:,1:nx,1);
C_sync_2  = sol.A(:,1:nx,2);

A0_sync_1  = sol.A(:,nx+(1:nvars),1);
A0_sync_2  = sol.A(:,nx+(1:nvars),2);

A1_sync_1  = sol.A(:,(nx+(1:nvars)+nvars):end,1);
A1_sync_2  = sol.A(:,(nx+(1:nvars)+nvars):end,2);
% A1_sync_1(3,3) = 0;
% A1_sync_2(3,3) = 0;

% 
A0_sync_1(2,1) = 0;
A0_sync_2(2,1) = 0;
% 
% A1_sync_1(1,2) = 0;
% A1_sync_1(1,3) = 0;
% 
% A1_sync_2(1,2) = 0;
% A1_sync_2(1,3) = 0;
% 
% A1_sync_1(2,1) = 0;
% A1_sync_1(2,3) = 0;
% 
% A1_sync_2(2,1) = 0;
% A1_sync_2(2,3) = 0;


% A1_sync_1(1,3) = 0;
% A1_sync_1(2,3) = 0;
% 
% A1_sync_2(1,3) = 0;
% A1_sync_2(2,3) = 0;


%% Reduced form matrices
D_sync_1 = A0_sync_1\C_sync_1;
D_sync_2 = A0_sync_2\C_sync_2;

% B_sync_1 is A0_sync_1\A1_sync_1(1:nvars,1:nvars) for all lags - verified
B_sync_1 = A0_sync_1\A1_sync_1;
B_sync_2 = A0_sync_2\A1_sync_2;

O_sync_1 = A0_sync_1\SIG_sync_1;
O_sync_2 = A0_sync_2\SIG_sync_2;

%% Collect matrices

params.A0_sync_1 = A0_sync_1;
params.A0_sync_2 = A0_sync_2;
params.C_sync_1 = C_sync_1;
params.C_sync_2 = C_sync_2;
params.A1_sync_1 = A1_sync_1;
params.A1_sync_2 = A1_sync_2;
params.SIG_sync_1 = SIG_sync_1;
params.SIG_sync_2 = SIG_sync_2;
params.D_sync_1 = D_sync_1;
params.D_sync_2 = D_sync_2;
params.B_sync_1 = B_sync_1;
params.B_sync_2 = B_sync_2;
params.O_sync_1 = O_sync_1;
params.O_sync_2 = O_sync_2;
params.a12 = params_in.a12(dd);
params.a21 = params_in.a21(dd);
if isfield(params_in,'b12')==1
    params.b12 = params_in.b12(dd);
    params.b21 = params_in.b21(dd);
    params.c12 = params_in.c12(dd);
    params.c21 = params_in.c21(dd);
end

