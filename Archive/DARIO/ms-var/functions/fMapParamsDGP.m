function [param] = fMapParamsDGP(param,dd)

% Switching matrices
A0_sync_1 = [1 param.a0_1_2_sync_1(dd) param.a0_1_3_sync_1(dd);
    param.a0_2_1_sync_1(dd) 1 param.a0_2_3_sync_1(dd);
    param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];

A0_sync_2 = [1 param.a0_1_2_sync_2(dd) param.a0_1_3_sync_2(dd);
    param.a0_2_1_sync_2(dd) 1 param.a0_2_3_sync_2(dd);
    param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];

C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];

C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];

A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd) param.a1_1_3_sync_1(dd);
    param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)  param.a1_2_3_sync_1(dd);
    param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];

A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd) param.a1_1_3_sync_2(dd);
    param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)  param.a1_2_3_sync_2(dd);
    param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];

SIG_sync_1 = [param.s_1_1_sync_1(dd) 0 0;
    0 param.s_2_2_sync_1(dd) 0;
    0 0 param.s_3_3_sync_1(dd)];

SIG_sync_2 = [param.s_1_1_sync_2(dd) 0 0;
    0 param.s_2_2_sync_2(dd) 0;
    0 0 param.s_3_3_sync_2(dd)];

%% Reduced form matrices
D_sync_1 = A0_sync_1\C_sync_1;
D_sync_2 = A0_sync_2\C_sync_2;

B_sync_1 = A0_sync_1\A1_sync_1;
B_sync_2 = A0_sync_2\A1_sync_2;

O_sync_1 = A0_sync_1\SIG_sync_1;
O_sync_2 = A0_sync_2\SIG_sync_2;

%% Collect matrices

param.A0_sync_1 = A0_sync_1;
param.A0_sync_2 = A0_sync_2;
param.C_sync_1 = C_sync_1;
param.C_sync_2 = C_sync_2;
param.A1_sync_1 = A1_sync_1;
param.A1_sync_2 = A1_sync_2;
param.SIG_sync_1 = SIG_sync_1;
param.SIG_sync_2 = SIG_sync_2;
param.D_sync_1 = D_sync_1;
param.D_sync_2 = D_sync_2;
param.B_sync_1 = B_sync_1;
param.B_sync_2 = B_sync_2;
param.O_sync_1 = O_sync_1;
param.O_sync_2 = O_sync_2;

