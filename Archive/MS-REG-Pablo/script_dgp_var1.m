clear; clc
% Structural form
% A(st)*X(t) = C(st) + B(st)*X(t-1) + Z(st)*eps(t)

% Reduced form
% X(t) = D(st) + E(st)*X(t-1) + O(st)*eta(t)


A_0 = [ 1 a_0_12 a_0_13; 0 1 0; 0 0 1];
A_1 = [ 1 a_1_12 a_1_13; 0 1 0; 0 0 1];

C_0 = [mu_0; 0 ; 0];
C_1 = [mu_1; 0 ; 0];

B_0 = [ 0 b_0_12 b_0_13; 0 b22 0; 0 0 b33];
B_1 = [ 0 b_1_12 b_1_13; 0 b22 0; 0 0 b33];
Z = [sig1 0 0 ; 0 sig2 0; 0 0 sig3];

% Reduced form 
D = inv(A)*C;
E = inv(A)*B;
O = inv(A)*Z;

