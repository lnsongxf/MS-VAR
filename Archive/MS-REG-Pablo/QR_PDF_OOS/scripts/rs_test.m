function rs_test(z, rvec)
% INPUT
% z - T by 1 vector of the probability integral transforms
% rvec - 1 by x vector with values between 0 and 1 guiding the
% discretization of the uniform
% *************************************************************************
% Implementing the test
% *************************************************************************
% Creating the test statistics
% cumcumz corresponds to the xi and phi on page 8, respectively
T = size(z,1);
cumcumz = (repmat(z,1,size(rvec,2)) < repmat(rvec,size(z,1),1)) - repmat(rvec,size(z,1),1);
v = sum(cumcumz,1)/sqrt(T);

KS = max(abs(v)); % Kolmogorov-Smirnov test statistics
CVM = mean(v.^2); % Cramer-von-Mises test statistics

% If implementing the test as in Theorem 2, which is applicable for
% one-step-ahead density forecasts, the following two lines of code apply.
% The critical values are as in Table 1 of the paper and correspond to 1% 5% and 10%,
% respectively
QVrej_one_step_ahead = (KS > [1.61 1.34 1.21]);   % Kolmogorov-Smirnov test 
CVMrej_one_step_ahead = (CVM > [0.74 0.46 0.35]);% Cramer-von-Mises test 

% Simulate the critical values applicable for multi-step-ahead density
% forecast calibration
% Implements the test based on the Theorem 3
% setup some preliminaries for the bootstrap
el = floor(T.^(1/3)); % block length
bootMC = 200;

tableboot = CVfinalbootstrapInoue(el, bootMC, z, rvec);
QVrej_multi_step_ahead = (KS > tableboot(:,1))'; 
CVMrej_multi_step_ahead = (CVM > tableboot(:,2))';

% *************************************************************************
% Report the results
% *************************************************************************
disp('Test statistics [KS CVM]');
disp([KS CVM]);
disp('Critical Values, 1% 5% 10%')
disp('One-step-ahead density calibration test')
disp('Kolmogorov-Smirnov Test');
disp([1.61 1.34 1.21]);
disp('Cramer-von-Mises Test');
disp([0.74 0.46 0.35]);
disp('Multi-step-ahead density calibration test');
disp('Kolmogorov-Smirnov Test');
disp(tableboot(:,1)');
disp('Cramer-von-Mises Test');
disp(tableboot(:,2)');
