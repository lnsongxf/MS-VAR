% This function implements the bootstrap proposed in Theorem 4
function result = CVfinalbootstrapInoue(el, bootMC, pit, rvec)
% INPUTS
% el - block length for the bootstrap
% bootMC - number of bootstrap replications
% pit - a vector of probability integral transforms, P by 1 vector, P being
% the effective sample size
% rvec - distritization used for the domain of the uniform distribritution
% (1 by x vector, x determined by the discretization used)
% OUTPUT
% results, a vector 3 by 2, the first column has the critical values of the
% Kolmogorov-Smirnov test, while the second one for the Cramer-von-Mises
% test. The rows are corresponding to 1%, 5% and 10% critical values,
% respectively.

KSv = zeros(bootMC,1); CVMv = zeros(bootMC,1);  

P = size(pit,1);
size_rvec = length(rvec);

for bootrep = 1:bootMC
    z = 1/sqrt(el)*randn(P-el+1,1);
    emp_cdf = (repmat(pit,1,size_rvec) <= repmat(rvec,P,1));
    mean_emp_cdf = mean(emp_cdf,1);
    
    K_star = zeros(1,size_rvec);
    for j = 1:P-el+1
         K_star = K_star + (P^(-1/2)*z(j,1)).*(sum(emp_cdf(j:j+el-1,:) - repmat(mean_emp_cdf,el,1),1));
    end
    
    KSv(bootrep,1) = max(abs(K_star'));
    CVMv(bootrep,1) = mean(K_star'.^2);
end

KSv = sort(KSv,'ascend');       cvKv = KSv(bootMC*[0.99 0.95 0.90]);
CVMv = sort(CVMv,'ascend');     cvMv = CVMv(bootMC*[0.99 0.95 0.90]); 

result = [cvKv cvMv];