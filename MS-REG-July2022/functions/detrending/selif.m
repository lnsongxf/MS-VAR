function xs=selif(x,t)
% SELIF 
% SELIF(x,t) selects the elements of x for which t=1
% x: NxK matrix, t: Nx1 matrix of 0's and 1's
nt=repmat(t,1,rows(x'));
nx=x.*nt;
nx_keep = nx(nt==1);
xs=reshape(nx_keep,sum(t),rows(x'));
