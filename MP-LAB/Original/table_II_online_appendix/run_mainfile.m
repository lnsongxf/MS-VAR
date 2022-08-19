clear variables;
clc
close all






%% Table II Appendix
disp('----------------------------------------------------------')
disp('**********************************************************')
disp('************************* TABLE II Appendix **************')
disp('**********************************************************')
disp('----------------------------------------------------------')



disp('********************** Restrictions 2 and 3 *********************')
cd ..
gcd = pwd;
eval(['load ',gcd,'/figure_XV_panel_c/results/results.mat']);
cd table_II_online_appendix

A0tilde_11 = squeeze(A0tilde(1,1,:));
A0tilde_21 = squeeze(A0tilde(2,1,:));
A0tilde_31 = squeeze(A0tilde(3,1,:));
A0tilde_41 = squeeze(A0tilde(4,1,:));
A0tilde_51 = squeeze(A0tilde(5,1,:));
A0tilde_61 = squeeze(A0tilde(6,1,:));

%response to output
ppsiy    = -A0tilde_11./A0tilde_61;
%response to prices
ppsip    = -A0tilde_21./A0tilde_61;
%response to commodity prices
ppsipc   = -A0tilde_31./A0tilde_61;
%response total reserves
ppsitr   = -A0tilde_41./A0tilde_61;
%response non-borrowed reserves
ppsinbr  = -A0tilde_51./A0tilde_61;


% probability that ppsiy<0
probnegppsiy = size(find(ppsiy<0),1)/size(ppsiy,1);
% probability that ppsippi<0
probnegppsip = size(find(ppsip<0),1)/size(ppsip,1);
% probability that ppsipc<0
probnegpc    = size(find(ppsipc<0),1)/size(ppsipc,1);


disp('Probability of Satisfying Restriction 1: \psi_{tr}~=0 or \psi_{nbr}~=0')
disp(1-(size(find([abs(ppsitr)<1e-12 & abs(ppsitr)<1e-12]),1)/size(ppsitr,1)));

disp('Probability that \psi_{y}<0')
disp(size(find(ppsiy<0),1)/size(ppsiy,1));

disp('Probability that \psi_{p}<0')
disp(size(find(ppsip<0),1)/size(ppsip,1));

disp('Probability that \psi_{y}<0 | \psi_{p}<0')
disp(size(find([ppsiy<0 | ppsip<0]),1)/size(ppsiy,1));

