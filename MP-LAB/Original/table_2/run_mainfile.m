clear variables;
clc
close all





%% Table 2
disp('----------------------------------------------------------')
disp('**********************************************************')
disp('************************* TABLE 2 ************************')
disp('**********************************************************')
disp('----------------------------------------------------------')
cd ..
gcd = pwd;
eval(['load ',gcd,'/figure_3/results/results.mat']);
cd table_2

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

disp('********************** Restriction 3 *********************')
disp('Fed Funds Rate Contemporaneous Response to Output: 50-th Quantile')
disp(quantile(ppsiy,0.5))
disp('Fed Funds Rate Contemporaneous Response to Output: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsiy,0.16),quantile(ppsiy,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Output: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsiy,0.025),quantile(ppsiy,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Prices: 50-th Quantile')
disp(quantile(ppsip,0.5))
disp('Fed Funds Rate Contemporaneous Response to Prices: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsip,0.16),quantile(ppsip,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Prices: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsip,0.025),quantile(ppsip,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 50-th Quantile')
disp(quantile(ppsipc,0.5))
disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsipc,0.16),quantile(ppsipc,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsipc,0.025),quantile(ppsipc,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 50-th Quantile')
disp(quantile(ppsitr,0.5))
disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsitr,0.16),quantile(ppsitr,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsitr,0.025),quantile(ppsitr,0.975)])


disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 50-th Quantile')
disp(quantile(ppsinbr,0.5))
disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsinbr,0.16),quantile(ppsinbr,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsinbr,0.025),quantile(ppsinbr,0.975)])






cd ..
gcd = pwd;
eval(['load ',gcd,'/figure_4_panel_a/results/results.mat']);
cd table_2

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


disp('********************** Restrictions 1 and 3 *********************')
disp('Fed Funds Rate Contemporaneous Response to Output: 50-th Quantile')
disp(quantile(ppsiy,0.5))
disp('Fed Funds Rate Contemporaneous Response to Output: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsiy,0.16),quantile(ppsiy,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Output: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsiy,0.025),quantile(ppsiy,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Prices: 50-th Quantile')
disp(quantile(ppsip,0.5))
disp('Fed Funds Rate Contemporaneous Response to Prices: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsip,0.16),quantile(ppsip,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Prices: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsip,0.025),quantile(ppsip,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 50-th Quantile')
disp(quantile(ppsipc,0.5))
disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsipc,0.16),quantile(ppsipc,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Commodity Prices: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsipc,0.025),quantile(ppsipc,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 50-th Quantile')
disp(quantile(ppsitr,0.5))
disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsitr,0.16),quantile(ppsitr,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Total Reserves: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsitr,0.025),quantile(ppsitr,0.975)])

disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 50-th Quantile')
disp(quantile(ppsinbr,0.5))
disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 16-th Quantile and 84-th Quantile')
disp([quantile(ppsinbr,0.16),quantile(ppsinbr,0.84)])
disp('Fed Funds Rate Contemporaneous Response to Non-Borrowed Reserves: 2.5-th Quantile and 97.5-th Quantile')
disp([quantile(ppsinbr,0.025),quantile(ppsinbr,0.975)])



