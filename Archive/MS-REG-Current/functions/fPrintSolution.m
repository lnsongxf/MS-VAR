
function [] = fPrintSolution(outparams)
fprintf('\n *** A0 Matrix: State = 1 *** \n')
disp(outparams.A0_sync_1)    
fprintf('\n *** A0 Matrix: State = 2 *** \n')
disp(outparams.A0_sync_2)    

fprintf('\n *** A1 Matrix: State = 1 *** \n')
disp(outparams.A1_sync_1)    
fprintf('\n *** A1 Matrix: State = 2 *** \n')
disp(outparams.A1_sync_2)    


fprintf('\n *** C Matrix: State = 1 *** \n')
disp(outparams.C_sync_1)    
fprintf('\n *** C Matrix: State = 2 *** \n')
disp(outparams.C_sync_2)    

fprintf('\n *** SIG Matrix: State = 1 *** \n')
disp(outparams.SIG_sync_1)    
fprintf('\n *** SIG Matrix: State = 2 *** \n')
disp(outparams.SIG_sync_2) 


fprintf('\n *** TP Params: State = 1 *** \n')
disp([outparams.a12 outparams.b12 outparams.c12])    
fprintf('\n *** TP Params: State = 2 *** \n')
disp([outparams.a21 outparams.b21 outparams.c21])    
