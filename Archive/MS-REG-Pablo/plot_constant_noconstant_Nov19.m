clc;
%%
load data_foreign_no_constant.mat

outMar2020 = fGetDensityMS('2020-Mar',data_use_nc,a2tilde_to_a_nc,results_nc,pnames_nc,nDraws_nc,nParamDraws_nc);
outSep2020 = fGetDensityMS('2020-Sep',data_use_nc,a2tilde_to_a_nc,results_nc,pnames_nc,nDraws_nc,nParamDraws_nc);
outOct2020 = fGetDensityMS('2020-Oct',data_use_nc,a2tilde_to_a_nc,results_nc,pnames_nc,nDraws_nc,nParamDraws_nc);

