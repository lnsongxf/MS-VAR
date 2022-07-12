% Export required pmode structure elements to csv and then run R script that creates tex file.
% Writen by: Dylan Munson and Pablo Cuba-Borda. 
% Washington DC, Feb-23-2021
% ==============================================================================================

function pmcmc_tex(pmode,paramsCI,targetDir,targetName,tableTitle)

names=fieldnames(pmode);
pmode_cell=cell(4,length(names));

for ii=1:length(names)
    pmode_cell{1,ii}=names{ii};
    pmode_cell{2,ii}=round(paramsCI(ii,1),4);
    pmode_cell{3,ii}=round(paramsCI(ii,2),2);
    pmode_cell{4,ii}=round(paramsCI(ii,3),2);
end

filename=[targetDir targetName '.csv'];
writecell(pmode_cell,filename);

%Run R script that compiles Tex document
eval(['!Rscript ./textools/pmcmc_create_table2.R ' targetDir ' ' targetName ' ' tableTitle]);

% Clean up
delete([targetDir '*.log']);
delete([targetDir '*.aux']);
delete([targetDir '*.gz']);
delete([targetDir '*.fls']);

delete([targetDir '*.fdb_latexmk']);
