% Export required pmode structure elements to csv and then run R script that creates tex file.
% Writen by: Dylan Munson and Pablo Cuba-Borda. 
% Washington DC, Feb-23-2021
% ==============================================================================================

function pmode_tex(pmode,targetDir,targetName,tableTitle)

names=fieldnames(pmode);
pmode_cell=cell(2,length(names));

for i=1:length(names)
    pmode_cell{1,i}=names{i};
    pmode_cell{2,i}=round(pmode.(names{i}),4);
end

filename=[targetDir targetName '.csv'];
writecell(pmode_cell,filename);

%Run R script that compiles Tex document
eval(['!Rscript ./textools/pmode_create_table.R ' targetDir ' ' targetName ' ' tableTitle]);

% Clean up
delete([targetDir '*.log']);
delete([targetDir '*.aux']);
delete([targetDir '*.gz']);
delete([targetDir '*.fls']);
delete([targetDir '*.fdb_latexmk']);
