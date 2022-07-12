function OOS_slurm_RS(i)

    %% Important paths
    if ~isdeployed
        addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         % RISE Toolbox
        addpath(genpath('functions'));      % Main functions
        addpath(genpath('scripts'));        % Main scripts
        addpath(genpath('textools'));       % Tools to produce tables
        addpath(genpath('auxtools'));       % Tools for plotting and other
        addpath(genpath('cbrewer'));        % Color mixing tools
    end
    
    %% Load RISE
    rise_startup()
    
    %% Execute program
    msreg_2reg_endoprob_SVAR_FULL_OOS(i);
    
    %% Exit RISE
    rise_exit
end
