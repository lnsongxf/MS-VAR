        
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         % RISE Toolbox


load inputirf.mat

addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts


% [quantiles, pdensity] = fPredDensityIteratedIRFs(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype);

hh = opt.hh;
nlags = sv.nlags;
[~,~,~,f]=filter(sv);
p_reg2_updated  = f.updated_regime_probabilities.regime2.data;


%------------------------
% INITIAL CONDITIONS: Y0
%------------------------
tt = 100;

% Macro and financial factor in period t
f0 = FF(tt);
m0 = MF(tt);
gdp0 = GDP(tt);


%------------------------
% INITIAL STATE: s0
%------------------------

s0 = 1;

% udraw = rand(1); % draw a coin that determines the regime

% if udraw < (1-pbad(tt))
%     st_sim = 1;
% else
%     st_sim = 2; % bad regime
% end


%---------------------------
% DRAW SEQUENCE OF UNIFORMS
%---------------------------


upath = rand(1,hh); % draw a coin that determines the regime


%-----------------------
% SIMULATE BASELINE PATH 
%-----------------------
init.s0 = s0;
init.f0 = f0;
init.m0 = m0;
init.y0 = gdp0;

% Draws EPS

[y_out,s_out] = simulate_IteratedFull(init,shocks,param,opt);

