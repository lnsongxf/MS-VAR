%==========================================================================
%% housekeeping
%==========================================================================
clear variables;
close all;
userpath('clear');
clc;

rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed

currdir=pwd;
cd ..
get_help_dir_currdir=pwd;
addpath([get_help_dir_currdir,'/helpfunctions']); % set path to helper functions
cd(currdir)

prod_test = 1;
plotfigures = 1;
savefigures = 0;
bands = [50 16 84]; % Coverage bands
opt.force_estimate = 0;      % 1 = force to reestimate model


%==========================================================================
%% load the data and priors
%==========================================================================
% variables are in the following order: the description corresponds to the updated dataset
% 1: Real GDP               (Billions of Chained 2009 Dollars Seasonally Adjusted Annual Rate)
% 2: GDP Deflator           (Index 2009=100)
% 3: PCOM                   (Dow Jones Spot Average (Symbol _DJSD) Commodity Price Index, commercially available from Global Financial Data Index 1982-1984=100)
% 4: Total Reserves         (Billions of Dollars)
% 5: Non-borrowed Reserves  (Billions of Dollars)
% 6: Federal Funds Rate     (Annual Rate in Percentages)
data  = readtable([get_help_dir_currdir,'/data/csvfiles/dataset.csv']);


% data   = data(find(strcmp(data.dates,'1965-01-01')==1):find(strcmp(data.dates,'2007-06-01')==1),:);
data   = data(find(data.dates=='1965-01-01'):find(data.dates=='2007-06-01'),:);

%% data
% all variables are in log times 100 except for the federal funds rate that enters the SVAR in annualized percentages
num  = [100*log(data.monthly_GDP) 100*log(data.monthly_GDPDEF) 100*log(data.CPRINDEX) 100*log(data.TRARR) 100*log(data.BOGNONBR) data.FEDFUNDS];


%==========================================================================
%% model setup
%==========================================================================
nlag      = 12;                  % number of lags
nvar      = 6;                   % number of endogenous variables
nex       = 0;                   % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex;     % number of exogenous variables
nd        = 1e5;                 % number of orthogonal-reduced-form (B,Sigma,Q) draws
iter_show = 1e4;                 % display iteration every iter_show draws
horizon   = 60;                  % maximum horizon for IRFs
index     = [4 8 16 24 32 40];                  % define  horizons for the FEVD
horizons  = 0:0;                 % horizons to restrict
NS        = 1 + numel(horizons); % number of objects in F(A_{0},A_{+}) to which we impose sign and zero restrictios: F(THETA)=[A_{0};L_{0},...,L_{horizons}]
e         = eye(nvar);           % create identity matrix
maxdraws  = 1e4;                 % max number of importance sampling draws
conjugate = 'structural';        % structural or irfs or empty

%==========================================================================
%% identification: declare Ss and Zs matrices
%==========================================================================
% restrictions on A0  and/or IRFs

% sign restrictions
S = cell(nvar,1);
for ii=1:nvar
    S{ii}=zeros(0,nvar*NS);
end
ns1  = 4;
S{1} = zeros(ns1,nvar*NS);
S{1}(1,1)   = -1;
S{1}(2,2)   = -1;
S{1}(3,6)   =  1;
S{1}(4,nvar+6)   =  1;

% zero restrictions
Z=cell(nvar,1);
for i=1:nvar
    Z{i}=zeros(0,nvar*NS);
end

nz1  = 2;
Z{1} = zeros(nz1,nvar*NS);
Z{1}(1,4)=1;
Z{1}(2,5)=1;


%==========================================================================
%% Setup info
%==========================================================================
info=SetupInfo(nvar,m,Z,@(x)chol(x));

% ZF(A_{0},A_{+})
info.nlag     = nlag;
info.horizons = horizons;
info.ZF       = @(x,y)ZF(x,y);

% functions useful to compute the importance sampler weights
iw_info = info;
fs      = @(x)ff_h(x,iw_info);
r       = @(x)ZeroRestrictions(x,iw_info);

if strcmp(conjugate,'irfs')==1
    fo              = @(x)f_h(x,iw_info);
    fo_str2irfs     = @(x)StructuralToIRF(x,iw_info);
    fo_str2irfs_inv = @(x)IRFToStructural(x,iw_info);
    r_irfs          = @(x)IRFRestrictions_more_general(x,iw_info); 
end


% function useful to check the sign restrictions
fh_S_restrictions  = @(x)SF(x,iw_info,S);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [yt(t-1), ... , yt(t-nlag), constant];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
yt = num(nlag+1:end,:);
T  = size(yt,1);
xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nvar*(i-1)+1:nvar*i) = num((nlag-(i-1)):end-i,:) ;
end
if nex>=1
    xt(:,nvar*nlag+nex)=ones(T,1);
end
% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors


%% prior for reduced-form parameters
nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,nvar);
PphiBar             = zeros(nvar);



%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;




%% useful definitions
% definitions used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([nd,1]); % reduced-form lag parameters
Sigmadraws     = cell([nd,1]); % reduced-form covariance matrices
Qdraws         = cell([nd,1]); % orthogonal matrices
storevefh      = zeros(nd,1);  % volume element f_{h}
storevegfhZ    = zeros(nd,1);  % volume element g o f_{h}|Z
uw             = zeros(nd,1);  % unnormalized importance sampler weights

if strcmp(conjugate,'irfs')==1
    storevephi      = zeros(nd,1);  % volume element f_{h}
    storevegphiZ    = zeros(nd,1);  % volume element g o f_{h}|Z
end

% definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
F      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    F((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end

% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

if exist('results/results.mat','file')==0 || opt.force_estimate==1
%% initialize counters to track the state of the computations
counter = 1;
record  = 1;
count   = 0;
tStart = tic;
while record<=nd
    
    
    %% step 1 in Algorithm 2
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(PpsiTilde,nvar*m,1);
    Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;
    
   
    %% steps 2:4 of Algorithm 2
    w           = DrawW(iw_info);   
    x           = [vec(Bdraw); vec(Sigmadraw); w];
    structpara  = ff_h_inv(x,iw_info);
    
    % store the matrix Q associated with step 3
    Qdraw            = SpheresToQ(w,iw_info,Bdraw,Sigmadraw);
    Qdraws{record,1} = reshape(Qdraw,nvar,nvar);
    
   
    %% check if sign restrictions hold
    signs      = fh_S_restrictions(structpara);
    
    
    if (sum(signs{1}*e(:,1)>0))==size(signs{1}*e(:,1),1)
        
        count=count+1;
        
        %% compute importance sampling weights
       
        switch conjugate
            
            case 'structural'
                
                
                storevefh(record,1)   = (nvar*(nvar+1)/2)*log(2)-(2*nvar+m+1)*LogAbsDet(reshape(structpara(1:nvar*nvar),nvar,nvar));
                storevegfhZ(record,1) = LogVolumeElement(fs,structpara,r); 
                uw(record,1)          = exp(storevefh(record,1) - storevegfhZ(record,1));
                
            case 'irfs'
                
                irfpara                = fo_str2irfs(structpara);
                storevephi(record,1)   = LogVolumeElement(fo,structpara)   + LogVolumeElement(fo_str2irfs_inv,irfpara);
                storevegphiZ(record,1) = LogVolumeElement(fs,structpara,r) + LogVolumeElement(fo_str2irfs_inv,irfpara,r_irfs); 
                uw(record,1)           = exp(storevephi(record,1) - storevegphiZ(record,1));
                
            otherwise
                
                uw(record,1) = 1;
                
        end
        
    else
        
        uw(record,1) = 0;
        
    end
    
    if counter==iter_show
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(nd-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    record=record+1;
    
end


tElapsed = toc(tStart);
imp_w  = uw/sum(uw);
ne = floor(1/sum(imp_w.^2));


%% store draws
Ltilde        = zeros(horizon+1,nvar,nvar,ne); % define array to store IRF
A0tilde       = zeros(nvar,nvar,ne); % define array to store A0
Aplustilde    = zeros(m,nvar,ne); % define array to store Aplus
FEVD          = zeros(nvar,size(index,2),ne);                 % define array to store FEVD
hist_is_draws = zeros(ne,1);   % define array to store draws from importance sampler


for s=1:min(ne,maxdraws)
    
    %% draw: B,Sigma,Q
    is_draw = randsample(1:size(imp_w,1),1,true,imp_w);
    hist_is_draws(s,1)=is_draw;

    Bdraw       = Bdraws{is_draw,1};
    Sigmadraw   = Sigmadraws{is_draw,1};
    Qdraw       = Qdraws{is_draw,1};
    
    
    % compute matrix F: useful for FEVD
    hSigmadraw = hh(Sigmadraw);
    A0         = hSigmadraw\e;
    Aplus      = Bdraw*A0;
    % Page 8 ARRW
    for l=1:nlag-1
        A{l} = Aplus((l-1)*nvar+1:l*nvar,1:end);
        F((l-1)*nvar+1:l*nvar,1:nvar)=A{l}/A0;
    end
    A{nlag} = Aplus((nlag-1)*nvar+1:nlag*nvar,1:end);
    F((nlag-1)*nvar+1:nlag*nvar,:)=[A{nlag}/A0 extraF];
    
    x=[reshape(Bdraw,m*nvar,1); reshape(Sigmadraw,nvar*nvar,1); Qdraw(:)];
    structpara = f_h_inv(x,info);
    
    
    LIRF =IRF_horizons(structpara, nvar, nlag, m, 0:horizon);
    
    
    for h=0:horizon
        Ltilde(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
    end

    %     for j=1:size(index,2)
    %     FEVD(:,j,s)         = variancedecomposition(F',J,Sigmadraw, hh(Sigmadraw)'*Qdraw(:,1),nvar,index(1,j));
    %     end

    A0tilde(:,:,s)    = reshape(structpara(1:nvar*nvar),nvar,nvar);
    Aplustilde(:,:,s) = reshape(structpara(nvar*nvar+1:end),m,nvar);



end

A0tilde    = A0tilde(:,:,1:s);
Aplustilde = Aplustilde(:,:,1:s);
Ltilde     = Ltilde(:,:,:,1:s);
% FEVD       = FEVD(:,:,1:s);



cd results
savefile='results.mat';
save(savefile,'Ltilde','A0tilde','Aplustilde','imp_w','ne','FEVD');
cd ..
else
    load('results/results.mat')
end

%% IRFs
if prod_test ==1


    rps = size(A0tilde,3);
    varlist = [{'GDP' };    {'PGDP'};    {'PCOM'};    {'TR'  };    {'NBR' };    {'FFR' }];
    shocklist = {'FFR' };
    nimp = 60;
    nvar = size(varlist,1);
    nlags = nlag;
    outs = zeros(1,rps,nimp,nvar);
    nshocks = 1;

    % Equation is:
    % (A0')*y_t = (Aplus')*y_{t-1} + e_t
    % IRF is:
    % y_t = inv(A0')*(Aplus')*y_{t-1} + inv(A0')*e_t

    nbetas = NaN(nvar,nvar*nlags,rps);
    imps =  NaN(nvar,nvar,size(A0tilde,3));
    parfor iii=1:rps
        A0_dr = squeeze(A0tilde(:,:,iii));
        Aplus_dr = squeeze(Aplustilde(:,:,iii));
        invA0 = (A0_dr')\eye(nvar);  % CORRECT
        IMP_dr= (A0_dr')\eye(nvar); % CORRECT
        nbeta_dr = (invA0)*(Aplus_dr'); % CORRECT
        nbetas(:,:,iii) = nbeta_dr;
        imps(:,:,iii) = IMP_dr;
        A_companion_dr=zeros(nvar*nlags,nvar*nlags);
        A_companion_dr(1:nvar,:)=nbeta_dr(:,1:nvar*nlags);
        A_companion_dr(nvar+1:nvar*nlags,1:nvar*nlags-nvar)=eye(nvar*nlags-nvar);
        U1_dr=[IMP_dr; zeros(nvar*nlags-nvar,nvar)];
        nnn_dr = size(U1_dr,1);
        Eye_comp_dr = eye(nnn_dr)-A_companion_dr;
        invEye_comp_dr = (Eye_comp_dr)\eye(size(Eye_comp_dr,2)); % inv(Eye_comp)
        Zk1_dr = zeros(nvar,nimp,nvar*nlags);
        impulse1_dr = zeros(nvar,nimp,nvar);
        for r=nshocks
            for k=1:nimp
                Zk1_dr(r,k,:)=(A_companion_dr^(k-1)*U1_dr(:,r))';
            end
            impulse1_dr(r,:,1:end)=Zk1_dr(r,:,1:nvar);
            if impulse1_dr(6,1,6)<0; impulse1_dr(r,:,:) = impulse1_dr(r,:,:)*(-1);end
            for ss=1:nvar
                outs(r,iii,:,ss)=impulse1_dr(r,:,ss);
            end
        end
    end
    nbeta_mode = squeeze(quantile(nbetas,0.5,3));
    imp_mode =  squeeze(quantile(imps,0.5,3));
    % Figures
    if plotfigures ==1

        %%***************************************************************
        %                         HOUSEKEEPING                          *
        %%***************************************************************
        Ltildeq50=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 50th
        Ltildeq16=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 16th
        Ltildeq84=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 84th

        Ltildeq025=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 025th
        Ltildeq975=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 975th


        for ii=1:size(Ltilde,1)
            for jj=1:size(Ltilde,2)
                for kk=1:size(Ltilde,3)
                    Ltildeq50(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.5);
                    Ltildeq16(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.16);
                    Ltildeq84(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.84);
                    Ltildeq025(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.025);
                    Ltildeq975(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.975);
                end
            end
        end

        save("Ltildes.mat","Ltildeq50","Ltildeq025","Ltildeq975","Ltildeq84","Ltildeq16");


        Yname = varlist;
        for ee=1:length(nshocks)
            pp = nshocks(ee);
            SIZE_VARS = size(Yname,1); p_lines = floor(SIZE_VARS/3); p_cols = ceil(SIZE_VARS/p_lines); selec = 1:nvar; xaxis = 1:nimp; xaxis = xaxis';
            h=figure('Units','normalized','Color',[0.9412 0.9412 0.9412],'outerposition',[0,0,0.8,0.8],'Name',eval(['''' cell2mat(shocklist(selec(pp))) '-shock''']));
            figure(h);
            for uuu=1:nvar
                subplot(p_lines,p_cols,uuu)
                temp=squeeze(outs(pp,:,:,uuu));
                temp1=squeeze(prctile(temp,bands,1))';
                grpyat = [(1:nimp)' temp1(1:nimp,2); (nimp:-1:1)' temp1(nimp:-1:1,3)];
                patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]);  hold on;
                ss(1) = plot(xaxis,temp1(1:nimp,1),'k-','LineWidth',3); hold on;
                ss(2) = plot(xaxis,squeeze(Ltildeq50(1:nimp,uuu,1)),'-b','LineWidth',1);
                plot(xaxis,zeros(nimp),':k');
                set(gca,'FontSize',16)
                title(Yname(uuu),'FontSize',16) %'FontWeight','bold',
                ylabel('percent','FontSize',16)
                xlabel('months','FontSize',16)
                set(gca,'XTick',0:4:nimp)
                axis tight
                hold off;
            end
        end
        %         run results/run_plot_prod.m
    end

    %%% SAVE FIGURES
    if savefigures ==1
        h = get(0,'children'); h = sort(h);
        for i=1:length(h); saveas(h(i), get(h(i),'Name'), 'png'); end
    end
end


