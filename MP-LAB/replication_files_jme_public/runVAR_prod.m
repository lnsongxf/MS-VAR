%% --------------------- constant-parameter model --------------------- %%

%% Housekeeping
clear; clc; tic;
% close all;

% RISE Toolbox
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');

% Load RISE
rise_startup()

prod_test = 1;
plotfigures = 1;
savefigures = 0;
bands = [50 16 84]; % Coverage bands

%% Opts
opt.mp_restriction = 2;      % 1 = restrictions of as in Arias et al. (2019), 2 = coeffs fixed at their estimations, 3 = test
opt.paramsUse      = 'mode'; % mode = posterior mode, mcmc = parameter uncertainty
opt.force_estimate = 0;      % 1 = force to reestimate model

%% Create dataset

[db,varlist0]=createData;
varlist=fieldnames(varlist0);

%% Create the VAR

clc
% from JME code:
% nlag = 12; % number of lags
% nex  = 0;  % set equal to 1 if a constant is included; 0 otherwise
nlags=12;
exog={};
constant=false;
sv0=svar(varlist,exog,nlags,constant);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
if opt.mp_restriction ==1
    lin_restr={
        % Restriction 1
        % systematic component of MP
        'a0(6,TR)=0'
        'a0(6,NBR)=0'
        };
    nonlin_restr={
        % Restriction 2
        % systematic component of MP
        'a0(6,GDP)<0' % Negative, as it is on the LHS
        'a0(6,PGDP)<0' % Negative, as it is on the LHS
        };
elseif opt.mp_restriction ==2
    lin_restr={
        % Restriction 1
        % systematic component of MP
        'a0(6,TR)=0'
        'a0(6,NBR)=0'
        's(6,6)=1/1.06'
        'a0(6,GDP)=-0.85' % Negative, as it is on the RHS
        'a0(6,PGDP)=-2.33' % Negative, as it is on the RHS
%         'a0(6,GDP)=-0.89' % Negative, as it is on the RHS
%         'a0(6,PGDP)=-2.48' % Negative, as it is on the RHS
        };
    nonlin_restr={
%         'a0(6,GDP)>s(6,6)' % Negative, as it is on the RHS
%         'a0(6,PGDP)>s(6,6)' % Negative, as it is on the RHS
        };
elseif opt.mp_restriction ==3
    lin_restr={
        % Restriction 1
        % systematic component of MP
        'a0(6,TR)=0'
        'a0(6,NBR)=0'
        };
    nonlin_restr={
        % Restriction 2
        % systematic component of MP
        'a0(6,GDP)<0' % Negative, as it is on the LHS
        'a0(6,PGDP)<0' % Negative, as it is on the LHS

        % Arias etal. (2019) setup:
        % s(6,FFR)*FFR_t + a0(6,GDP)*GDP_t + a0(6,PGDP)*PGDP_t  = e_t
        % 1*FFR_t = - s(6,FFR)^(-1)*a0(6,GDP)*GDP_t - s(6,FFR)^(-1)*a0(6,PGDP)*PGDP_t  + s(6,FFR)^(-1)*e_t

        % RISE setup:
        % 1*FFR_t + a0(6,GDP)*GDP_t + a0(6,PGDP)*PGDP_t  = s(6,FFR)* e_t
        % 1*FFR_t = - a0(6,GDP)*GDP_t - a0(6,PGDP)*PGDP_t  + s(6,FFR)* e_t
        };
end

restrictions=[lin_restr;nonlin_restr];

%% set priors

var_prior=svar.prior_template();
var_prior.type='sz';
prior=struct('var',var_prior);

%% Find posterior mode
clc
sv=sv0;
startdb = '1965M1';
enddb   = '2007M6';
if strcmp(opt.paramsUse,'mode')
    if exist(['sPMode' num2str(opt.mp_restriction) '.mat'],'file')==0 || opt.force_estimate==1
        sv=estimate(sv,db,{startdb,enddb},prior,restrictions);
        % sv=estimate(sv,db,{startdb,enddb},prior,restrictions,'fmincon');

        % Structure to store posterior mode results
        sPMode.prior        = prior;
        sPMode.restrictions = restrictions;
        sPMode.startdb      = startdb;
        sPMode.enddb        = enddb;
        sPMode.sv           = sv;
        sPMode.nlags        = nlags;
        sPMode.exog         = exog;
        sPMode.constant     = constant;

        % Store solution structure
        save(['sPMode' num2str(opt.mp_restriction) '.mat'],'sPMode');
    else
        load(['sPMode' num2str(opt.mp_restriction) '.mat']);
        sv = sPMode.sv;
    end
else
    %%%%%%%%
    % MCMC %
    %%%%%%%%

    % Load or estimate mcmc
    if exist(['sPMCMC' num2str(opt.mp_restriction) '.mat'],'file')==0 || opt.force_estimate==1
        % Estimate posterior distribution
        run scriptEstimateMCMC_prod.m

        % Store solution structure
        save(['sPMCMC' num2str(opt.mp_restriction) '.mat'],'sPmcmc')

    else

        % Load stored results
        load(['sPMCMC' num2str(opt.mp_restriction) '.mat']);

        % Map objects:
        sv = sPmcmc.sv;
        pmcmc = sPmcmc.params;
        results = sPmcmc.results;
        [ff,lb,ub,x0,vcov,self]=pull_objective(sv);

        % Collect Posterior Estimates and Confidence Bands
        pmode=posterior_mode(sv);

        pnames=fieldnames(pmode);

        a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

        params=[results.pop.x];

        params_M = prctile(params,50,2);
        params_UB = prctile(params,95,2);
        params_LB = prctile(params,5,2);

        print_structural_form(sv)

        F = fieldnames(pmode);
        C = num2cell(struct2cell(pmode),3);
        C = cellfun(@(c)[c{:}],C,'uni',0);
        params_pmode = vertcat(C{:});

        params_M_Struct=a2tilde_to_a(params_M);
        params_UB_Struct=a2tilde_to_a(params_UB);
        params_LB_Struct=a2tilde_to_a(params_LB);

        paramsAll = a2tilde_to_a(params);
        for ii=1:numel(F)
            eval(['pmcmc.' F{ii} '= paramsAll(ii,:) ;'])
        end

        %% Collect Parameters
        paramsCI_mode = [params_pmode params_pmode-(params_M_Struct-params_LB_Struct)  params_pmode+(params_UB_Struct-params_M_Struct)];
        paramsCI_mcmc = [params_M_Struct params_LB_Struct  params_UB_Struct];

        %% Print Table of Estimates
        %         pmcmc_tex(pmode,paramsCI_mcmc,texfolder,modelname,'MCMC');
        %
        %         pmcmc_tex(pmode,paramsCI_mode,texfolder,[modelname '_mode'],'Mode');


    end

end


%% Printing estimates
clc

% Printing structural form
print_structural_form(sv)

% % Printing reduced form
% print_solution(sv)

%% Impulse responses

myirfs=irf(sv,[],61);

snames=fieldnames(myirfs);

ci=[16,50,68];

figure('name','Simple impulse response functions')
iter=0;
% focus on MP shock
shk=snames{end};
for iv=1:numel(varlist)
    v=varlist{iv};
    iter=iter+1;
    subplot(3,2,iter)
    ax = gca;
    plot(myirfs.(shk).(v),'linewidth',2)
    %     d=myirfs.(shk);
    %     out=fanchart(d,ci);
    %     plot_fanchart(out)
    set(ax,'XTick',0:12:60)
    set(ax,'XTickLabel',0:5)
    xlabel('Years')
    title(v)
end

%% IRFs
if prod_test ==1

    if strcmp(opt.paramsUse,'mode')
        pmode=posterior_mode(sv);
        params_temp = pmode;
        rps = size(pmode.a0_1_1,2);
    else
        params_temp = pmcmc;
        rps = size(pmcmc.a0_1_1,2);
    end
    nimp = 60;
    nvar = size(varlist,1);
    outs = zeros(1,rps,nimp,nvar);
    nshocks = 6;

    % Equation is:
    % (A0)*y_t = (Aplus)*y_{t-1} + SIG*e_t
    % IRF is:
    % y_t = inv(A0)*(Aplus)*y_{t-1} + inv(A0)*SIG*e_t

    parfor iii=1:rps
%         waitbar(iii/(rps),h,'Calculating IRFs...')
        params_dr = scriptParams_FULL_companion_fixed(sv,params_temp,iii);
        A0_dr = params_dr.A0_sync_1;
        Aplus_dr = params_dr.A1_sync_1;
        invA0 = (A0_dr)\eye(nvar);  % CORRECT
        IMP_dr= invA0*(params_dr.SIG_sync_1); % CORRECT
        nbeta_dr = (invA0)*(Aplus_dr); % CORRECT
        A_companion_dr=zeros(nvar*nlags,nvar*nlags);
        A_companion_dr(1:nvar,:)=nbeta_dr(:,1:nvar*nlags);
        A_companion_dr(nvar+1:nvar*nlags,1:nvar*nlags-nvar)=eye(nvar*nlags-nvar);
        invIMP_dr = (IMP_dr)\eye(size(IMP_dr,2)); % inv(IMP_dr)
        U1_dr=[IMP_dr; zeros(nvar*nlags-nvar,nvar)];
        nnn_dr = size(U1_dr,1);
        Eye_comp_dr = eye(nnn_dr)-A_companion_dr;
        Zk1_dr = zeros(nvar,nimp,nvar*nlags);
        impulse1_dr = zeros(nvar,nimp,nvar);
        for r=nshocks
            for k=1:nimp
                Zk1_dr(r,k,:)=(A_companion_dr^(k-1)*U1_dr(:,r))';
            end
            impulse1_dr(r,:,1:end)=Zk1_dr(r,:,1:nvar);
            if impulse1_dr(r,1,r)<0; impulse1_dr(r,:,:) = impulse1_dr(r,:,:)*(-1);end
            for ss=1:nvar
                outs(r,iii,:,ss)=impulse1_dr(r,:,ss);
            end
        end
    end
%     parfor iii=1:rps
%         params_dr = scriptParams_FULL_companion_fixed(sv,params_temp,iii);
%         IMP_dr=params_dr.O_sync_1;
%         A_companion_dr=zeros(nvar*nlags,nvar*nlags);
%         nbeta_dr = params_dr.B_sync_1;
%         A_companion_dr(1:nvar,:)=nbeta_dr(:,1:nvar*nlags);
%         A_companion_dr(nvar+1:nvar*nlags,1:nvar*nlags-nvar)=eye(nvar*nlags-nvar);
%         invIMP_dr = (IMP_dr)\eye(size(IMP_dr,2)); % inv(IMP_dr)
%         U1_dr=[IMP_dr; zeros(nvar*nlags-nvar,nvar)];
%         nnn_dr = size(U1_dr,1);
%         Eye_comp_dr = eye(nnn_dr)-A_companion_dr;
%         invEye_comp_dr = (Eye_comp_dr)\eye(size(Eye_comp_dr,2)); % inv(Eye_comp)
%         Zk1_dr = zeros(nvar,nimp,nvar*nlags);
%         impulse1_dr = zeros(nvar,nimp,nvar);
%         for r=nshocks
%             for k=1:nimp
%                 Zk1_dr(r,k,:)=(A_companion_dr^(k-1)*U1_dr(:,r))';
%             end
%             impulse1_dr(r,:,1:end)=Zk1_dr(r,:,1:nvar);
%             if impulse1_dr(r,1,r)<0; impulse1_dr(r,:,:) = impulse1_dr(r,:,:)*(-1);end
%             for ss=1:nvar
%                 outs(r,iii,:,ss)=impulse1_dr(r,:,ss);
%             end
%         end
%     end

    % Figures
    if plotfigures ==1
        Yname = varlist;
        for ee=1:length(nshocks)
            pp = nshocks(ee);
            SIZE_VARS = size(Yname,1); p_lines = floor(SIZE_VARS/3); p_cols = ceil(SIZE_VARS/p_lines); selec = 1:nvar; xaxis = 1:nimp; xaxis = xaxis';
            h=figure('Units','normalized','Color',[0.9412 0.9412 0.9412],'outerposition',[0,0,0.8,0.8],'Name',eval(['''' cell2mat(Yname(selec(pp))) '-shock''']));
            figure(h);
            for uuu=1:nvar
                subplot(p_lines,p_cols,uuu)
                temp=squeeze(outs(pp,:,:,uuu));
                if strcmp(opt.paramsUse,'mcmc')
                    temp1=squeeze(prctile(temp,bands,1))';
                    grpyat = [(1:nimp)' temp1(1:nimp,2); (nimp:-1:1)' temp1(nimp:-1:1,3)];
                    patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]);  hold on;
                else
                    temp1=temp;
                end
                ss(1) = plot(xaxis,temp1(1:nimp,1),'k-','LineWidth',3); hold on;
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
    end

    %%% SAVE FIGURES
    if savefigures ==1
        h = get(0,'children'); h = sort(h);
        for i=1:length(h); saveas(h(i), get(h(i),'Name'), 'png'); end
    end
end