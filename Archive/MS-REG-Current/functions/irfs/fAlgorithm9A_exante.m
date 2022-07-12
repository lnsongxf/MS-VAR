function [irf, pctl] = fAlgorithm9A_exante(init_in, shocks_in, param_in, opt_in, ps_bad_in,sd_base,sd_shock,varargin)


% Number of iteratios for regime shocks
niterU = opt_in.nDrawsU;

% Number of iterations for structural shocks
niterEPS = opt_in.nDrawsEPS;

% Map regime shocks
if isempty(varargin)
    u_shocks = rand(niterEPS,opt_in.hh);
else
    u_shocks = varargin{:};
end

%-----------------------------------
% Set shock matrices for baseline
%-----------------------------------
shocks.baseline.baseline = shocks_in;

% Baseline shocks for f-irf
shocks.f.baseline = shocks.baseline.baseline;

% Baseline shocks for m-irf
shocks.m.baseline = shocks.baseline.baseline;

% Baseline shocks for y-irf
shocks.y.baseline = shocks.baseline.baseline;


%-----------------------------------
% Set shocked paths of innovations
%-----------------------------------

shocks.f.shock = shocks.f.baseline;
shocks.f.shock(1,1,:) = shocks.f.shock(1,1,:) + sd_shock.f;
if opt_in.scenario
    shocks.f.shock(2,1,:) = shocks.f.shock(2,1,:) + sd_shock.m;
    shocks.f.shock(3,1,:) = shocks.f.shock(3,1,:) + sd_shock.y;
end

shocks.m.shock = shocks.m.baseline;
shocks.m.shock(2,1,:) = shocks.m.shock(2,1,:) + sd_shock.m;
if opt_in.scenario
    shocks.m.shock(1,1,:) = shocks.m.shock(1,1,:) + sd_shock.f;
    shocks.m.shock(3,1,:) = shocks.m.shock(3,1,:) + sd_shock.y;
end

shocks.y.shock = shocks.y.baseline;
shocks.y.shock(3,1,:) = shocks.y.shock(3,1,:) + sd_shock.y;
if opt_in.scenario
    shocks.y.shock(1,1,:) = shocks.y.shock(1,1,:) + sd_shock.f;
    shocks.y.shock(2,1,:) = shocks.y.shock(2,1,:) + sd_shock.m;
end


%Create matrix of s paths corresponding to the probability of the bad state
%that is passed to the function (ps_bad).  The probability in each column
%should correspond to the probability of being in the bad state that
%period.  Do this by looping over periods (the number of columns).



% Pre-allocate the memory
total_l = niterEPS*niterU;
path_f_baseline.f = NaN(total_l,opt_in.hh);
path_f_baseline.m = NaN(total_l,opt_in.hh);
path_f_baseline.y = NaN(total_l,opt_in.hh);
path_f_baseline.p12 = NaN(total_l,opt_in.hh);
path_f_baseline.p21 = NaN(total_l,opt_in.hh);
path_m_baseline.f = NaN(total_l,opt_in.hh);
path_m_baseline.m = NaN(total_l,opt_in.hh);
path_m_baseline.y = NaN(total_l,opt_in.hh);
path_m_baseline.p12 = NaN(total_l,opt_in.hh);
path_m_baseline.p21 = NaN(total_l,opt_in.hh);
path_y_baseline.f = NaN(total_l,opt_in.hh);
path_y_baseline.m = NaN(total_l,opt_in.hh);
path_y_baseline.y = NaN(total_l,opt_in.hh);
path_y_baseline.p12 = NaN(total_l,opt_in.hh);
path_y_baseline.p21 = NaN(total_l,opt_in.hh);
path_f_shock.f = NaN(total_l,opt_in.hh);
path_f_shock.m = NaN(total_l,opt_in.hh);
path_f_shock.y = NaN(total_l,opt_in.hh);
path_f_shock.p12 = NaN(total_l,opt_in.hh);
path_f_shock.p21 = NaN(total_l,opt_in.hh);
path_m_shock.f = NaN(total_l,opt_in.hh);
path_m_shock.m = NaN(total_l,opt_in.hh);
path_m_shock.y = NaN(total_l,opt_in.hh);
path_m_shock.p12 = NaN(total_l,opt_in.hh);
path_m_shock.p21 = NaN(total_l,opt_in.hh);
path_y_shock.f = NaN(total_l,opt_in.hh);
path_y_shock.m = NaN(total_l,opt_in.hh);
path_y_shock.y = NaN(total_l,opt_in.hh);
path_y_shock.p12 = NaN(total_l,opt_in.hh);
path_y_shock.p21 = NaN(total_l,opt_in.hh);
path_s_baseline.f = NaN(total_l,opt_in.hh);
path_s_baseline.m = NaN(total_l,opt_in.hh);
path_s_baseline.y = NaN(total_l,opt_in.hh);
path_s_shock.f = NaN(total_l,opt_in.hh);
path_s_shock.m = NaN(total_l,opt_in.hh);
path_s_shock.y = NaN(total_l,opt_in.hh);

%Loop over iterations and pass to fSimulateIterated each time.  Use one row
%of s_mat each time
ccount = 1;
count = 1;
for j = 1:niterEPS
    for i = 1:niterU
        
        % Pass regime shocks
        u_iter_shocks(1,1) = 0*u_shocks(j,1);
        u_iter_shocks(1,2:opt_in.hh) = 0*u_shocks(j,2:opt_in.hh);
        
        %Get "baseline" paths using initial shocks
        [Y_out_f_baseshock,  s_out_f_baseshock,p12_out_f_baseshock,p21_out_f_baseshock]   = fSimulateIterated(init_in,shocks.f.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_m_baseshock,  s_out_m_baseshock,p12_out_m_baseshock,p21_out_m_baseshock]   = fSimulateIterated(init_in,shocks.m.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_y_baseshock,  s_out_y_baseshock,p12_out_y_baseshock,p21_out_y_baseshock]   = fSimulateIterated(init_in,shocks.y.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
        
        %Get shocked path
        u_iter_shocks(1,1) = 0*u_shocks(j,1);
        u_iter_shocks(1,2:opt_in.hh) = u_shocks(j,2:opt_in.hh);

        % For switch in period 2
%          u_iter_shocks(1,2:opt_in.hh) = 1; %u_shocks(j,2:opt_in.hh);
%          u_iter_shocks(1,3:opt_in.hh) = 0; %u_shocks(j,2:opt_in.hh);
        
        [Y_out_f_shock,  s_out_f_shock,p12_out_f_shock,p21_out_f_shock]   = fSimulateIterated(init_in,shocks.f.shock(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_m_shock,  s_out_m_shock,p12_out_m_shock,p21_out_m_shock]   = fSimulateIterated(init_in,shocks.m.shock(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_y_shock,  s_out_y_shock,p12_out_y_shock,p21_out_y_shock]   = fSimulateIterated(init_in,shocks.y.shock(:,:,j),param_in,opt_in,u_iter_shocks);
        
        %Collect iterations into matrices and structures
        path_f_baseline.f(count,:) = Y_out_f_baseshock(1,:);
        path_f_baseline.m(count,:) = Y_out_f_baseshock(2,:);
        path_f_baseline.y(count,:) = Y_out_f_baseshock(3,:);
        path_f_baseline.p12(count,:) = p12_out_f_baseshock;
        path_f_baseline.p21(count,:) = p21_out_f_baseshock;
        
        path_m_baseline.f(count,:) = Y_out_m_baseshock(1,:);
        path_m_baseline.m(count,:) = Y_out_m_baseshock(2,:);
        path_m_baseline.y(count,:) = Y_out_m_baseshock(3,:);
        path_m_baseline.p12(count,:) = p12_out_m_baseshock;
        path_m_baseline.p21(count,:) = p21_out_m_baseshock;
        
        
        path_y_baseline.f(count,:) = Y_out_y_baseshock(1,:);
        path_y_baseline.m(count,:) = Y_out_y_baseshock(2,:);
        path_y_baseline.y(count,:) = Y_out_y_baseshock(3,:);
        path_y_baseline.p12(count,:) = p12_out_y_baseshock;
        path_y_baseline.p21(count,:) = p21_out_y_baseshock;
        
        path_f_shock.f(count,:) = Y_out_f_shock(1,:);
        path_f_shock.m(count,:) = Y_out_f_shock(2,:);
        path_f_shock.y(count,:) = Y_out_f_shock(3,:);
        path_f_shock.p12(count,:) = p12_out_f_shock;
        path_f_shock.p21(count,:) = p21_out_f_shock;
        
% %         if (path_f_shock.y(count,3) - path_f_baseline.y(count,3))>0 && ccount <10
% %         if ccount ==9
%             
%             citer = abs(rand(1))*[1 0.9 0.8];
%             figure(100);
%             subplot(2,2,1); hold on;
%             plot(path_f_shock.f(count,:) - path_f_baseline.f(count,:),'color',citer);
%             
%             subplot(2,2,2); hold on;
%             plot(path_f_shock.m(count,:) - path_f_baseline.m(count,:),'color',citer);
%             
%             subplot(2,2,3); hold on;
%             plot(path_f_shock.y(count,:) - path_f_baseline.y(count,:),'color',citer);
%             
%             subplot(2,2,4); hold on;
%             plot(s_out_f_shock,'color',citer);
%             legend(num2str(count)); legend boxoff
%             ccount = ccount+1;
% 
% %         end
%         
        path_m_shock.f(count,:) = Y_out_m_shock(1,:);
        path_m_shock.m(count,:) = Y_out_m_shock(2,:);
        path_m_shock.y(count,:) = Y_out_m_shock(3,:);
        path_m_shock.p12(count,:) = p12_out_m_shock;
        path_m_shock.p21(count,:) = p21_out_m_shock;
        
        path_y_shock.f(count,:) = Y_out_y_shock(1,:);
        path_y_shock.m(count,:) = Y_out_y_shock(2,:);
        path_y_shock.y(count,:) = Y_out_y_shock(3,:);
        path_y_shock.p12(count,:) = p12_out_y_shock;
        path_y_shock.p21(count,:) = p21_out_y_shock;
        
        % Collect path of s
        path_s_baseline.f(count,:) = s_out_f_baseshock';
        path_s_baseline.m(count,:) = s_out_m_baseshock';
        path_s_baseline.y(count,:) = s_out_y_baseshock';
        path_s_shock.f(count,:) = s_out_f_shock';
        path_s_shock.m(count,:) = s_out_m_shock';
        path_s_shock.y(count,:) = s_out_y_shock';
        
        count = count + 1;
    end
end

%Calculate IRFs as shocked - baseline
fields = fieldnames(path_f_baseline);
for i=1:length(fields)
    irf.f_shock.(fields{i}) = path_f_shock.(fields{i}) - path_f_baseline.(fields{i});
    irf.m_shock.(fields{i}) = path_m_shock.(fields{i}) - path_m_baseline.(fields{i});
    irf.y_shock.(fields{i}) = path_y_shock.(fields{i}) - path_y_baseline.(fields{i});
    
    % Compute percentiles
    pctl.f_shock.(fields{i}) = prctile(irf.f_shock.(fields{i}),[1 5 10 25 50 75 90 95 99]);
    pctl.m_shock.(fields{i}) = prctile(irf.m_shock.(fields{i}),[1 5 10 25 50 75 90 95 99]);
    pctl.y_shock.(fields{i}) = prctile(irf.y_shock.(fields{i}),[1 5 10 25 50 75 90 95 99]);
end

%%
% figure(1);clf;
% plot(path_s_shock.f'); hold on;
% plot(path_s_baseline.f','LineWidth',2,'Color',[0 0 1]); hold on; box off; ylim([0.9 2.1]);
% 
% %%
% figure(2); clf;
% plot(irf.f_shock.f','color',[0.9 0.9 0.9]); hold on;
% plot(prctile(irf.f_shock.f,1));
% plot(prctile(irf.f_shock.f,99),'--');
%%

irf.path_s_baseline = path_s_baseline; %will get passed to count_s_deviations to try to determine source of Algorithm 6 problem with reproducing inputted path of s

% Compute prob of bad regime
fields = fieldnames(path_s_shock);

for i=1:length(fields)
    
    pctl.p_bad_path_shock.(fields{i}) = sum(path_s_shock.(fields{i})==2)/(niterU*niterEPS);
    pctl.p_bad_path_shockbase.(fields{i}) = sum(path_s_baseline.(fields{i})==2)/(niterU*niterEPS);
    
    pctl.s_shock.(fields{i}) = prctile(path_s_shock.(fields{i}),[1 5 10 25 50 75 90 95 99]);
    
    
end

end