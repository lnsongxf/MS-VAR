function [irf, pctl] = fAlgorithm9_exante(init_in, shocks_in, param_in, opt_in, ps_bad_in,sd_base,sd_shock,varargin)


% Number of iteratios for regime shocks
niterU = opt_in.nDrawsU;

% Number of iterations for structural shocks
niterEPS = opt_in.nDrawsEPS;

% Map regime shocks
if isempty(varargin)
    u_shocks = rand(niterU,opt_in.hh);
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

s_mat = NaN(niterU,1);
for i=1:1
    ps_bad_period = ps_bad_in(i);
    bad_count = round(ps_bad_period*niterU); %number of elements in column i that should be state 2; round to nearest integer
    good_vec = ones(niterU-bad_count,1); %vector with correct number of "good" elements
    bad_vec = 2*ones(bad_count,1); %vector with correct number of "bad" elements
    full_vec = [good_vec ; bad_vec];
    s_mat(:,i) = full_vec; %replace correct column in s_mat with full concatenated vector
end

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
count = 1;
tic;
for j = 1:niterEPS
    for i = 1:niterU % counter
%         ttt = (j-1)*niterU+i;
        if floor(100*(count)/(total_l))==100*(count)/(total_l)
            disp(['percentage completed:' num2str(100*(count)/(total_l)) '%']);
        end
        % Recover path for u from s for first period
        [u_shocks_temp, ~, ~] = fGetUFromS(init_in,shocks.baseline.baseline(:,:,j),s_mat(i,:),opt_in,param_in);     
        % Pass regime shocks
        u_iter_shocks(1,1) = u_shocks_temp(1,1);
        u_iter_shocks(1,2:opt_in.hh) = u_shocks(i,2:opt_in.hh);
            
        %Get "shocked baseline" paths using initial shocks
        [Y_out_f_baseshock,  s_out_f_baseshock,p12_out_f_baseshock,p21_out_f_baseshock]   = fSimulateIterated(init_in,shocks.f.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_m_baseshock,  s_out_m_baseshock,p12_out_m_baseshock,p21_out_m_baseshock]   = fSimulateIterated(init_in,shocks.m.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
        [Y_out_y_baseshock,  s_out_y_baseshock,p12_out_y_baseshock,p21_out_y_baseshock]   = fSimulateIterated(init_in,shocks.y.baseline(:,:,j),param_in,opt_in,u_iter_shocks);
           
        %Actual shocked paths
    
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
% plot(path_s_shock.y'); hold on;
% plot(path_s_baseline','LineWidth',2,'Color',[0 0 1]); hold on; box off; ylim([0.9 2.1]);
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
    
    
    
end
toc;
end