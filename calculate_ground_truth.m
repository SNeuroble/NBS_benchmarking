% Estimate ground truth

%% Setup
% make sure config files in NBS_benchmarking are correct

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
setup_benchmarking;

UI.perms.ui='1';
UI.matrices.ui=m;

if ~do_TPR
    error('Doing FPR but ground truth not defined for FPR.');
end

% make special design matrix - TODO: consider replacing this part of setup
n_subs_subset=n_subs;
if do_TPR
    % set up design matrix for one-sample t-test
    % data should be organized: s1_gr1,s2_gr1, ... , sn-1_group2, sn_group2
    dmat=zeros(n_subs_subset*2,n_subs_subset+1);
    dmat(1:(n_subs_subset),1)=1;
    dmat((n_subs_subset+1):end,1)=-1;
    for i=1:n_subs_subset
        dmat(i,i+1)=1;
        dmat(n_subs_subset+i,i+1)=1;
    end

    nbs_contrast=zeros(1,n_subs_subset+1);
    nbs_contrast(1)=1;

    nbs_exchange=[1:n_subs_subset, 1:n_subs_subset];
else
    % set up design matrix for two-sample t-test
    % data should be organized: s1_gr1, ... sn_gr1, sn+1_gr2, ... s2*n_gr2
    dmat=zeros(n_subs_subset,2);
    dmat(1:(n_subs_subset/2),1)=1;
    dmat((n_subs_subset/2+1):end,2)=1;

    nbs_contrast=[1,-1];

    nbs_exchange='';
end

UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.exchange.ui=nbs_exchange;

% Run

nbs=NBSrun_smn(UI);

% Record results (not saving p-values - just interested in effect sizes)
edge_stats=nbs.NBS.edge_stats;
cluster_stats=full(nbs.NBS.cluster_stats);

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end
if do_TPR; condition_str=['_',rep_params.task_condition]; else condition_str=['_',rep_params.non_task_condition]; end

output_filename=[output_dir,'nbs_ground_truth__',UI.statistic_type.ui,size_str,condition_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
fprintf('Saving results in %s\n',output_filename)
save(output_filename,'edge_stats','cluster_stats','UI');

