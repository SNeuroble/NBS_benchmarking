% Estimate ground truth test-statistic for subsequent estimation of Cohen's D
% Note that significance is later estimated from t-statistics (see supplemental_effect_size_results), except for clusters, for which sig can be estimated here

%% Setup
% make sure config files in NBS_benchmarking are correct
clear all;

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
do_ground_truth=1;
stat_type_gt='Size';
n_perms_gt='1';
setup_benchmarking; 

% Calculate significant clusters? Uses full permutations, thus takes time --  only set if you have a good amount of times
calc_significant_clusters=0;
stat_type_gt__clust='Size'; % 'Size' or 'TFCE'


%% Estimate statistics
% Really just using NBS to solve GLM, except in the first case where we're also interested in recording any clusters

% Copy a few variables for saving
UI_light=UI;
rep_params_copy=rep_params;

% Estimate edge-level t-stat and any suprathreshold clusters
UI.matrices.ui=m;
nbs=NBSrun_smn(UI);

if ~use_preaveraged_constrained
    % Estimate network-level t-stat
    % note that NBS will transform this vector to a matrix while assuming no entries on diagonal. Thus the edge_stats are correct but nothing beyond that (cluster stats, pvals) 
    UI.matrices.ui=m_net;
    nbs_net=NBSrun_smn(UI);

    % Estimate whole brain-level t-stat (simple single statistic pooled over connectome) 
    UI.matrices.ui=m_pool_all;
    nbs_pool_all=NBSrun_smn(UI);
end

% Optional: estimate significant clusters
if calc_significant_clusters
    % redo setup with sufficient perms
    stat_type_gt=stat_type_gt__clust;
    n_perms_gt='5000';
    setup_benchmarking; 

    % copy a few vars for saving
    UI_light_clust=UI;
    rep_params_copy_clust=rep_params;
    
    % estimate significant clusters
    UI.matrices.ui=m;
    nbs_clust=NBSrun_smn(UI);

    % re-run with negative
    UI.contrast.ui=nbs_contrast_neg;
    nbs_clust_neg=NBSrun_smn(UI);

    % record results
    edge_stats_clust=nbs_clust.NBS.edge_stats;
    pvals_clust=nbs_clust.NBS.pval(:); % TODO: had to vectorize for TFCE... should give all outputs in same format tho 

    edge_stats_clust_neg=nbs_clust_neg.NBS.edge_stats;
    pvals_clust_neg=nbs_clust_neg.NBS.pval(:); % TODO: same as above
    cluster_stats=full(nbs_clust.NBS.cluster_stats);
    cluster_stats_neg=full(nbs_clust_neg.NBS.cluster_stats);
end


%% Record results (not saving p-values - just interested in effect sizes)

edge_stats=nbs.NBS.edge_stats;
%cluster_stats=full(nbs.NBS.cluster_stats);
if ~use_preaveraged_constrained
    edge_stats_net=nbs_net.NBS.edge_stats;
    edge_stats_pool_all=nbs_pool_all.NBS.edge_stats;
end

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end
if use_both_tasks; condition_str=[rep_params.task1,'_v_',rep_params.task2]; else condition_str=rep_params.task1; end

output_filename=[output_dir,'ground_truth__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
fprintf('Saving results in %s\n',output_filename)
if use_preaveraged_constrained
    save(output_filename,'edge_stats','UI_light','rep_params_copy');
else
    save(output_filename,'edge_stats','edge_stats_net','edge_stats_pool_all','UI_light','rep_params_copy');
end
% also save matrices
output_matrices_filename=[output_dir,'ground_truth_matrices_ONLY__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
save(output_matrices_filename,'m','UI_light','subIDs');


if calc_significant_clusters
    output_filename=[output_dir,'ground_truth_sig_clust__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
    save(output_filename,'edge_stats_clust','edge_stats_clust_neg','pvals_clust','pvals_clust_neg','cluster_stats','cluster_stats_neg','UI_light_clust','rep_params_copy_clust');
end

