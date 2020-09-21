% Estimate ground truth

%% Setup
% make sure config files in NBS_benchmarking are correct
clear all;
[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
do_ground_truth=1;
setup_benchmarking; % temporary - TODO: merge

% calculate subject measures for network and simple sum over full connectome
m_reordered=reorder_matrix_by_atlas(m);
m_net=zeros(10,10,nsubs);
m_pool_all=zeros(1,nsubs);
for i=1:nsubs
    m_net(:,:,i)=summarize_matrix_by_atlas(m_reordered(:,:,i));
    m_pool_all(i)=mean(tril(m(:,:,i),-1));
end
clearvars m_reordered

%% Estimate statistics
% Really just using NBS to solve GLM, except in the first case where we're also interested in recording any clusters

% Estimate edge-level t-stat and any suprathreshold clusters
UI_light=UI;
UI.matrices.ui=m;
nbs=NBSrun_smn(UI);

% Estimate network-level t-stat
UI.matrices.ui=m_net;
nbs_net=NBSrun_smn(UI);

% Estimate simple connectome-wide t-stat
UI.matrices.ui=m_pool_all;
nbs_pool_all=NBSrun_smn(UI);

%% Record results (not saving p-values - just interested in effect sizes)

edge_stats=nbs.NBS.edge_stats;
cluster_stats=full(nbs.NBS.cluster_stats);

edge_stats_net=nbs_net.NBS.edge_stats;
edge_stats_pool_all=nbs_pool_all.NBS.edge_stats;

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end
if use_both_tasks; condition_str=[rep_params.task1,'_v_',rep_params.task2]; else condition_str=rep_params.task1; end

output_filename=[output_dir,'ground_truth__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
fprintf('Saving results in %s\n',output_filename)
save(output_filename,'edge_stats','edge_stats_net','edge_stats_pool_all','cluster_stats','UI_light','rep_params');

