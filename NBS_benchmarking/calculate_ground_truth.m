% Estimate ground truth

%% Setup
% make sure config files in NBS_benchmarking are correct
clear all;
[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
do_ground_truth=1;
setup_benchmarking; % temporary - TODO: merge

% Run
UI_light=UI;
UI.matrices.ui=m;
nbs=NBSrun_smn(UI);

% Record results (not saving p-values - just interested in effect sizes)
edge_stats=nbs.NBS.edge_stats;
cluster_stats=full(nbs.NBS.cluster_stats);

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end
if use_both_tasks; condition_str=[rep_params.task1,'_v_',rep_params.task2]; else condition_str=rep_params.task1; end

output_filename=[output_dir,'ground_truth__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
fprintf('Saving results in %s\n',output_filename)
save(output_filename,'edge_stats','cluster_stats','UI_light','rep_params');

