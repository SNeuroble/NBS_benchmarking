% Summarize and visualize ground truth

%% User-defined

previous_results_basename='nbs_ground_truth__TFCE_01082020_0824.mat';
do_visualization=1;

%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;

% results file full path
previous_results_filename=[output_dir,previous_results_basename];

%% Load results
load(previous_results_filename);
previous_results_filename__already_loaded=previous_results_filename;

n_nodes=size(cluster_stats_summary,1);

%% Summarize for interpretation
if do_visualization;

    % histogram
    h=histogram(edge_results);

    % network-level results
    summarize_matrix_by_atlas((edge_stats/n_repetitions)');

end


