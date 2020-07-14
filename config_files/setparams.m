%% User-defined parameters for running NBS benchmarking

% Scripts
system_dependent_paths=1; % if true, will update paths based on system using setpaths.m, replacing the paths defined here
nbs_dir='/Users/steph/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
other_scripts_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/support_scripts/'; % misc scripts used for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas

% Data parameters
%   Will be used in setup_benchmarking.m to parse files in the form:
%     task_IDs_file=[data_dir,task,subIDs_suffix];
%     non_task_IDs_file=[data_dir,non_task,subIDs_suffix];
%     this_file_task = [data_dir,task,'/',<ID>,'_',task,data_type_suffix];
%     this_file_non_task = [data_dir,non_task,'/',<ID>,'_',non_task,data_type_suffix];
%   Each data text file is assumed to be n_nodes x n_nodes x n_subjects
data_dir='/Users/steph/Documents/data/mnt/hcp_1200/matrices/';
output_dir='/Users/steph/Documents/data/mnt/NBS_benchmarking_results/';
do_TPR=1;
task='SOCIAL'; % for TPR
task_gt='SOCIAL'; % for ground truth 
non_task='REST'; % for FPR or TPR contrast
subIDs_suffix='_subIDs.txt'; 
data_type_suffix='_GSR_matrix.txt';

% Resampling parameters
n_workers=8; % num parallel workers for parfor, best to use # workers = # cores
mapping_category='subnetwork'; % for cNBS
n_repetitions=500;
n_subs_subset=40; % size of subset is full group size (N=n*2 for two sample t-test or N=n for one-sample)

% NBS parameters
% TODO: right now dmat/contrast only designed for t-test, and edge_groups can only be Shen atlas
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_test_stat='t-test'; % alternatives are one-sample and F-test - don't change for now, bc dmat and contrast based on this
n_perms='1000'; % previously: '5000'
tthresh_first_level='3.1'; % corresponds with p=0.005-0.001 (DOF=10-1000)
pthresh_second_level='0.05';
all_cluster_stat_types={'Size','TFCE','Constrained'}; % NBS stats to be benchmarked: {'Size', 'TFCE', 'Constrained', 'SEA'}
%cluster_stat_type='Constrained'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA' % smn - commented out bc looping in script
cluster_size_type='Extent'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'


%%%%% DEVELOPERS ONLY %%%%%
% Use a small subset of perms for faster development - inappropriate for inference

testing=0;
test_n_perms='20';
