%% User-defined parameters for running NBS benchmarking

% Scripts
system_dependent_paths=1; % if true, will update paths based on system using setpaths.m, replacing the paths defined here
nbs_dir='/Users/steph/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
other_scripts_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/support_scripts/'; % misc scripts used for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas

% Data parameters
%   Will be used in setup_benchmarking.m to parse files in the form:
%     task1_IDs_file=[data_dir,task1,subIDs_suffix];
%     task2_IDs_file=[data_dir,task2,subIDs_suffix];
%     this_file_task1 = [data_dir,task1,'/',<ID>,'_',task1,data_type_suffix];
%     this_file_task2 = [data_dir,task2,'/',<ID>,'_',task2,data_type_suffix];
%   Each data text file is assumed to be n_nodes x n_nodes x n_subjects
data_dir='/Users/steph/Documents/data/mnt/hcp_1200/matrices/';
output_dir='/Users/steph/Documents/data/mnt/NBS_benchmarking_results/';
do_TPR=1;
use_both_tasks=1; % for a paired-sample test
paired_design=1; % for now, if using both tasks must use a paired design
task1='EMOTION'; % for TPR
task_gt='SOCIAL'; % for ground truth 
task2='REST'; % for FPR or TPR contrast
subIDs_suffix='_subIDs.txt'; 
data_type_suffix='_GSR_matrix.txt';

% Resampling parameters
n_workers=8; % num parallel workers for parfor, best to use # workers = # cores
mapping_category='subnetwork'; % for cNBS
n_repetitions=500;
n_subs_subset=40; % size of subset is full group size (N=n*2 for two sample t-test or N=n for one-sample)

% NBS parameters
% TODO: right now dmat/contrast only designed for t-test, and edge_groups can only be Shen atlas
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR 'Run FDR'
nbs_test_stat='t-test'; % alternatives are one-sample and F-test - don't change for now, bc dmat and contrast based on this
n_perms='1000'; % previously: '5000'
tthresh_first_level='3.1'; % corresponds with p=0.005-0.001 (DOF=10-1000)
pthresh_second_level='0.05';
all_cluster_stat_types={'Omnibus'}; % NBS stats to be benchmarked: {'Size', 'TFCE', 'Constrained', 'SEA', 'FDR', Omnibus} % note that FDR is EDGE-LEVEL
%all_cluster_stat_types={'Size','TFCE','Constrained','FDR'}; % NBS stats to be benchmarked: {'Size', 'TFCE', 'Constrained', 'SEA', 'FDR'} % note that FDR is EDGE-LEVEL
%cluster_stat_type='Constrained'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA' % smn - commented out bc looping in script
cluster_size_type='Extent'; % 'Intensity' | 'Extent' - only relevant if stat_type is 'Size'
all_omnibus_types={'Threshold_Both_Dir', 'Multidimensional_cNBS', 'Multidimensional_all_edges'}; % 'Threshold_Positive' | 'Threshold_Both_Dir' | 'Average_Positive' | 'Average_Both_Dir' | 'Multidimensional_cNBS' - only relevant if omnibus_type is 'Omnibus' 

% Under development:  'Between_minus_within_cNBS' | 'Multidimensional_all_edges'

%%%%% DEVELOPERS ONLY %%%%%
% Use a small subset of perms for faster development - inappropriate for inference

testing=1;
test_n_perms=10;
test_n_repetitions=10;
test_n_workers=0;

