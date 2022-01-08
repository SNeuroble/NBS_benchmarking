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
task1='WM'; % for TPR - EMOTION GAMBLING LANGUAGE MOTOR RELATIONAL SOCIAL WM
task_gt='GAMBLING'; %'REST_176frames'; % for ground truth
task2='REST'; % for FPR or TPR contrast (REST2 for FPR)
subIDs_suffix='_subIDs.txt'; 
data_type_suffix='_GSR_matrix.txt';

%   GROUND TRUTH ONLY: Specify whether to use resting runs which have been trimmed to match each task's scan duration
%   Number of frames per task run (single encoding direction) obtained from https://protocols.humanconnectome.org/HCP/3T/imaging-protocols.html
%   Note that all scans were acquired with the same TR so scan duration can be equivalently specified across scans using number of frames
use_trimmed_rest=0; % 0 by default
n_frames.EMOTION=176;
n_frames.GAMBLING=253;
n_frames.LANGUAGE=316;
n_frames.MOTOR=284;
n_frames.RELATIONAL=232;
n_frames.SOCIAL=274;
n_frames.WM=405;
n_frames.REST=1200;
n_frames.REST2=1200;

% Resampling parameters
n_workers=15; % num parallel workers for parfor, best to use # workers = # cores
mapping_category='subnetwork'; % for cNBS
n_repetitions=500; %n_repetitions=500;
n_subs_subset=40; % n = [40, 80, 120]; size of subset is full group size (N=n*2 for two sample t-test or N=n for one-sample)

% NBS parameters
% TODO: right now dmat/contrast only designed for t-test, and edge_groups can only be Shen atlas
nbs_method='Run NBS'; 
nbs_test_stat='t-test'; % alternatives are one-sample and F-test - don't change for now, bc dmat and contrast based on this
n_perms='1000'; % default to: 1000, more conservative: 5000
tthresh_first_level='3.1'; % for Size - corresponds with p=0.005-0.001 (DOF=10-1000)
pthresh_second_level='0.05'; 
do_Constrained_FWER_second_level=0; % for Constrained: 0->FDR, 1->Bonferroni % TODO: remove here and in setup script - no longer used (check)
all_cluster_stat_types={'Parametric_Bonferroni'};
%   inferential approach to use for benchmarking: {'Size', 'TFCE', 'Constrained', 'SEA', 'Omnibus', 'Parametric_FDR', 'Parametric_Bonferroni', 'FDR'} % note that FDR (Nonparametric), Parametric FDR, and Parametric_Bonferroni are all edge-level
% TODO: rename to "stat_type"
cluster_size_type='Extent'; % 'Intensity' | 'Extent' - only relevant if stat_type is 'Size'
all_omnibus_types={'Multidimensional_cNBS'}; % 'Threshold_Positive' | 'Threshold_Both_Dir' | 'Average_Positive' | 'Average_Both_Dir' | 'Multidimensional_cNBS' - only relevant if omnibus_type is 'Omnibus' 
%all_omnibus_types={'Threshold_Both_Dir', 'Multidimensional_cNBS', 'Multidimensional_all_edges'}; % 'Threshold_Positive' | 'Threshold_Both_Dir' | 'Average_Positive' | 'Average_Both_Dir' | 'Multidimensional_cNBS' - only relevant if omnibus_type is 'Omnibus' 
omnibus_type_gt='Multidimensional_all_edges';

%%%%% DEVELOPERS ONLY %%%%%
% Use a small subset of perms for faster development - inappropriate for inference

testing=0;
test_n_perms=100;
%test_n_perms=n_perms;
test_n_repetitions=100;
%test_n_workers=4;
test_n_workers=n_workers;

