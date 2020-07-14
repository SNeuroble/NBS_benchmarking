%% User-defined parameters for running NBS benchmarking

testing=0; % developers only - speeds up analyses for troubleshooting but inappropriate for inference

% Data parameters
do_TPR=1;
task='SOCIAL'; % for TPR
task_gt='SOCIAL'; % for ground truth 
non_task='REST'; % for FPR or TPR contrast

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
%cluster_stat_type='Constrained'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA' % smn - commented out bc looping in script
cluster_size_type='Extent'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'

% DEVELOPERS ONLY - Use a small subset of reps and perms to speed up troubleshooting - inappropriate for inference
if testing 
    n_repetitions=70;   
    n_perms='20';
end
