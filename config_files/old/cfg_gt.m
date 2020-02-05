%% Ground Truth Definitions

% script directories
nbs_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
nbs_addon_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS';
% dir for other scripts used (structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas)
other_scripts_directory='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/general_mri';

% input data - assumes square matrix - ground truth should use HCP toy data
data_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/ConstableLab/MRRC Neuroinformatics/resources/scripts/cpm_hackathon/test_data/HCP900_rest_n50.mat';

% ground truth directory
ground_truth_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS/nbs_results/results_ground_truth';

% user-defined - resampling parameters
create_new_standard=0; % if want to build ground truth for regression testing
ground_truth_date='11202019'; % date to compare with (mmddyyyy) for regression testing
% testing=0; % does not apply to regression tests
% do_simulated_effect=0; % does not apply to regression tests
networks_with_effects=[1,5]; % networks to add simulated effects into - only relevant if adding effect
mapping_category='subnetwork'; % for cNBS
% n_repetitions=100; % does not apply to regression tests
n_workers=4; % num parallel workers for parfor, best to use # workers = # cores - does not apply to regression tests

% user-defined - NBS parameters
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_contrast='[1,-1]'; % Do not change - right now this is the only option
% dmat - TBD - right now the above contrast is the only option
nbs_test_stat='t-test'; % alternatives are one-sample and F-test
n_perms='1000'; % previously: '5000'
zthresh_first_level='3.1'; % p=0.01
pthresh_second_level='0.05';
cluster_stat_type='Size'; % 'Size' | 'Constrained' | 'TFCE' - smn
cluster_size_type='Intensity'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
% edge_groups - TBD - right now uses default Shen grouping
nbs_exchange='';
