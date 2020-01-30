%% Definitions for running NBS benchmarking and summarization

% script directories (other_scripts are misc stuff used for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas)
nbs_path='/Users/steph/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
other_scripts_directory='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/general_mri';

%nbs_path='/mridata2/home2/smn33/scripts/NBS_benchmarking/NBS1.2'; % if on server
%other_scripts_directory='/mridata2/home2/smn33/scripts/matlab/myscripts/general_mri_new/general_mri/'; % if on server

% data paths - input data is assumed to be n_nodes x n_nodes x n_subjects typically HCP toy data)
data_path='/Users/steph/Steph-Lab/Misc/ConstableLab/MRRC Neuroinformatics/resources/scripts/cpm_hackathon/test_data/HCP900_rest_n50.mat';
output_dir='/Users/steph/Steph-Lab/NBS_benchmarking/results_benchmarking/';

% data_path='/mnt/smn33_data/hcp/data_01ffd_v7_3.mat'; % if mounted locally
% data_path='/mnt/store1/mridata2/mri_group/smn33_data/hcp/data_01ffd_v7_3.mat'; % if on server 
% output_dir='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/'; % if on server

% user-defined - resampling parameters (TODO: consider putting in one variable like UI)
testing=1; % developers only - speeds up analyses for troubleshooting but inappropriate for inference
n_workers=4; % num parallel workers for parfor, best to use # workers = # cores
do_simulated_effect=0;
networks_with_effects=[1,5]; % networks to add simulated effects into - only relevant if adding effect
mapping_category='subnetwork'; % for cNBS
n_repetitions=1000;
n_subs_subset=40; % size of subset is full group size (N=n*2 for two sample t-test or N=n for one-sample)

% user-defined - NBS parameters
% TODO: right now dmat/contrast only designed for t-test, and edge_groups can only be Shen atlas
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_contrast='[1,-1]'; % Do not change - right now this is the only option
nbs_test_stat='t-test'; % alternatives are one-sample and F-test
n_perms='1000'; % previously: '5000'
zthresh_first_level='3.1'; % p=0.01
pthresh_second_level='0.05';
cluster_stat_type='Size'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA' - smn
cluster_size_type='Intensity'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
nbs_exchange='';

% user-defined: visualization parameters
do_visualization=1;

% DEVELOPERS ONLY - Use a small subset of reps and perms to speed up troubleshooting but inappropriate for inference
if testing 
    n_repetitions=70;   
    n_perms='20';
end

