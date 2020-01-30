%% Benchmarking Definitions

% script directories
nbs_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
% nbs_addon_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS';
% dir for other scripts used (structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas)
other_scripts_directory='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/general_mri';

% input data - assumes square matrix (typically HCP toy data)
% TBD: full dataset
data_path='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/ConstableLab/MRRC Neuroinformatics/resources/scripts/cpm_hackathon/test_data/HCP900_rest_n50.mat';

% input previous results
use_previous_results=0; % if this is true, then will skip processing and load previously calculated results 
previous_results_filename='nbs_benchmark_results__Constrained_12132019_1500.mat';
% previous_results_filename='nbs_benchmark_results__Size_Extent_12122019_2128.mat';
% previous_results_filename='nbs_benchmark_results__Size_Intensity_12182019_1906.mat';
% previous_results_filename='nbs_benchmark_results__TFCE_12132019_1042.mat';

% output_dir
output_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS/results_benchmarking/';

% user-defined - mainly resampling parameters (TODO: consider putting in one variable like UI)
testing=1;
do_simulated_effect=0;
networks_with_effects=[1,5]; % networks to add simulated effects into - only relevant if adding effect
mapping_category='subnetwork'; % for cNBS
n_repetitions=1000;
n_streams=4; % num parallel streams for parfor, best to use # streams = # cores
do_visualization=1;

% user-defined - NBS parameters
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_contrast='[1,-1]'; % Do not change - right now this is the only option
% dmat - TBD - right now the above contrast is the only option
nbs_test_stat='t-test'; % alternatives are one-sample and F-test
n_perms='1000'; % previously: '5000'
zthresh_first_level='3.1'; % p=0.01
pthresh_second_level='0.05';
cluster_stat_type='SEA'; % 'Size' | 'TFCE' | 'SEA' | 'Constrained' - smn
cluster_size_type='Intensity'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
% edge_groups - TBD - right now uses default Shen grouping
nbs_exchange='';
