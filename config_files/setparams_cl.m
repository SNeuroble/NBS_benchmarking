%% User-defined parameters for running NBS via command line
% Note: can define all numerical arguments as numerical workspace variables
% or strings
% E.g., edge_groups_file can be defined as 
%     edge_groups_file='~/NBS_benchmarking/NBS_addon/SchizophreniaExample/Example_74node_map.mat';
%     OR
%     edge_groups_file=edge_groups;
% Strings are permissible because the GUI-based NBS parses string data
% entered by the user

% Scripts
nbs_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2/'; % NBS toolbox
% nbs_addon_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/NBS_addon/';

% Data
data_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2/SchizophreniaExample/matrices.mat';

% Model
design_matrix_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2/SchizophreniaExample/designMatrix.mat'; % 2D design matrix
contrast=[-1,1];
exchange=[];

% Edge groups - REQUIRED FOR CNBS
edge_groups_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/NBS_addon/SchizophreniaExample/Example_74node_map.mat'; % n_nodes x n_nodes edge matrix mask with nonzeros as follows: 1=subnetwork 1, 2=subnetwork 2, etc.

% Parameters
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_test_stat='t-test'; % 't-test' | 'one-sample' | 'F-test'
n_perms=1000; % recommend n_perms=5000 to appreciable reduce uncertainty of p-value estimation (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise/Theory
tthresh_first_level=3.1; % corresponds with p=0.005-0.001 (DOF=10-1000)
pthresh_second_level=0.05;
cluster_stat_type='Constrained'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA'
cluster_size_type='Extent'; % REQUIRED FOR STAT_TYPE=SIZE - 'Intensity' | 'Extent'


%%%%% DEVELOPERS ONLY %%%%%
% Use a small subset of perms for faster development - inappropriate for inference

testing=0;
test_n_perms='20';
