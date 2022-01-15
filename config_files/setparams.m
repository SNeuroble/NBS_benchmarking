%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% User-defined parameters for running NBS with design matrix via command line
%
% Can define all numerical arguments as numeric or string types
% (Original NBS designed to parse string data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NBS toolbox
nbs_dir='/home/smn33/scripts/NBS1.2';

% Data: n_nodes x n_nodes x n_subs
% for testing, get from NBS toolbox "SchizophreniaExample" dir
data_file='/home/smn33/scripts/NBS1.2/SchizophreniaExample/matrices.mat';

% Model: 2D matrix
% for testing, get from NBS toolbox "SchizophreniaExample" dir
design_matrix_file='/home/smn33/scripts/NBS1.2/SchizophreniaExample/designMatrix.mat';
contrast=[-1,1];
exchange=[];

% Edge groups -- required for cNBS
% file contains n_nodes x n_nodes edge matrix mask with nonzeros as follows: 1=subnetwork 1, 2=subnetwork 2, etc.
% for testing, get a map from this NBS_benchmarking toolbox "NBS_addon" dir
edge_groups_file='/home/smn33/scripts/NBS_benchmarking/NBS_addon/SchizophreniaExample/Example_74node_map.mat'; 

% NBS parameters
nbs_method='Run NBS';
nbs_test_stat='t-test';     % 't-test' | 'one-sample' | 'F-test'
                            % Current model (see above design matrix) only designed for t-test
n_perms=1000;               % recommend n_perms=5000 to appreciably reduce uncertainty of p-value estimation (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise/Theory)
tthresh_first_level=3.1;    % t=3.1 corresponds with p=0.005-0.001 (DOF=10-1000)
                            % Only used if cluster_stat_type='Size'
pthresh_second_level=0.05;  % FWER or FDR rate
cluster_stat_type='Size';    % 'Parametric_Bonferroni' | 'Parametric_FDR' | 'Size' | 'TFCE' | 'Constrained' | 'Constrained_FDR' | 'Omnibus' | 'SEA' (under construction)
cluster_size_type='Extent'; % 'Intensity' | 'Extent'
                            % Only used if cluster_stat_type='Size'
omnibus_type='Threshold_Positive';  % 'Threshold_Positive' | 'Threshold_Both_Dir' | 'Average_Positive' | 'Average_Both_Dir' | 'Multidimensional_cNBS' | 'Multidimensional_all_edges' 
                                    % Only used if cluster_stat_type='Omnibus'

%%%%% DEVELOPERS ONLY %%%%%
% Use a small subset of permutations for faster development -- inappropriate for inference

testing=0;
test_n_perms='20';

