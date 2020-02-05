%% Definitions for running NBS benchmarking and summarization

% script directories
nbs_path='/mnt/store1/mridata2/mri_group/smn33_data/hcp/cNBS/NBS1.2';

% input data - assumes square matrix (n_nodes x n_nodes x n_subs; typically HCP toy data)
data_path='/mnt/store1/mridata2/mri_group/smn33_data/hcp/data_01ffd_v7_3.mat';

% output_dir
%output_dir='/mnt/store1/mridata2/mri_group/smn33_data/hcp/cNBS/results_benchmarking/';

% setup
addpath(genpath(nbs_path))
m=struct2array(load(data_path,'data'));
n_subs=size(m,3);
%m=reorder_matrix_by_atlas(m,rep_params.mapping_category); % TODO: reorder at end

% make design matrix for 2-sample t-test (equal group size)
dmat=zeros(n_subs,2);
dmat(1:(n_subs/2),1)=1;
dmat((n_subs/2+1):end,2)=1;

% user-defined - NBS parameters
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_contrast='[1,-1]'; % Do not change - right now this is the only option
nbs_test_stat='t-test'; % alternatives are one-sample and F-test
n_perms='5000'; % min: '5000'
zthresh_first_level='3.1'; % p=0.01
pthresh_second_level='0.05';
cluster_stat_type='Size'; % 'Size' | 'Constrained' | 'TFCE' - smn
cluster_size_type='Intensity'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'

% assign NBS parameters to UI structure (standard specified in NBS.m)
UI.method.ui=nbs_method; % TODO: revise to include vanilla FDR
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms;
UI.thresh.ui=zthresh_first_level;
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; % 'Size' | 'Constrained' | 'TFCE' - smn
UI.size.ui=cluster_size_type; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
UI.edge_groups.ui=''; % for Constrained only
UI.exchange.ui='';

nbs=NBSrun_smn(UI);

