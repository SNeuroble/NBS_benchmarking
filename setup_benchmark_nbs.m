% Setup for running NBS benchmarking
% This will set up the parameters needed to run NBS benchmarking
% A true use_previous_results flag will skip this setup

% add script directories
addpath(genpath(nbs_path))
addpath(genpath(nbs_addon_path))
addpath(genpath(other_scripts_directory))

%if ~use_previous_results

% assign repetition parameters to structure
% goal: should be able to load this and UI after loading config file and get similar results as reference
rep_params.data_path=data_path;
rep_params.testing=testing;
rep_params.do_simulated_effect=do_simulated_effect;
rep_params.networks_with_effects=networks_with_effects;
rep_params.mapping_category=mapping_category;
rep_params.n_repetitions=n_repetitions;
rep_params.n_subs_subset=n_subs_subset;

% load data
load_data='y';
if exist('m','var');
    load_data=input('Data is already loaded in workspace. Load anew? (y/n)','s');
end    
if strcmp(load_data,'y')
    fprintf('Loading data.\n');
    m=struct2array(load(rep_params.data_path,'data'));
    m=reorder_matrix_by_atlas(m,rep_params.mapping_category);
elseif strcmp(load_data,'n')
    fprintf('Using previously loaded data and assuming already reordered.\n');
else; error('Input must be y or n.');
end
n_nodes=size(m,1); % assuming square
n_subs=size(m,3);

% make design matrix for 2-sample t-test (equal group size)
dmat=zeros(rep_params.n_subs_subset,2);
dmat(1:(rep_params.n_subs_subset/2),1)=1;
dmat((rep_params.n_subs_subset/2+1):end,2)=1;

% make edge groupings (for cNBS)
edge_groups=load_atlas_edge_groups(n_nodes,rep_params.mapping_category);
edge_groups=tril(edge_groups,-1);
% TODO: in NBS function, should we require zero diag? Automatically clear diag? Something else?

% assign NBS parameters to UI structure (standard specified in NBS.m)
UI.method.ui=nbs_method; % TODO: revise to include vanilla FDR
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms; % previously: '5000'
UI.thresh.ui=zthresh_first_level; % p=0.01
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; % 'Size' | 'Constrained' | 'TFCE' - smn
UI.size.ui=cluster_size_type; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
UI.edge_groups.ui=edge_groups; % smn
UI.exchange.ui='';

%end








