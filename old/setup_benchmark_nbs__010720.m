% Setup for NBS benchmarking
% This will set up the parameters needed to run NBS benchmarking
% A true use_previous_results flag will skip this setup

% add script directories
addpath(genpath(nbs_path))
addpath(genpath(nbs_addon_path))
addpath(genpath(other_scripts_directory))

if ~use_previous_results

    if testing
        % if testing, use only a small number of perms and reps
        n_repetitions=70;
        n_perms='20';
    end

    % load data
    %data_obj=matfile(data_path); 
    %size_data=size(data_obj,'data'); % assuming square
    %n_nodes=size_data(1); % assuming square
    %n_subs=size_data(3);
    load_data='y';
    if exist('data','var') || exist('m','var')
        load_data=input('Data is already loaded in workspace. Load anew? (y/n)','s');
    end    
    if strcmp(load_data,'y')
        fprintf('Loading data');
        m=struct2array(load(data_path,'data'));
        m=reorder_matrix_by_atlas(m,mapping_category);
    elseif strcmp(load_data,'n')
        fprintf('Using previously loaded data and assuming already reordered.');
    else; error('Input must be y or n.');
    end
    n_nodes=size(m,1); % assuming square
    n_subs=size(m,3);

    % make design matrix for 2-sample t-test (equal group size)
    dmat=zeros(n_subs_subset,2);
    dmat(1:(n_subs_subset/2),1)=1;
    dmat((n_subs_subset/2+1):end,2)=1;

    % make edge groupings (for cNBS)
    edge_groups=load_atlas_edge_groups(n_nodes,mapping_category);
    edge_groups=tril(edge_groups,-1);
    % TODO: in NBS function, should we require zero diag? Automatically clear diag? Something else?

    % assign NBS parameters (see NBS.m)
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

end
