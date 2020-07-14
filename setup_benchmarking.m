%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Setup for running NBS benchmarking
% This will load data and set up the parameters needed to run NBS benchmarking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% set paths and params
setparams_bench;
if system_dependent_paths
    setpaths;
end
addpath(genpath(nbs_dir));
addpath(genpath(other_scripts_dir));

% Developers: parameter changes
if testing
    n_perms=test_n_perms;
end

% get task based on whether ground truth
if exist('do_ground_truth')
    if do_ground_truth
        task=task_gt;
        preload_data=1;
        do_TPR=1;
        cluster_stat_type='Size';
        n_perms='1';

        fprintf('Setting params for ground truth.\nDoing standard NBS for %s vs REST (1 perm) - will preload data.\n',task);
    end
else
    do_ground_truth=0;
    preload_data=0;
end


% get non-task IDs
non_task_IDs_file=[data_dir,non_task,subIDs_suffix];
non_task_IDs=fileread(non_task_IDs_file);
non_task_IDs=strsplit(non_task_IDs,newline);

% get subject IDs
if do_TPR

    % compare IDs (thanks https://www.mathworks.com/matlabcentral/answers/358722-how-to-compare-words-from-two-text-files-and-get-output-as-number-of-matching-words)

    % get task IDs
    task_IDs_file=[data_dir,task,subIDs_suffix];
    task_IDs=fileread(task_IDs_file);
    task_IDs=strsplit(task_IDs,newline);

    % compare w non-task
    fprintf(['Comparing subject IDs from from task (',task_IDs_file,') with IDs from non-task (',non_task_IDs_file,').\n']);
    [subIDs,~,~]=intersect(task_IDs,non_task_IDs);
    subIDs=subIDs(2:end); % bc the first find is empty - TODO add a check here first

else
    subIDs=non_task_IDs;
end

% Developers: if testing limit to subset of data
if testing
    n_subs=n_subs_subset;
else
    n_subs=length(subIDs);
end

% use all data is doing ground truth
if do_ground_truth
    n_subs_subset=n_subs;
end

% Setup up template for data

template_file=[data_dir,non_task,'/',subIDs{1},'_',non_task,data_type_suffix];
fprintf('Template file is: %s.\n',template_file);
template=importdata(template_file);
n_nodes=size(template,1); % assuming square
trimask=logical(triu(ones(size(template)),1));

% Load data

if preload_data
   
    load_data='y';
    if exist('m','var')
        load_data=input('Some data is already loaded in the workspace. Replace? (y/n)','s');
    end

    if strcmp(load_data,'y')
    
        fprintf('Loading %d subjects. Progress:\n',n_subs);

        % load data differently for TPR or FPR
        if do_TPR

            m=zeros(n_nodes*(n_nodes-1)/2,n_subs*2);
            for i = 1:n_subs
                this_file_task = [data_dir,task,'/',subIDs{i},'_',task,data_type_suffix];
                d=importdata(this_file_task);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m(:,i) = d(trimask);
                this_file_non_task = [data_dir,non_task,'/',subIDs{i},'_',non_task,data_type_suffix];
                d=importdata(this_file_non_task);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m(:,n_subs+i) = d(trimask);
                
                % print every 50 subs x 2 tasks
                if mod(i,50)==0; fprintf('%d/%d  (x2 tasks)\n',i,n_subs); end
            end
        
        else % for FPR

            m=zeros(size(template,1),size(template,2),n_subs);
            for i = 1:n_subs
                this_file_non_task = [data_dir,non_task,'/',subIDs{i},'_',non_task,data_type_suffix];
                d=importdata(this_file_non_task);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m(:,i) = d(trimask);
                % print every 100
                if mod(i,100)==0; fprintf('%d/%d\n',i,n_subs); end
            end
        
        end

    else
        fprintf('Okay, keeping previous data and assuming already reordered.\n');
    end

else
    fprintf('This version does not pre-load data.\n');
end

% make design matrix
if do_TPR
    % set up design matrix for one-sample t-test
    % data should be organized: s1_gr1,s2_gr1, ... , sn-1_group2, sn_group2
    dmat=zeros(n_subs_subset*2,n_subs_subset+1);
    dmat(1:(n_subs_subset),1)=1;
    dmat((n_subs_subset+1):end,1)=-1;
    for i=1:n_subs_subset
        dmat(i,i+1)=1;
        dmat(n_subs_subset+i,i+1)=1;
    end

    % set up contrasts - positive and negative
    nbs_contrast=zeros(1,n_subs_subset+1);
    nbs_contrast(1)=1;

    nbs_contrast_neg=nbs_contrast;
    nbs_contrast_neg(1)=-1;

    % set up exchange
    nbs_exchange=[1:n_subs_subset, 1:n_subs_subset];
else
    % set up design matrix for two-sample t-test
    % data should be organized: s1_gr1, ... sn_gr1, sn+1_gr2, ... s2*n_gr2
    dmat=zeros(n_subs_subset,2);
    dmat(1:(n_subs_subset/2),1)=1;
    dmat((n_subs_subset/2+1):end,2)=1;
   
    % set up contrasts - positive and negative
    nbs_contrast=[1,-1];
    nbs_contrast_neg=[-1,1];
    
    % set up exchange
    nbs_exchange='';
end

% make edge groupings (for cNBS)
edge_groups=load_atlas_edge_groups(n_nodes,mapping_category);
edge_groups=tril(edge_groups,-1);
% TODO: in NBS function, should we require zero diag? Automatically clear diag? Something else?

% assign params to structures
% should be able to run config file, load rep_params and UI from reference, and replicate reference results

% assign repetition parameters to rep_params
rep_params.data_dir=data_dir;
rep_params.testing=testing;
%rep_params.do_simulated_effect=do_simulated_effect;
%rep_params.networks_with_effects=networks_with_effects;
rep_params.mapping_category=mapping_category;
rep_params.n_repetitions=n_repetitions;
rep_params.n_subs_subset=n_subs_subset;
rep_params.do_TPR=do_TPR;
if do_TPR; rep_params.task=task; end
rep_params.non_task=non_task;

% assign NBS parameters to UI (see NBS.m)
UI.method.ui=nbs_method; % TODO: revise to include vanilla FDR
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms; % previously: '5000'
UI.thresh.ui=tthresh_first_level; % p=0.01
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; % 'Size' | 'TFCE' | 'Constrained' | 'SEA'
UI.size.ui=cluster_size_type; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
UI.edge_groups.ui=edge_groups; % smn
UI.exchange.ui=nbs_exchange;



