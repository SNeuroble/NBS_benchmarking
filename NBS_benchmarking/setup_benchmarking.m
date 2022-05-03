%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Setup for running NBS benchmarking
% This will load data and set up the parameters needed to run NBS benchmarking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Set paths and params
setparams_bench;
if system_dependent_paths
    setpaths;
end
addpath(genpath(nbs_dir));
addpath(genpath(other_scripts_dir));

% Developers: parameter changes
if testing
    n_perms=test_n_perms;
    n_repetitions=test_n_repetitions;
    n_workers=test_n_workers;
end

% set ground truth parameters
if exist('do_ground_truth')
    if do_ground_truth
        task1=task_gt;
        do_TPR=1;
        %use_both_tasks=1;
        cluster_stat_type=stat_type_gt;
        %cluster_stat_type='Size';
        omnibus_type=omnibus_type_gt;
        n_perms=n_perms_gt;
        %n_perms='1';
      
        fprintf(['Setting params for ground truth (',cluster_stat_type,' + ',n_perms,' perm (1 perm is just to estimate t-stats)).  Will preload data from ',task1,'.\n']);
        if ~use_both_tasks % temporary check that both tasks are being used
            warning('Not performing a contrast between two tasks. Specify ''use_both_tasks=1'' to match results from Levels of Inference study.')
        end
    end
else
    do_ground_truth=0;
end

% special setup for nonparametric FDR (part of classic NBS toolbox) and parametric FDR and Bonferroni (newly added here)  since can't run it using cluster_stat_type
if strcmp(cluster_stat_type,'FDR')
    nbs_method='Run FDR';
elseif strcmp(cluster_stat_type,'Parametric_FDR') || strcmp(cluster_stat_type,'Parametric_Bonferroni')
    nbs_method='Run Parametric Edge-Level Correction';
end

%% Check expected inclusion of resting tasks
if use_both_tasks
    if do_TPR
        if ~contains(task2,'REST') 
        %if ~(strcmp(task2,'REST') || strcmp(task2,'REST2'))
            warning('Not using REST or REST2 for contrasting with task; this differs from the intended use of these benchmarking scripts.');
        else
            if use_trimmed_rest
                task2=[task2,'_',num2str(n_frames.(task1)),'frames'];
                if do_TPR
                    warning('Using trimmed rest for power calculation.');
                end
            end
        end
    else
        if ~(contains(task1,'REST') && contains(task2,'REST'))
        %if ~((strcmp(task1,'REST') && strcmp(task2,'REST2')) || (strcmp(task1,'REST2') && strcmp(task2,'REST')))
            error('Not using REST vs. REST2 or REST2 vs. REST for estimating false positives; this differs from the intended use of these benchmarking scripts.');
        else
            if use_trimmed_rest
                task1=[task1,'_',num2str(n_frames.(task1)),'frames'];
                task2=[task2,'_',num2str(n_frames.(task2)),'frames'];
            end
        end
    end
end
    
%% Get subjects IDs, num subjects, num nodes

% get subject task1 IDs
task1_IDs_file=[data_dir,task1,subIDs_suffix];
task1_IDs=fileread(task1_IDs_file);
task1_IDs=strsplit(task1_IDs,newline);
task1_IDs(strcmp('',task1_IDs)) = [];

% if do_TPR
if use_both_tasks % doing a paired sample test

    % compare IDs (thanks https://www.mathworks.com/matlabcentral/answers/358722-how-to-compare-words-from-two-text-files-and-get-output-as-number-of-matching-words)

    % get task2 IDs
    task2_IDs_file=[data_dir,task2,subIDs_suffix];
    task2_IDs=fileread(task2_IDs_file);
    task2_IDs=strsplit(task2_IDs,newline);
    task2_IDs(strcmp('',task2_IDs)) = [];
    
    % compare task1 w task2
    fprintf(['Comparing subject IDs from from task1:\n(',task1_IDs_file,')\nwith IDs from task2:\n(',task2_IDs_file,').\n']);
    [subIDs,~,~]=intersect(task1_IDs,task2_IDs);
    subIDs(strcmp('',subIDs)) = [];
    %subIDs=subIDs(2:end); % in case the first find is empty - TODO add a check here first - NEW should be fine now bc checking for empty in task1_IDs and task2_IDs AND explicitly checking for empty in intersection

else % doing a one-sample test
    fprintf(['Using subject IDs from from task1 (',task1_IDs_file,').\n']);
    task1_IDs(strcmp('',task1_IDs)) = [];
    subIDs=task1_IDs;
end

% Developers: if testing limit to subset of data
if testing
    n_subs=n_subs_subset;
else
    n_subs=length(subIDs);
end

% if doing ground truth, use all data
if do_ground_truth
    n_subs_subset=n_subs;
end
% return; %TODO: TEMPORARY

% Load one example matrix to get number of nodes and nets
template_file=[data_dir,task1,'/',subIDs{1},'_',task1,data_type_suffix];
fprintf('Template file is: %s.\n',template_file);
template=importdata(template_file);

n_var=size(template,1);
if ~use_preaveraged_constrained
    n_nodes=n_var; % assuming square
    triumask=triu(true(size(template)),1);

    template_net=summarize_matrix_by_atlas(template(:,:,1),'suppressimg',1);
    n_node_nets=size(template_net,1); % square
    trilmask_net=tril(true(n_node_nets)); % I know it's annoying to switch to tril here instead of triu, but summarize_matrix_by_atlas gives results in lower triangle (tril is my habit) while NBS uses triu...
end

%% GROUND TRUTH ONLY: Pre-load and reorder data unless already specified to be loaded
% note that for benchmarking (not ground truth), we load/reorder data during each repetition

if do_ground_truth
    
    load_data='y';
    if exist('m','var')
        load_data=input('Some data is already loaded in the workspace. Replace? (y/n)','s');
    end

    if strcmp(load_data,'y')
    
        fprintf('Loading %d subjects. Progress:\n',n_subs);

        % load data differently for paired or one-sample test
        if use_both_tasks % for paired
            
            if paired_design
                if use_preaveraged_constrained % no need to reorder/pool
                    m=zeros(n_var,n_subs*2);
                    for i = 1:n_subs
                        % task 1
                        this_file_task1 = [data_dir,task1,'/',subIDs{i},'_',task1,data_type_suffix];
                        m(:,i)=importdata(this_file_task1);
                        % task 2
                        this_file_task2 = [data_dir,task2,'/',subIDs{i},'_',task2,data_type_suffix];
                        m(:,n_subs+i)=importdata(this_file_task2);
                        
                        % print every 50 subs x 2 tasks
                        if mod(i,50)==0; fprintf('%d/%d  (x2 tasks)\n',i,n_subs); end
                    end
                else % go ahead and reorder by atlas and pool by nets/all for ground truth
                    m=zeros(n_nodes*(n_nodes-1)/2,n_subs*2);
                    m_net=zeros(n_node_nets*(n_node_nets-1)/2+n_node_nets,n_subs*2);
                    m_pool_all=zeros(1,n_subs*2);
                    for i = 1:n_subs
                        
                        % task 1
                        this_file_task1 = [data_dir,task1,'/',subIDs{i},'_',task1,data_type_suffix];
                        d=importdata(this_file_task1);
                        d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                        m(:,i) = d(triumask);

                        d_net=summarize_matrix_by_atlas(d,'suppressimg',1);
                        m_net(:,i)=d_net(trilmask_net);
                      
                        m_pool_all(i)=mean(d(triumask));

                        % task 2
                        this_file_task2 = [data_dir,task2,'/',subIDs{i},'_',task2,data_type_suffix];
                        d=importdata(this_file_task2);
                        d=reorder_matrix_by_atlas(d,mapping_category);
                        m(:,n_subs+i) = d(triumask);
                        
                        d_net=summarize_matrix_by_atlas(d,'suppressimg',1);
                        m_net(:,n_subs+i)=d_net(trilmask_net);
                        
                        m_pool_all(n_subs+i)=mean(d(triumask));

                        % print every 50 subs x 2 tasks
                        if mod(i,50)==0; fprintf('%d/%d  (x2 tasks)\n',i,n_subs); end
                    end
                end
            else
                error('This script hasn''t been fully updated/tested for two-sample yet.');
            end
        
        else % for single task v 0

            if use_preaveraged_constrained % no need to reorder/pool
                m=zeros(n_var,n_subs);
                for i = 1:n_subs
                    this_file_task1 = [data_dir,task1,'/',subIDs{i},'_',task1,data_type_suffix];
                    m(:,i)=importdata(this_file_task1);
                    % print every 100
                    if mod(i,100)==0; fprintf('%d/%d\n',i,n_subs); end
                end
            else % go ahead and reorder by atlas and pool by nets/all for ground truth
                m=zeros(n_nodes*(n_nodes-1)/2,n_subs);
                m_net=zeros(n_node_nets*(n_node_nets-1)/2+n_node_nets,n_subs);
                m_pool_all=zeros(1,n_subs);
                for i = 1:n_subs
                    this_file_task1 = [data_dir,task1,'/',subIDs{i},'_',task1,data_type_suffix];
                    d=importdata(this_file_task1);
                    d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                    m(:,i) = d(triumask);
                    d_net=summarize_matrix_by_atlas(d,'suppressimg',1);
                    m_net(:,i)=d_net(trilmask_net);
                    m_pool_all(i)=mean(d(triumask));
                    % print every 100
                    if mod(i,100)==0; fprintf('%d/%d\n',i,n_subs); end
                end
            end 
        end

    else
        fprintf('Okay, keeping previous data and assuming already reordered.\n');
    end

end


%% Set up design matrix, contrasts, and (for cNBS) edge groups

if use_both_tasks
    
    if paired_design
        
        % set up design matrix for paired one-sample t-test
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
        
    else % two-sample test
            
            error('This script hasn''t been fully updated/tested for two-sample yet.');
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

else
    
    % set up design matrix for one-sample t-test
    % data should be organized: s1, ..., sn
    dmat=ones(n_subs_subset,1);

    % set up contrasts - positive and negative
    nbs_contrast=[1];
    nbs_contrast_neg=[-1];

    % set up exchange
    nbs_exchange='';
end


% make edge groupings (for cNBS/SEA)
if ~use_preaveraged_constrained % no need to reorder/pool
    edge_groups=load_atlas_edge_groups(n_nodes,mapping_category);
    edge_groups=tril(edge_groups,-1);
    % TODO: in NBS function, should we require zero diag? Automatically clear diag? Something else?
else
    edge_groups=1:length(template)+1; % plus 1 because the script expects there to be a "0" (and will subsequently ignore..."
end


%% Assign params to structures
% Goal: should be able to run config file, load rep_params and UI from reference, and replicate reference results

% assign repetition parameters to rep_params
rep_params.data_dir=data_dir;
rep_params.testing=testing;
%rep_params.do_simulated_effect=do_simulated_effect;
%rep_params.networks_with_effects=networks_with_effects;
rep_params.mapping_category=mapping_category;
rep_params.n_repetitions=n_repetitions;
rep_params.n_subs_subset=n_subs_subset;
rep_params.do_TPR=do_TPR;
rep_params.use_both_tasks=use_both_tasks;
rep_params.paired_design=paired_design;
rep_params.task1=task1;
if use_both_tasks; rep_params.task2=task2; end

% assign NBS parameters to UI (see NBS.m)
UI.method.ui=nbs_method;
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms;
UI.thresh.ui=tthresh_first_level;
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; 
UI.size.ui=cluster_size_type;
UI.omnibus_type.ui=omnibus_type; 
UI.edge_groups.ui=edge_groups;
UI.use_preaveraged_constrained.ui=edge_groups;
UI.exchange.ui=nbs_exchange;
%UI.do_Constrained_FWER_second_level.ui=do_Constrained_FWER_second_level; 


