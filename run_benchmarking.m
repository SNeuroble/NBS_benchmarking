% Do NBS-based method benchmarking (cNBS, TFCE, etc)
%
% main outputs:
% edge_stats_all: mean and sd of edge_stats_all
% cluster_stats_all: mean and sd of cluster_stats_all
% pvals_all: total # positives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
% make sure config files in NBS_benchmarking are correct

clear all
[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
setparams_bench;

for stat_id=1:length(all_cluster_stat_types)

tic
cluster_stat_type=all_cluster_stat_types{stat_id};
setup_benchmarking;

%% Initialize counters and pre-randomize data

%resume_from_previous=0;
%if exist('this_repetition') % resuming from another process
%    response=input(sprintf('%0.0d perms already completed. Resume? ',this_repetition),'s');
%    if strcmp(response,'y')
%        reps_completed_previously=this_repetition;
%        resume_from_previous=1;
%    end
%end
%if resume_from_previous==0 % starting from the beginning
    %reps_completed_previously=0;
FWER=0;
FWER_neg=0;
edge_stats_all=zeros(n_nodes*(n_nodes-1)/2,rep_params.n_repetitions);
edge_stats_all_neg=zeros(n_nodes*(n_nodes-1)/2,rep_params.n_repetitions);
if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
    cluster_stats_all=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
    cluster_stats_all_neg=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
    pvals_all=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
    pvals_all_neg=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
else
    cluster_stats_all=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
    cluster_stats_all_neg=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
    pvals_all=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
    pvals_all_neg=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
end
%end

% randomly subsample subject IDs into groups 
for r=1:rep_params.n_repetitions
    if use_both_tasks
        if paired_design
            ids=randperm(n_subs,n_subs_subset)';
            ids=[ids;ids+n_subs];
        else
            error('This script hasn''t been fully updated/tested for two-sample yet.');
            ids=randperm(n_subs,n_subs_subset*2)';
        end
    else
        ids=randperm(n_subs,n_subs_subset)';
    end
    ids_sampled(:,r)=ids;
end

% if FPR, set up random task order
% Note that we don't want to use balanced perms (cf. Southworth et al., Properties of Balanced Permutations)
if ~do_TPR
    if (use_both_tasks && paired_design) || ~use_both_tasks
        switch_task_order=randi([0,1],n_subs_subset,rep_params.n_repetitions);
    else
        error('This script hasn''t been fully updated/tested for two-sample yet.');
    end
end


%% Run NBS repetitions
% using parfor which requires Parallel Computing Toolbox, but if can't get it set to 1 worker in setparams

c = parcluster('local');
if n_workers>c.NumWorkers
    fprintf('Specified %d workers but only %d available. Setting to max available.\n',n_workers,c.NumWorkers);
     n_workers=c.NumWorkers;
end
if isempty(gcp('nocreate')); my_pool = parpool(n_workers); end % set from here bc doesn't limit to the specified n streams on server
if rep_params.testing; fprintf('*** TESTING MODE ***\n'); end
%if rep_params.do_simulated_effect; fprintf('*** SYNTHETIC EFFECT ADDED ***\n'); end
fprintf('Starting benchmarking repetitions.\n');



parfor (this_repetition=1:rep_params.n_repetitions)
%parfor (this_repetition=(1+reps_completed_previously):rep_params.n_repetitions)
%for this_repetition=(1+reps_completed_previously):rep_params.n_repetitions
    fprintf('* Repetition %d - positive contrast\n',this_repetition)

    ids_thisrep=ids_sampled(:,this_repetition);
    m_test=zeros(n_nodes*(n_nodes-1)/2,n_subs_subset*2);
    
    if use_both_tasks
        
        if paired_design
            
            for i = 1:n_subs_subset
                
                %if FPR, use the predefined task order
                if ~do_TPR
                    if switch_task_order(i,this_repetition); this_task1=task1; this_task2=task2;
                    else; this_task1=task2; this_task2=task1;
                    end
                else
                    this_task1=task1;
                    this_task2=task2;
                end
                
                this_file_task1 = [data_dir,this_task1,'/',subIDs{ids_thisrep(i)},'_',this_task1,'_GSR_matrix.txt'];
                d=importdata(this_file_task1);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m_test(:,i) = d(trimask);
                this_file_task2 = [data_dir,this_task2,'/',subIDs{ids_thisrep(i)},'_',this_task2,'_GSR_matrix.txt'];
                d=importdata(this_file_task2);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m_test(:,n_subs_subset+i) = d(trimask);
                
            end
            
        else
            error('This script hasn''t been fully updated/tested for two-sample yet.');
            
            %if FPR, use the predefined task order
            if ~do_TPR
                if switch_task_order(i,this_repetition); this_task1=task2; this_task2=task1;
                else; this_task1=task1; this_task2=task2;
                end
            else
                this_task1=task1;
                this_task2=task2;
            end
            
            for i = 1:n_subs_subset
                this_file_task1 = [data_dir,this_task1,'/',subIDs{ids_thisrep(i)},'_',this_task1,'_GSR_matrix.txt'];
                d=importdata(this_file_task1);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m_test(:,i) = d(trimask);
            end
            for i = n_subs_subset+1:n_subs_subset*2
                this_file_task2 = [data_dir,this_task2,'/',subIDs{ids_thisrep(i)},'_',this_task2,'_GSR_matrix.txt'];
                d=importdata(this_file_task2);
                d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
                m_test(:,i) = d(trimask);
            end
        end
    else
        
        %if FPR, use the predefined task order
        if ~do_TPR
            if switch_task_order(i,this_repetition); task_flipper=-1;
            else; task_flipper=1;
            end
        else
            task_flipper=1;
        end
        
        for i = 1:n_subs_subset
            this_file_task1 = [data_dir,task1,'/',subIDs{ids_thisrep(i)},'_',task1,'_GSR_matrix.txt'];
            d=importdata(this_file_task1);
            d=reorder_matrix_by_atlas(d,mapping_category); % reorder bc proximity matters for SEA and cNBS
            m_test(:,i) = task_flipper * d(trimask);
        end
    end
    %m_test=m(:,ids_sampled(:,this_repetition)); % TEST

    % simulate effects
%    if rep_params.do_simulated_effect
%        effect=+ismember(edge_groups,rep_params.networks_with_effects);
%        effect=effect+effect';
%        m_with_effect=m_test;
%        m_with_effect(:,:,1:1:n_subs/2)=m_with_effect(:,:,1:n_subs/2)+effect;
%        m_test=m_with_effect;
%    end

    UI_new=UI;
    UI_new.matrices.ui=m_test;
    

    nbs=NBSrun_smn(UI_new);
    
    % re-run with negative
    fprintf('* Repetition %d - negative contrast\n',this_repetition)
    UI_new_neg=UI_new;
    UI_new_neg.contrast.ui=nbs_contrast_neg;

    nbs_neg=NBSrun_smn(UI_new_neg);
    

    % check for any positives (if there was no ground truth effect, this goes into the FWER calculation)
    if nbs.NBS.n>0
        FWER=FWER+1;
    end
    if nbs_neg.NBS.n>0
        FWER_neg=FWER_neg+1;
    end

    % record everything
    if strcmp(cluster_stat_type,'FDR')
        edge_stats_all(:,this_repetition)=nbs.NBS.test_stat(trimask);
        pvals_all(:,this_repetition)=nbs.NBS.con_mat{1}(:); % Note: this represents significant edges, not p-values
        
        edge_stats_all_neg(:,this_repetition)=nbs_neg.NBS.test_stat(trimask);
        pvals_all_neg(:,this_repetition)=nbs_neg.NBS.con_mat{1}(:);  % Note: this represents significant edges, not p-values
        
    else
        edge_stats_all(:,this_repetition)=nbs.NBS.edge_stats;
        cluster_stats_all(:,:,this_repetition)=full(nbs.NBS.cluster_stats);
        pvals_all(:,this_repetition)=nbs.NBS.pval(:); % TODO: had to vectorize for TFCE... should give all outputs in same format tho 

        edge_stats_all_neg(:,this_repetition)=nbs_neg.NBS.edge_stats;
        cluster_stats_all_neg(:,:,this_repetition)=full(nbs_neg.NBS.cluster_stats);
        pvals_all_neg(:,this_repetition)=nbs_neg.NBS.pval(:); % TODO: same as above
    end

end

if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
   cluster_stats_all=squeeze(cluster_stats_all);
   cluster_stats_all_neg=squeeze(cluster_stats_all_neg);
end

run_time=toc;


%% Save

mkdir(output_dir)

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end

if testing; test_str='_testing'; else test_str=''; end

if use_both_tasks; condition_str=[rep_params.task1,'_v_',rep_params.task2];
else; condition_str=rep_params.task1;
end

output_filename=[output_dir,'results__',condition_str,'_',UI.statistic_type.ui,size_str,'_grsize',num2str(rep_params.n_subs_subset),test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];
fprintf('Saving results in %s\n',output_filename)
save(output_filename,'edge_stats_all','cluster_stats_all','pvals_all','FWER','edge_stats_all_neg','cluster_stats_all_neg','pvals_all_neg','FWER_neg','UI','rep_params','run_time');

% show that results are available in the workspace
previous_results_filename__already_loaded=output_filename;

end
