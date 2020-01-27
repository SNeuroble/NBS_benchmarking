%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NBS-based method benchmarking (cNBS, TFCE, etc)
% saves results
% edge_stats_all: mean and sd of edge_stats_all
% cluster_stats_all: mean and sd of cluster_stats_all
% pvals_all: total # positives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

config_file='/mnt/store1/mridata2/mri_group/smn33_data/hcp/cNBS/config_files/cfg.m';
run(config_file);
setup_benchmark_nbs
%% Run repetitions

resume_from_previous=0;
if exist('this_repetition') % resuming from another process
    response=input(sprintf('%0.0d perms already completed. Resume? ',this_repetition),'s');
    if strcmp(response,'y')
        reps_completed_previously=this_repetition;
        resume_from_previous=1;
    end
end
if resume_from_previous==0 % starting from the beginning
    reps_completed_previously=0;
    FWER=0;
    FP_mat=zeros(n_nodes);
    edge_stats_all=zeros(n_nodes*(n_nodes-1)/2,rep_params.n_repetitions);
    if strcmp(UI.statistic_type.ui,'Constrained')
        cluster_stats_all=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
        pvals_all=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
    else
        cluster_stats_all=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
        pvals_all=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
    end
end

% Do NBS - note that using parfor with more than one worker requires Parallel Computing Toolbox

if rep_params.testing; fprintf('*** TESTING MODE ***\n'); end
if rep_params.do_simulated_effect; fprintf('*** SYNTHETIC EFFECT ADDED ***\n'); end
my_pool = parpool(n_workers);
parfor (this_repetition=(1+reps_completed_previously):rep_params.n_repetitions)
%parfor (this_repetition=(1+reps_completed_previously):rep_params.n_repetitions, n_workers) % this didn't work to limit workers, so explicitly set up as above
%for this_repetition=(1+reps_completed_previously):rep_params.n_repetitions
    fprintf('* Repetition %d\n',this_repetition)

    % shuffle data
    ids=randperm(n_subs,rep_params.n_subs_subset);
    m_test=m(:,:,ids);
    m_test=reorder_matrix_by_atlas(m_test,rep_params.mapping_category); 

    % simulate effects
    if rep_params.do_simulated_effect
        effect=+ismember(edge_groups,rep_params.networks_with_effects);
        effect=effect+effect';
        m_with_effect=m_test;
        m_with_effect(:,:,1:1:n_subs/2)=m_with_effect(:,:,1:n_subs/2)+effect;
        m_test=m_with_effect;
    end

    UI_new=UI;
    UI_new.matrices.ui=m_test;

    nbs=NBSrun_smn(UI_new);

    % check for any positives (if there was no ground truth effect, this goes into the FWER calculation)
    if nbs.NBS.n>0
        FWER=FWER+1;
    end

    % record everything
    edge_stats_all(:,this_repetition)=nbs.NBS.edge_stats;
    cluster_stats_all(:,:,this_repetition)=full(nbs.NBS.cluster_stats);
    pvals_all(:,this_repetition)=nbs.NBS.pval(:); % TODO: had to vectorize for TFCE... should give all outputs in same format tho 

end

%sequential (4 perms): 196.559752 seconds
%parallel (4 perms): 123.899169 seconds

%% Save
if strcmp(UI.statistic_type.ui,'Constrained')
    cluster_stats_all=squeeze(cluster_stats_all);
end

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end

mkdir(output_dir)

output_filename=[output_dir,'nbs_benchmark_results__',UI.statistic_type.ui,size_str,'_',datestr(now,'mmddyyyy_HHMM')];
save(output_filename,'edge_stats_all','cluster_stats_all','pvals_all','FWER','UI','rep_params');

% show that results are available in the workspace
previous_results_filename__already_loaded=output_filename;


