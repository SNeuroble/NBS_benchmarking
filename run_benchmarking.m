%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Do NBS-based method benchmarking (cNBS, TFCE, etc)
%
% main outputs:
% edge_stats_all: mean and sd of edge_stats_all
% cluster_stats_all: mean and sd of cluster_stats_all
% pvals_all: total # positives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User-defined

%config_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/config_files/cfg.m';
config_file='/mridata2/home2/smn33/scripts/NBS_benchmarking/config_files/cfg.m'; % if server

%% Setup

% assuming current folder is NBS_benchmarking
[current_path,~,~]=fileparts(mfilename('fullpath'));
addpath(genpath(current_path));

% set up the rest from config file
run(config_file);
setup_benchmarking;

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
    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
        cluster_stats_all=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
        pvals_all=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
    else
        cluster_stats_all=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
        pvals_all=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
    end
end

% Do NBS
% using parfor which requires Parallel Computing Toolbox, but if don't have set to 1 worker

if isempty(gcp('nocreate')); my_pool = parpool(n_workers); end % set from here bc doesn't limit to the specified n streams on server
if rep_params.testing; fprintf('*** TESTING MODE ***\n'); end
if rep_params.do_simulated_effect; fprintf('*** SYNTHETIC EFFECT ADDED ***\n'); end

parfor (this_repetition=(1+reps_completed_previously):rep_params.n_repetitions)
% for this_repetition=(1+reps_completed_previously):rep_params.n_repetitions
    fprintf('* Repetition %d\n',this_repetition)

    % shuffle data
    ids=randperm(n_subs,rep_params.n_subs_subset);
    m_test=m(:,:,ids);

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

if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
   cluster_stats_all=squeeze(cluster_stats_all);
end

%% Save

mkdir(output_dir)

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end

output_filename=[output_dir,'nbs_benchmark_results__',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM')];
save(output_filename,'edge_stats_all','cluster_stats_all','pvals_all','FWER','UI','rep_params');

% show that results are available in the workspace
previous_results_filename__already_loaded=output_filename;
