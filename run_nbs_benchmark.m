%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NBS-based method benchmarking (cNBS, TFCE, etc)
% saves results
% edge_stats_all: mean and sd of edge_stats_all
% cluster_stats_all: mean and sd of cluster_stats_all
% pvals_all: total # positives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
[current_path,~,~]=fileparts(mfilename('fullpath'));
addpath(genpath(current_path));

config_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS/config_files/config_nbs.m';
run(config_file);
setup_benchmark_nbs;

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
    edge_stats_all=zeros(n_nodes*(n_nodes-1)/2,n_repetitions);
    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
        cluster_stats_all=zeros(length(unique(edge_groups))-1,1,n_repetitions); % minus 1 to not count "zero"
        pvals_all=zeros(length(unique(UI.edge_groups.ui))-1,n_repetitions); % minus 1 to not count "zero"
    else
        cluster_stats_all=zeros(n_nodes,n_nodes,n_repetitions); 
        pvals_all=zeros(n_nodes*n_nodes,n_repetitions);
    end
end

% Do NBS - note that using parfor with more than one worker requires Parallel Computing Toolbox

if testing; fprintf('*** TESTING MODE ***\n'); end
if do_simulated_effect; fprintf('*** SYNTHETIC EFFECT ADDED ***\n'); end
parfor (this_repetition=(1+reps_completed_previously):n_repetitions,n_streams)
% for this_repetition=(1+reps_completed_previously):n_repetitions
    fprintf('* Repetition %d\n',this_repetition)

    % shuffle data
    ids=randperm(n_subs);
    m_test=m(:,:,ids);

    % simulate effects
    if do_simulated_effect
        effect=+ismember(edge_groups,networks_with_effects);
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
        % record survivors, if any
%             for this_clust=1:nbs.NBS.n
%     %             component_mat=full(nbs.NBS.con_mat{this_clust});
%     %             % add to running summary matrix
%     %             FP_mat=FP_mat+component_mat;
%             end
        % add to running summary of FWER
    end

    % save everything
    edge_stats_all(:,this_repetition)=nbs.NBS.edge_stats;
    cluster_stats_all(:,:,this_repetition)=full(nbs.NBS.cluster_stats);
    pvals_all(:,this_repetition)=nbs.NBS.pval(:); % TODO: had to vectorize for TFCE... should give all outputs in same format tho 

end
%sequential (4 perms): 196.559752 seconds
%parallel (4 perms): 123.899169 seconds


%% Save
if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end

mkdir(output_dir)

output_filename=[output_dir,'nbs_benchmark_results__',UI.statistic_type.ui,size_str,'_',datestr(now,'mmddyyyy_HHMM')];
save(output_filename,'edge_stats_all','cluster_stats_all','pvals_all','FWER','UI');

benchmarking_results__already_loaded=1;
