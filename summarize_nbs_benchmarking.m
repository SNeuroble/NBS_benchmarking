%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Summarize results from NBS-based method benchmarking
%
% Should run run_nbs_benchmark first to get results
%
% outputs:
% edge_stats_summary: mean and sd of edge_stats_all
% cluster_stats_summary: mean and sd of cluster_stats_all
% positives_total: total # positives
%
% creates visualization
%
% TODO: make into function w user-defined as inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% user-defined

benchmarking_results='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS/results_benchmarking/nbs_benchmark_results__SEA_01272020_1236.mat';
% benchmarking_results='/Users/steph/Documents/data/mnt/nbs_benchmark_results__TFCE_01082020_0824.mat';
do_visualization=1;

%% Load saved results
% TODO: this is messy; clean up

% only load if not already loaded
if ~(exist('benchmarking_results__already_loaded') && strcmp(benchmarking_results__already_loaded,benchmarking_results))
    clearvars -except benchmarking_results do_visualization
    load(benchmarking_results,'edge_stats_all','cluster_stats_all','pvals_all','UI');
    results_vars = who('-file', benchmarking_results);
    already_summarized=ismember('positives_total',results_vars);
    output_filename=benchmarking_results;
    benchmarking_results__already_loaded=benchmarking_results;
end

%% Check whether summary already exists
% TODO: this can be checked before the above to avoid loading results then
% summary data

if already_summarized
    response=input('Already summarized. Redo? (y/n)','str');
    if strcmp(response,'y')
        fprintf('Okay, redoing.');
        already_summarized=0;
    elseif ~isstruct(edge_stats_summary) % clear previously saved summaries if too old
        fprintf('Previous summary data is too old. Redoing.');
        already_summarized=0;
    else
        fprintf('Okay, keeping previous results.');
    end
end


%% Summarize for interpretation (note that these matrices may have different sizes, so we sum over the last dimension)

summary_already_loaded=0;
if already_summarized
    load(benchmarking_results,'positives','positives_total','FWER_manual')
    n_repetitions=size(positives,length(size(positives)));
    fprintf('Depending on your version, FWER is either %0.0d (FWER) or %0.4d (FWER/n_reps).',FWER_manual,FWER_manual/n_repetitions);
    summary_already_loaded=1;
else
    % summarize edge and cluster stats
    edge_stats_summary.mean=mean(edge_stats_all,length(size(edge_stats_all)));
    edge_stats_summary.std=std(edge_stats_all,0,length(size(edge_stats_all)));

    cluster_stats_summary.mean=mean(cluster_stats_all,length(size(cluster_stats_all)));
    cluster_stats_summary.std=std(cluster_stats_all,0,length(size(cluster_stats_all)));

    % get positives
    positives=+(pvals_all<str2double(UI.alpha.ui));
    n_repetitions=size(positives,length(size(positives)));

    % before significance masking, make sure positives are in same space as cluster-level stats
    if ~isequal(size(positives),size(cluster_stats_all))
        
        if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
            if ~(size(cluster_stats_summary,1)==1 || size(cluster_stats_summary,2)==1)
                % the old cNBS saves cluster_stats_all as matrix - convert to
                % vector to match positives
                msk=UI.edge_groups.ui;
                % Force mask to be upper triangular only if not already
                tmp=tril(msk,-1); if any(tmp(:)); msk=triu(msk'); end

                % Set up to summarize all
                groups=unique(msk); groups=groups(2:end);
                n_groups=length(groups);
                n_nodes=size(msk,1);

                % Set up to summarize cluster stats
                cluster_stats_summary_mat=cluster_stats_summary;
                cluster_stats_summary.mean=zeros(n_groups,1);
                cluster_stats_summary.std=zeros(n_groups,1);
                cluster_stats_all_mat=cluster_stats_all;
                cluster_stats_all=zeros(n_groups,n_repetitions);

                % Summarize cluster stats, summarize positives in matrix form
                for i=1:n_groups
                    [mat_idx_i,mat_idx_j]=find(msk==groups(i),1);
                    cluster_stats_summary.mean(i)=cluster_stats_summary_mat.mean(mat_idx_i,mat_idx_j);
                    cluster_stats_summary.std(i)=cluster_stats_summary_mat.std(mat_idx_i,mat_idx_j);
                    for j=1:n_repetitions
                        cluster_stats_all(i,j)=cluster_stats_all_mat(mat_idx_i,mat_idx_j,j);
                    end

                end
            end
        elseif numel(positives)==numel(cluster_stats_all) % TODO: this is now true
            % reshape positives to matrix to match cluster_stats_all
            n_nodes=size(cluster_stats_all,1);
            n_repetitions=size(cluster_stats_all,3);
            positives=reshape(positives,n_nodes,n_nodes,n_repetitions);    
        else
            error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
        end
    end

    % Summarize positives, and mask with cluster_stats (all and
    % significant-only)
    cluster_stats_sig_all=cluster_stats_all.*positives; % why weight the positives by the effect size? don't we just care about the positives?
    cluster_stats_sig_summary.mean=mean(cluster_stats_sig_all,length(size(cluster_stats_sig_all)));
    cluster_stats_sig_summary.std=std(cluster_stats_sig_all,0,length(size(cluster_stats_sig_all)));
    positives_total=sum(positives,length(size(positives)));
    
    % double check FWER calculation
    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
        FWER_manual=sum(+any(positives))/n_repetitions;
    else
        positives_reshaped=reshape(positives,n_nodes^2,n_repetitions);
        FWER_manual=sum(+any(positives_reshaped))/n_repetitions;
    end

    if ~summary_already_loaded
        output_filename_summary=[output_filename(1:end-4),'_summary.mat'];
        save(output_filename_summary,'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual');
    end
end

%% Visualize

if do_visualization
    
    % convert cNBS positives to matrix form for visualization
    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
        
        % Get mask and force to be upper triangular only if not already (Manually verified this for lower tri mask and upper tri input)
        msk=UI.edge_groups.ui;
        tmp=tril(msk,-1); if any(tmp(:)); msk=triu(msk'); end
        
        % Summarize positives in matrix
        groups=unique(msk); groups=groups(2:end);
        n_groups=length(groups);
        n_nodes=size(msk,1);
        positives_total_vec=positives_total;
        positives_total=zeros(n_nodes);
        for i=1:n_groups
            positives_total(msk==i)=positives_total_vec(i);
        end
    end
    
    visualization_scaling_factor=100;
    summarize_matrix_by_atlas((visualization_scaling_factor*positives_total/n_repetitions)'); % TODO: confirm whether need the transpose - not needed for Size - Extent 
end
