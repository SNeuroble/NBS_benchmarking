%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Summarize and visualize results from NBS-based method benchmarking
%
% Should run run_nbs_benchmark first to get results
%
% main outputs:
% edge_stats_summary: mean and sd of edge_stats_all
% cluster_stats_summary: mean and sd of cluster_stats_all
% positives_total: total # positives
%
% TODO: make into function w benchmarking_results as inputs - also want to enable optional argument to pass results existing in workspace...
% function ['edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual']=summarize_results(?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User-defined

%previous_results_filename='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/nbs_benchmark_results__SEA_02032020_1306.mat';
previous_results_filename='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/nbs_benchmark_results__Constrained_01072020_1423.mat';
%previous_results_filename='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/nbs_benchmark_results__Size_Extent_01072020_1800.mat';
%previous_results_filename='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/nbs_benchmark_results__Size_Intensity_01082020_1416.mat';
%previous_results_filename='/mnt/store1/mridata2/mri_group/smn33_data/hcp/NBS_benchmarking_results/nbs_benchmark_results__TFCE_01082020_0824.mat';
%previous_results_filename='/Users/steph/Steph-Lab/NBS_benchmarking/results_benchmarking/nbs_benchmark_results__Size_Intensity_01292020_1534.mat';

do_visualization=1;

%% Load results

% check whether NBS results file previously loaded
load_new_data=0;
output_filename=previous_results_filename;
if exist('previous_results_filename__already_loaded')
    if ~strcmp(previous_results_filename__already_loaded,previous_results_filename)
        user_response=input(sprintf('Data already loaded into the workspace does not match the data specified in the config file. Keep already loaded data or replace with results specified in config file? (keep/replace)\nPreviously loaded: %s\nFrom config file:  %s\n>',previous_results_filename__already_loaded,previous_results_filename),'s');
        if strcmp(user_response,'replace')
            fprintf('Replacing previously loaded data with data from file.\n');
            load_new_data=1; 
        else
            fprintf('Keeping loaded data.\n');
            output_filename=previous_results_filename__already_loaded;
        end
    else
        fprintf('Using data already loaded into the workspace.\n');
        output_filename=previous_results_filename__already_loaded;
    end
else
    fprintf('Loading data.\n');
    load_new_data=1;
end

if load_new_data
    clearvars -except output_dir previous_results_filename do_visualization
    load(previous_results_filename);
    previous_results_filename__already_loaded=previous_results_filename;
end
size_cluster_stats_summary=size(cluster_stats_all);

n_nodes=size_cluster_stats_summary(1);
n_repetitions=size_cluster_stats_summary(end);

% remove any NBS summary data already loaded in workspace
if exist('edge_stats_summary') && ~isstruct(edge_stats_summary)
    clearvars edge_stats_summary cluster_stats_summary
end

%% Summarize for interpretation (note that these matrices may have different sizes, so we summarize over the last dimension)

% summarize edge and cluster stats
edge_stats_summary.mean=mean(edge_stats_all,length(size(edge_stats_all)));
edge_stats_summary.std=std(edge_stats_all,0,length(size(edge_stats_all)));

cluster_stats_summary.mean=mean(cluster_stats_all,length(size(cluster_stats_all)));
cluster_stats_summary.std=std(cluster_stats_all,0,length(size(cluster_stats_all)));

% get positives
positives=+(pvals_all<str2double(UI.alpha.ui));

% before significance masking, make sure positives are in same space as cluster-level stats
if ~isequal(size(positives),size(cluster_stats_all))

    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
    
        % the old cNBS saves cluster_stats_all as matrix - convert to
        % vector to match positives
        msk=UI.edge_groups.ui;
        % Force mask to be upper triangular only if not already
        tmp=tril(msk,-1); if any(tmp(:)); msk=triu(msk'); end
        groups=unique(msk); groups=groups(2:end);
        n_groups=length(groups);
        cluster_stats_summary_mat=cluster_stats_summary;
        cluster_stats_summary.mean=zeros(n_groups,1);
        cluster_stats_summary.std=zeros(n_groups,1);
        n_repetitions=size(positives,2);
        cluster_stats_all_mat=cluster_stats_all;
        cluster_stats_all=zeros(n_groups,n_repetitions);
        for i=1:n_groups
            [mat_idx_i,mat_idx_j]=find(msk==groups(i),1);
            cluster_stats_summary.mean(i)=cluster_stats_summary_mat.mean(mat_idx_i,mat_idx_j);
            cluster_stats_summary.std(i)=cluster_stats_summary_mat.std(mat_idx_i,mat_idx_j);
            for j=1:n_repetitions
                cluster_stats_all(i,j)=cluster_stats_all_mat(mat_idx_i,mat_idx_j,j);
            end
        end

    elseif numel(positives)==numel(cluster_stats_all)
        
        % reshape positives to matrix to match cluster_stats_all
        n_nodes=size(cluster_stats_all,1);
        n_repetitions=size(cluster_stats_all,3);
        positives=reshape(positives,n_nodes,n_nodes,n_repetitions);
    
    else
        error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
    end
    
end

% summarize positives, and mask with cluster_stats (all and
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

% save stuff but first check whether already exists
output_filename_summary=[output_filename(1:end-4),'_summary.mat'];
save_summary_file=1;
if exist(output_filename_summary, 'file') == 2
    user_response=input('Summary data already exists. Overwrite? (yes/no)\n','s');
    if strcmp(user_response,'yes')
        fprintf('Replacing previous summary.\n');
    else
        fprintf('Keeping existing summary.\n');
        save_summary_file=0;
    end
end

if save_summary_file
    save(output_filename_summary,'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual');
end

   
% make sure positives are in matrix form for cNBS and SEA visualization
if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')

    msk=UI.edge_groups.ui;
    % Force mask to be upper triangular only if not already (Manually verified this for lower tri mask and upper tri input)
    tmp=tril(msk,-1); if any(tmp(:)); msk=triu(msk'); end
    groups=unique(msk); groups=groups(2:end);
    n_groups=length(groups);

    n_nodes=size(msk,1);
    positives_total=zeros(n_nodes);
    for i=1:n_groups
        positives_total(msk==i)=sum(positives(i,:)); % TBD
    end
    
    positives_total_mat=positives_total;
    save(output_filename_summary,'positives_total_mat','-append');

end

%% Visualize
if do_visualization

    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
        visualization_scaling_factor=1000;
    else
        visualization_scaling_factor=100000;
    end

    fprintf('Making visualizations.\n');
    summarize_matrix_by_atlas((visualization_scaling_factor*positives_total/n_repetitions)'); % TODO: confirm whether need the transpose - not needed for Size - Extent 
%     caxis([0,0.01]);
    
end
