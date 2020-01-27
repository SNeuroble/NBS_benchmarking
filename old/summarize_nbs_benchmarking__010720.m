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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load saved results
config_file='/mnt/store1/mridata2/mri_group/smn33_data/hcp/cNBS/config_files/cfg.m';
run(config_file);

previous_results=[output_dir,previous_results_filename];
if ~(exist('previous_results__already_loaded') && strcmp(previous_results__already_loaded,previous_results))
    clearvars -except use_previous_results output_dir previous_results do_visualization
    load(previous_results);
    previous_results__already_loaded=previous_results;
end
if exist('edge_stats_summary') && ~isstruct(edge_stats_summary) % check for older saved results
    clearvars edge_stats_summary cluster_stats_summary
end


%% Summarize for interpretation (note that these matrices may have different sizes, so we sum over the last dimension)

% summarize edge and cluster stats
edge_stats_summary.mean=mean(edge_stats_all,length(size(edge_stats_all)));
edge_stats_summary.std=std(edge_stats_all,0,length(size(edge_stats_all)));

cluster_stats_summary.mean=mean(cluster_stats_all,length(size(cluster_stats_all)));
cluster_stats_summary.std=std(cluster_stats_all,0,length(size(cluster_stats_all)));

% get positives
positives=+(pvals_all<str2double(UI.alpha.ui));

% before significance masking, make sure positives are in same space as cluster-level stats
if ~isequal(size(positives),size(cluster_stats_all))
    if strcmp(UI.statistic_type.ui,'Constrained')
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
if strcmp(UI.statistic_type.ui,'Constrained')
    FWER_manual=sum(+any(positives));
else
    positives_reshaped=reshape(positives,n_nodes^2,n_repetitions);
    FWER_manual=sum(+any(positives_reshaped));
end

save(previous_results,'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual','-append');

%% Visualize

if do_visualization
    
    % make sure positives are in matrix form - esp for cNBS
    if strcmp(UI.statistic_type.ui,'Constrained')

        msk=UI.edge_groups.ui;
        % Force mask to be upper triangular only if not already (Manually verified this for lower tri mask and upper tri input)
        tmp=tril(msk,-1); if any(tmp(:)); msk=triu(msk'); end
        groups=unique(msk); groups=groups(2:end);
        n_groups=length(groups);

    %     % get matrix to group mapping
    %     n_groups=length(unique(msk))-1; % remove zero
    %     group_idx_first=zeros(n_groups,2);
    %     for i=1:n_groups
    %         [group_idx_first(i,1),group_idx_first(i,2)]=find(msk==i,1);
    %     end
    % 
        % Summarize positives in matrix
    %         [pos_idx(:,1),pos_idx(:,2)]=find(positives); % group by repetition
    %         pos_idx_mat=[group_idx_first(pos_idx(:,1),:),pos_idx(:,2)];
        n_nodes=size(msk,1);
        positives_total=zeros(n_nodes);
        for i=1:n_groups
    %             cluster_stats_summary__sig(msk==i)=cluster_stats_summary__sig__by_group(i);
            positives_total(msk==i)=sum(positives(i,:)); % TBD
        end
    end
    
    visualization_scaling_factor=100000;
    summarize_matrix_by_atlas((visualization_scaling_factor*positives_total/n_repetitions)'); % TODO: confirm whether need the transpose - not needed for Size - Extent 
%     caxis([0,0.01]);
    
    
    
%     % Get significant cluster stats
%     if strcmp(UI.statistic_type.ui,'Constrained')
%         % positives reported for each network, not the full mat. Create positives in full matrix for visualization
% 
%         % Force mask to be upper triangular only if not already (Manually verified this for lower tri mask and upper tri input)
%         msk=UI.edge_groups.ui;
%         tmp=tril(msk,-1);
%         if any(tmp(:)); msk=triu(msk'); end
% 
%         % get matrix to group mapping
%         n_groups=length(unique(msk))-1; % remove zero
%         group_idx_first=zeros(n_groups,2);
%         for i=1:n_groups
%             [group_idx_first(i,1),group_idx_first(i,2)]=find(msk==i,1);
%         end
% 
%         % Summarize positives in matrix
% %         [pos_idx(:,1),pos_idx(:,2)]=find(positives); % group by repetition
% %         pos_idx_mat=[group_idx_first(pos_idx(:,1),:),pos_idx(:,2)];
%         positives_
%         positives_summary=zeros(n_nodes);
%         for i=1:n_groups
% %             cluster_stats_summary__sig(msk==i)=cluster_stats_summary__sig__by_group(i);
%             positives_summary(msk==i)=positives(i); % TBD
%         end
% 
%     else
%         pvals_all_struct=reshape(pvals_all,n_nodes,n_nodes,n_repetitions);
%         positives=+(pvals_all_struct<str2double(UI.alpha.ui));
%         positives_summary=sum(positives,3);
%         cluster_stats_sig_all=cluster_stats_all.*positives; % why weight the positives by the effect size? don't we just care about the positives?
%         cluster_stats_sig_summary=mean(cluster_stats_sig_all,3);
%     end
     
%     edge_stats_summary_mat=structure_data(edge_stats_summary,'triangleside','upper');
%     draw_atlas_boundaries((edge_stats_summary_mat));
%     edge_stats_abs_summary_mat=structure_data(edge_stats_abs_summary,'triangleside','upper');
%     summarize_matrix_by_atlas((edge_stats_abs_summary_mat)*10);
    
%     draw_atlas_boundaries(cluster_stats_summary_sig/n_repetitions);    
%     summarize_matrix_by_atlas((cluster_stats_summary_sig/n_repetitions)'); % TODO: confirm whether need the transpose - not needed for Size - Extent 

    % summarize_matrix_by_atlas((cluster_stats_summary')*10);

    % visualize cluster stats
    % figure;
    % image(nbs_new.NBS.cluster_stats*50)

    % visualize edge stats
    %  figure;
    % image(structure_data(nbs.NBS.edge_stats*50,triu(ones(268),1)));
end
