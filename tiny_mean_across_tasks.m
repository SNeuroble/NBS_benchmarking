%tiny_mean_across_tasks
% Before starting locally, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;
stat_type_gt='Size_Extent';
tasks={'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};
n_tasks=length(tasks);

save_figs__gt=1;
clim=[-0.3,0.3]; % for N=40
clim_detailed=[-0.5,0.5]; % for N=40

%% Setup
% set summary prefixes
summary_output_dir_gt=[output_dir,'all_tasks_',stat_type_gt,'_summary/'];
ground_truth_summary_prefix=[summary_output_dir_gt,'nbs_ground_truth__all_tasks_',stat_type_gt,'_'];

% output files
summary_file=[ground_truth_summary_prefix,'_all_tasks.mat'];

% setup summary output dir
mkdir(summary_output_dir_gt);


%% Load

for i = 1:n_tasks
    
    task=tasks{i};
    
    switch task
        case 'EMOTION'
            date_time_str_ground_truth='03012020_1722';
        case 'GAMBLING'
            date_time_str_ground_truth='03012020_1652';
        case 'LANGUAGE'
            date_time_str_ground_truth='03012020_1704';
        case 'MOTOR'
            date_time_str_ground_truth='03012020_1717';
        case 'RELATIONAL'
            date_time_str_ground_truth='03012020_1736';
        case 'SOCIAL' %
            date_time_str_ground_truth='03012020_1733';
        case 'WM'
            date_time_str_ground_truth='03012020_1709';
        otherwise
            error('Undefined ground truth. Did you create ground truth or set the ''date_time_str_ground_truth''?')
    end

    % set files
    ground_truth_results_basename_prefix=['nbs_ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth];
    ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
    
    % load
    %% TODO: REPLACE UI LIGHT BC PROB inadvertently saved UI.matrices from previous run of calc_ground_truth
    %% THEN NEED TO REPLACE N SUBS
    if i==1
        load(ground_truth_filename,'edge_stats','cluster_stats'); %,'UI_light');
        n_nodes=size(cluster_stats,1);
    else
        load(ground_truth_filename,'edge_stats','rep_params'); %,'UI_light');
    end
    
    % t-stat -> d-coefficient - transpose because need for fitting spline
%     n_subs=size(UI_light.design.ui,2)-1;
    n_subs=rep_params.n_subs_subset;
    d=(edge_stats/sqrt(n_subs))';
    
    d_all(:,i)=d;
    
end


%% Summarize

mean_d=mean(d_all,2);
task_d_delta=d_all-repmat(mean_d,1,n_tasks);

%% Visualize

% re-create upper triangular mask
triu_msk=triu(true(n_nodes),1);
ids_triu=find(triu_msk);

% put mean back into upper triangle
mean_d_mat=zeros(n_nodes);
mean_d_mat(triu_msk)=mean_d;

% edge-level results
draw_atlas_boundaries(mean_d_mat');
colormap(bipolar([],0.1));
caxis(clim);
if save_figs__gt
    saveas(gcf,[ground_truth_summary_prefix,'_mean_esz_by_edges'],'png')
end

% network-level results
summarize_matrix_by_atlas(mean_d_mat');
colormap(bipolar([],0.1));
caxis(clim);

if save_figs__gt
    saveas(gcf,[ground_truth_summary_prefix,'_mean_esz_by_networks'],'png')
end

% task deltas
task_d_delta_mat=zeros(n_nodes);
for i=1:n_tasks
    
    task=tasks{i};
    
    % put delta back in upper tri
    task_d_delta_mat(triu_msk)=task_d_delta(:,i);
    
    % edge-level results
    draw_atlas_boundaries(task_d_delta_mat');
    colormap(bipolar([],0.1));
    caxis(clim_detailed);
    if save_figs__gt
        saveas(gcf,[ground_truth_summary_prefix,'_',task,'_esz_delta_by_edges'],'png')
    end

    % network-level results
    summarize_matrix_by_atlas(task_d_delta_mat');
    colormap(bipolar([],0.1));
    caxis(clim);

    if save_figs__gt
        saveas(gcf,[ground_truth_summary_prefix,'_',task,'_esz_delta_by_networks'],'png')
    end
end


%% Save

save(summary_file,'mean_d','d_all','task_d_delta','tasks');

    
    
