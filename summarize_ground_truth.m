function summarize_ground_truth(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting locally, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/
% This script summarizes ground truth effects and compares across tasks
% Must have calculated ground truth for specified tasks
% Summarization: calculates ground truth t->d
% Plot: binned effect size, binned TPR
% Usage: summarize_ground_truth({'LANGUAGE','WM'});
%   Task choices: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION; MOTOR; GAMBLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables
all_tasks={'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};
compare_all_tasks_default=1;
make_figs_default=1;
save_figs_default=1;
save_log_default=1;

p = inputParser;
addOptional(p,'tasks',all_tasks);
addOptional(p,'compare_all_tasks',compare_all_tasks_default);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'save_figs',save_figs_default);
addOptional(p,'save_log',save_log_default);
parse(p,varargin{:});

tasks=p.Results.tasks;
compare_all_tasks=p.Results.compare_all_tasks;
make_figs=p.Results.make_figs;
save_figs=p.Results.save_figs;
save_log=p.Results.save_log;


%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
% data_dir=''; output_dir=''; nbs_dir=''; other_scripts_dir='';
setpaths;
stat_type='NA'; % needed for setparams_summary
setparams_summary;

save_settings_for_all.asked.log=0;
save_settings_for_all.asked.figs=0;

if save_figs && ~make_figs
    resp=input('Specified to save figs but not make them. Do you want to make and save figs [y] or no figs [n]?\n','s');
    if strcmp(resp,'y')
        make_figs=1;
    else
        save_figs=0;
    end
end

%% Task-level ground truth

n_tasks=length(tasks);

fprintf('* Summarizing individual task ground truth.\n');

for i = 1:n_tasks
    
    task=tasks{i};
    fprintf(['Summarizing ',task,'\n']);
    date_time_str=date_time_str_ground_truth.(task);
    [dcoeff,save_settings_for_all]=do_summary(task,stat_type_gt,date_time_str,output_dir,make_figs,save_figs,save_log,save_settings_for_all); % TODO
    
    if compare_all_tasks
        d_all(:,i)=dcoeff;
    end
    
end

%% Cross-task ground truth

if compare_all_tasks
    
    % if want to compare all, check that as many tasks as full number of available tasks
    if length(tasks)==length(all_tasks)
        fprintf('* Combining effect size estimates across tasks.\n');
    else
        error(sprintf('Specified to compare all tasks but specified number of tasks (%d) does not match full number of available tasks (%d). Will not compare all tasks.\n',length(tasks),length(all_tasks)));
    end
    
    date_time_str=datestr(now,'mmddyyyy');
    [~,~]=do_summary('all_tasks',stat_type_gt,date_time_str,output_dir,make_figs,save_figs,save_log,save_settings_for_all,d_all,tasks); % TODO
    
end






%% Function to summarize tasks; here, task can be "all_tasks"
function [dcoeff,save_settings_for_all]=do_summary(task,stat_type_gt,date_time_str,output_dir,make_figs,save_figs,save_log,save_settings_for_all,varargin)

%% Setup
stat_type='NA';
setparams_summary;

% set input file names
results_basename_prefix=['nbs_ground_truth__',task,'_',stat_type_gt,'_',date_time_str];
matfilename=[output_dir,results_basename_prefix,'.mat'];

% set output file names
summary_output_dir=[output_dir,task,'_',stat_type_gt,'_summary/'];
summary_prefix=[summary_output_dir,results_basename_prefix];

% create summary output dir
if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end


%% Check whether already exist
logfile=[summary_prefix,'_log.txt'];
if save_log
    if exist(logfile,'file')
        if ~save_settings_for_all.asked.log || ~save_settings_for_all.use_same.log
            
            resp=input('Log file already exists. \nOverwrite? [y/n]\n> ','s');
            if strcmp(resp,'y')
                fprintf('Replacing log.\n');
            else
                save_log=0;
                fprintf('Okay, won''t overwrite.\n');
            end
            
            if ~save_settings_for_all.asked.log
                user_response=input(sprintf('Repeat for all? [yes/no]\n> '),'s');
                if strcmp(user_response,'yes')
                    fprintf('Using this setting for all.\n');
                    save_settings_for_all.use_same.log=1;
                    save_settings_for_all.log=save_log;
                else
                    fprintf('Okay, will ask each time.\n');
                    save_settings_for_all.use_same.log=0;
                end
                save_settings_for_all.asked.log=1;
            end
            
        else
            save_log=save_settings_for_all.log;
        end
    end
end


%% Calculate effect size

% load (TODO: only need to get n_nodes only from first file loaded)

if strcmp(task,'all_tasks') % get average dcoeff across tasks
    
    % Load mean data
    dcoeff_all=varargin{1};
    tasks=varargin{2};
    
    % Calculate mean and individual task deltas
    mean_dcoeff=mean(dcoeff_all,2);
    task_dcoeff_delta=dcoeff_all-repmat(mean_dcoeff,1,length(tasks));
    
    % Check for saved file
    % not asking about "setting_for_all" bc this only happens once
    mean_all_dcoeff_filename=[summary_prefix,'_dcoeff.mat'];
    save_task_comparison=1;
    if exist(mean_all_dcoeff_filename, 'file') == 2
        user_response=input(sprintf('Cross-task summary data already exists. Overwrite? [yes/no]\n> '),'s');
        if strcmp(user_response,'yes')
            fprintf('Replacing previous summary.\n');
        else
            fprintf('Won''t replace existing summary.\n');
            save_task_comparison=0;
        end
    end
    
    % Save
    if save_task_comparison
        save(mean_all_dcoeff_filename,'dcoeff_all','mean_dcoeff','task_dcoeff_delta','tasks');
    end
    
    % visualization stuff
    dcoeff=[mean_dcoeff, task_dcoeff_delta]; % concatenate all results for visualization
    tasks=[{'mean'}, tasks];
    edge_results_suffix={'_esz_by_edges'};
    network_results_suffix={'_esz_by_edges'};
    
    % create upper triangular mask
    n_nodes=268; % TODO: automate
    triu_msk=triu(true(n_nodes),1);
    ids_triu=find(triu_msk);
    
    n_subs=NaN;
    
else % assuming a single task
    
    tasks={task};
    load(matfilename,'edge_stats','cluster_stats','rep_params'); %,'UI_light');
    n_nodes=size(cluster_stats,1);
    
    % create upper triangular mask
    triu_msk=triu(true(n_nodes),1);
    ids_triu=find(triu_msk);
    
    % t-stat -> d-coefficient - transpose applied for fitting spline
    n_subs=rep_params.n_subs_subset;
    dcoeff=(edge_stats/sqrt(n_subs))';
    
    % file naming stuff
    edge_results_suffix={'_esz_by_edges'};
    
end


%% Mean effect size within thresholds

% only take the first bc not interested in deltas here
d=dcoeff(:,1);

% mean TPR "at" (around) thresholds
thresholds=[thresh_small, thresh_med, thresh_large];
n_edges=size(d,1);
for t=1:length(thresholds)
    
    % get edges within thresholds (edges < thresh_high; edges < thresh_high & > thresh_low)
    ids_lt_thr=abs(d) <= thresholds(t);
    if t~=1
        ids_btw_thr_and_thr_below=abs(d) <= thresholds(t) & abs(d) >= thresholds(t-1);
    end
    
    % calc percent edges within dcoeff thresholds
    perc_edges_lt_thr(t)=sum(+ids_lt_thr) * 100 / n_edges;
    if t~=1
        perc_edges_btw_thr_and_thr_below(t-1)=sum(+ids_btw_thr_and_thr_below) * 100 / n_edges;
    end
end


%% Visualize

if make_figs
    
    % reset save_settings_for_all when doing mean and deltas, bc may
    % want to save these differently than individual tasks
    if strcmp(task,'all_tasks')
        save_settings_for_all.asked.figs=0;
    end
    
    for i=1:length(tasks)
        
        d=dcoeff(:,i);
        
        % special setup for 'all_tasks' - outfiles, new hist ymax
        if strcmp(task,'all_tasks')
            if i==1; special_str='_mean';
            else
                special_str=['_',tasks{i},'_delta'];
                ax_xmin=ax_xmin_delta; ax_xmax=ax_xmax_delta; ax_ymax_esz=ax_ymax_esz_delta;
            end
        else
            special_str='';
        end
        
        esz_hist_file=[summary_prefix,special_str,'_esz_hist.png'];
        results_file_by_edge=[summary_prefix,special_str,'_esz_by_edges'];
        results_file_by_network=[summary_prefix,special_str,'_esz_by_networks'];
        
        if save_figs
            % test whether any figs already saved - check histogram
            
            if exist(esz_hist_file,'file')
                if ~save_settings_for_all.asked.figs || ~save_settings_for_all.use_same.figs
                    
                    resp=input(sprintf('Ground truth figures already exist in %s. \nOverwrite? (Otherwise will plot without saving.) [y/n]\n> ',esz_hist_file),'s');
                    if strcmp(resp,'y')
                        fprintf('Replacing ground truth figures.\n');
                    else
                        save_figs=0;
                        fprintf('Okay, won''t overwrite.\n');
                    end
                    
                    if ~save_settings_for_all.asked.figs
                        user_response=input(sprintf('Repeat for all? [yes/no]\n> '),'s');
                        if strcmp(user_response,'yes')
                            fprintf('Using this setting for all.\n');
                            save_settings_for_all.use_same.figs=1;
                            save_settings_for_all.figs=save_figs;
                        else
                            fprintf('Okay, will ask each time.\n');
                            save_settings_for_all.use_same.figs=0;
                        end
                        save_settings_for_all.asked.figs=1;
                    end
                    
                else
                    save_figs=save_settings_for_all.figs;
                end
                
            end
        end
        
        
        % 1. Plot effect size histograms
        
        figure;
        bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
        h=histogram(d,bin_edges,'Normalization','probability');
        hold on;
        plot(h.BinEdges(1:end-1) + h.BinWidth/2, h.BinCounts/length(d))
        hold off;
        
        % add stuff to hist
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
        set(gca,'fontsize',fontsz)
        % highlight
        hold on
        rectangle('Position',[-thresh_large,ax_ymin,2*thresh_large,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
        hold off
        
        if save_figs
            saveas(gcf,esz_hist_file,'png')
        end
        
        
        % 2. Plot effect size spatial distributions
        % if single task, plots effect size
        % if all, plots mean effect size and deltas for each task
        
        % put stuff back into upper triangle
        d_mat=zeros(n_nodes);
        d_mat(triu_msk)=d;
        
        % edge-level results
        draw_atlas_boundaries(d_mat');
        colormap(bipolar([],0.1));
        caxis(clim);
        
        if save_figs
            saveas(gcf,results_file_by_edge,'png')
        end
        
        % network-level results
        summarize_matrix_by_atlas(d_mat');
        colormap(bipolar([],0.1));
        caxis(clim);
        
        if save_figs
            saveas(gcf,results_file_by_network,'png')
        end
    end
end

if save_log
    fprintf('Saving log in %s.\n',logfile);
    
    fid=fopen(logfile,'w');
    fprintf(fid,'Percent between d=+/-%1.1f: %1.3f\n',[thresholds; perc_edges_lt_thr]);
    fprintf(fid,'Percent between d=%1.1f and %1.1f: %1.3f\n',[thresholds(2:end); thresholds(1:end-1); perc_edges_btw_thr_and_thr_below]);
    fprintf(fid,'N=%d total subjects\n',n_subs);
    fclose(fid);
end



