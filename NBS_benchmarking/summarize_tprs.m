function summarize_tprs(summary_type,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script summarizes results data to ultimately compare true positive
% rates across levels of inference.
% 
% specify first arg ('summary_type') as follows:
% 
%       'dcoeff':           calculate dcoefficients for each task/ground truth level
% 
%       'positives':        calculate positives for each task/stat
% 
%       'tpr':              calculate true positive rates for each task/stat
%                           relies on completion of 'dcoeff' & 'positives'   
% 
%       'visualize_gt':     visualize ground truth, comparing levels of inference (ground truth)
%                           relies on completion of 'dcoeff'
% 
%       'visualize_tpr':    visualize true positive rates, comparing levels of inference
%                           relies on completion of 'tpr'
%
% Usage 1. 'positives'
%   summarize_tprs('positives','tasks',{'SOCIAL_v_REST'},'stat_types',{'Size_Extent'},'grsize',40,'make_figs',0,'save_logs',0);
%   Task choices: SOCIAL_v_REST; WM_v_REST; GAMBLING_v_REST; RELATIONAL_v_REST; EMOTION_v_REST; MOTOR_v_REST; GAMBLING_v_REST
%
% Usage 2. 'visualize_tpr'
%   summarize_tprs('visualize_tpr','grsize',120,'save_figs',0,'save_logs',0, 'do_combined',0);
%
% Required for 'dcoefficients': ground truth test statistics (see calculate_ground_truth)
% Required for 'positives': benchmarking results
% Required for 'tpr': dcoefficients and positives
% Required for 'visualize: 'tpr' for all tasks and stat types
%
% Relies on: summary_tools.m (defined functions), setparams_summary.m (defines misc params),
%   set_datetimestr_and_files.m (defines file parts)
%
% Recommended to first run with local access to data for calculate_positives
% --this intermediate file is slow to create but can then be
% used to recreate any summaries/visualizations. Note that this step
% doesn't rely on the ground truth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET PATHS/FUNCTIONS & DEFAULTS

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarking
addpath(genpath(current_path));
setpaths;
setparams_summary;
summary_tools; % contains summary functions

%% PARSE PARAMETERS

% parse user input
% something is up with the validator when I added this required arg, so I removed (i.e.,   ,@(x)validateattributes(x,{'char','cell'},{'nonempty'})   )

p = inputParser;
% addOptional(p,'summary_type',summary_type_default);
addRequired(p,'summary_type');
addOptional(p,'tasks',all_tasks);
addOptional(p,'stat_types',all_stat_types);
addOptional(p,'grsize',grsize_default);
addOptional(p,'save_dcoeff',save_settings.defaults.save_dcoeff);
addOptional(p,'save_benchmarking_summary',save_settings.defaults.save_benchmarking_summary);
addOptional(p,'save_figs',save_settings.defaults.save_figs);
addOptional(p,'save_logs',save_settings.defaults.save_logs);
addOptional(p,'save_summarized_data',save_settings.defaults.save_summarized_data);
% addOptional(p,'combine_all_tasks',combine_all_tasks_default);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'do_combined',do_combined_default);
% addOptional(p,'plot_ground_truth_only',plot_ground_truth_only_default);
parse(p,summary_type,varargin{:});

% summary_type=p.Results.summary_type;
tasks=p.Results.tasks;
stat_types=p.Results.stat_types;
grsize=p.Results.grsize; % used in setting files and param
save_settings.do.save_dcoeff=p.Results.save_dcoeff;
save_settings.do.save_benchmarking_summary=p.Results.save_benchmarking_summary;
save_settings.do.save_figs=p.Results.save_figs;
save_settings.do.save_logs=p.Results.save_logs;
save_settings.do.save_summarized_data=p.Results.save_summarized_data;
% combine_all_tasks=p.Results.combine_all_tasks;
make_figs=p.Results.make_figs;
pp.do_combined=p.Results.do_combined;
if pp.do_combined; make_figs=0; else; make_figs=p.Results.make_figs; end
% plot_ground_truth_only=p.Results.plot_ground_truth_only;


%% SETUP

% check requirements and silence some warnings
if make_figs || save_settings.do.save_logs % Check for curve fitting toolbox, required for figs or log
    v_=ver;
    [installedToolboxes{1:length(v_)}] = deal(v_.Name);
    curve_toolbox_exists = all(ismember('Curve Fitting Toolbox',installedToolboxes));
    if ~curve_toolbox_exists
        error('Curve fitting toolbox required for making summary figures, logs, or combining tasks but not installed. Please re-run with the toolbox installed OR with options ''make_figs'' and ''ings.do.save_logs'' set to 0.');
    end
end
warning('off', 'STRUCTURE_DATA:ASSUME_LOWER_TRI');
warning('off', 'SPLINES:CHCKXYWP:NaNs');

% set up naming/indexing vars based on levels of inference and ground truth
% stat_level_map.do_overlay=1; %TODO: set as user-defined param

<<<<<<< HEAD
stat_level_map.stat_types=all_stat_types;
stat_level_map.stat_levels_str=all_stat_types;
for s=1:size(stats_levelstr_map,1)
    stat_level_map.stat_levels_str{contains(stat_level_map.stat_levels_str,stats_levelstr_map{s,1})}=stats_levelstr_map{s,2};
end
=======
        fprintf(['Summarizing TPRs - ',task,'::',stat_type,'\n'])
        setparams_summary;
        
        ground_truth_results_basename_prefix=['ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
        bench_results_basename_prefix=['results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
        
        % set results filenames
        ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
        results_filename=[output_dir,bench_results_basename_prefix,'.mat'];
        benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];
        
        % set summary prefixes
        summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
        summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
        
        % make summary output dir
        if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        
        % check whether to save, incl overwriting existing
        summarize_benchmarking=1;
        [save_settings_for_all,summarize_benchmarking] = check_whether_to_save(save_settings_for_all,summarize_benchmarking,'summarize','Summary data',benchmarking_summary_filename);
        
        
        %% GROUND TRUTH EFFECT SIZES
        
        try
            load(ground_truth_filename,'edge_stats','edge_stats_net','edge_stats_pool_all','cluster_stats','rep_params');
            [warnmsg,~] = lastwarn;
            if contains(warnmsg,'Variable ') && contains(warnmsg,'not found.')
                error(['Unable to load all necessary variables from ',ground_truth_filename,'.\nPlease check these exist and try again.']);
            end
        catch
            error('Looks like ground truth data needed for calculating TPR does not exist for %s. Try running calculate_ground_truth.m\n',task);
        end
        
        % t-stat -> d-coefficient - transpose because need for fitting spline
        n_subs_total=rep_params.n_subs_subset;
        dcoeff=(edge_stats/sqrt(n_subs_total))';
        dcoeff_net=(edge_stats_net/sqrt(n_subs_total))';
        dcoeff_omnibus=(edge_stats_pool_all/sqrt(n_subs_total))';
        
        % get num nodes, edges, networks, and network pairs
        n_edges=length(dcoeff);
        n_nodes=size(cluster_stats,1);
        n_net_pairs=length(dcoeff_net);
        n_nets=sort(roots([1 1 -2*n_net_pairs])); % assuming n_nets x n_nets, x = n*(n+1)/2 -> n^2 + n - 2x
        n_nets=n_nets(end);
        
        % make upper triangular masks
        triu_msk=triu(true(n_nodes),1);
        ids_triu=find(triu_msk);
        triu_msk_net=triu(true(n_nets));
        
        % define num features as num edges/net pairs for histogram vis
        pp.n_features{1}=n_edges;
        pp.n_features{2}=n_net_pairs; % should be 55 for shen268
        
        % convert ground truth network results from lower to upper triangle - this sad mismatch is an unfortunate consequence of my summat scripts using the lower tri but NBS using upper tri
        tmp=tril(true(n_nets));
        dcoeff_net=structure_data(dcoeff_net,'mask',tmp);
        dcoeff_net=dcoeff_net';
        dcoeff_net=dcoeff_net(triu_msk_net);
        
        
        
        %% BENCHMARKING RESULTS - TRUE POSITIVE RATES
        
        % Load and summarize benchmarking results: 'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual'
        if summarize_benchmarking
            
            load(results_filename);
            n_repetitions=rep_params.n_repetitions;
                            
            % get positives and summarize
            positives=+(pvals_all<str2double(UI.alpha.ui));
            positives_neg=+(pvals_all_neg<str2double(UI.alpha.ui));
            positives_total=sum(positives,length(size(positives)));
            positives_total_neg=sum(positives_neg,length(size(positives)));
            
            % summarize edge and cluster stats (saved but not currently used in visualization/log)
            edge_stats_summary.mean=mean(edge_stats_all,length(size(edge_stats_all)));
            edge_stats_summary.std=std(edge_stats_all,0,length(size(edge_stats_all)));
            edge_stats_summary_neg.mean=mean(edge_stats_all_neg,length(size(edge_stats_all_neg)));
            edge_stats_summary_neg.std=std(edge_stats_all_neg,0,length(size(edge_stats_all_neg)));
            
            if strcmp(UI.statistic_type.ui,'FDR')
                cluster_stats_summary.mean=0;
                cluster_stats_summary.std=0;
                cluster_stats_summary_neg.mean=0;
                cluster_stats_summary_neg.std=0;
            else
                cluster_stats_summary.mean=mean(cluster_stats_all,length(size(cluster_stats_all)));
                cluster_stats_summary.std=std(cluster_stats_all,0,length(size(cluster_stats_all)));
                cluster_stats_summary_neg.mean=mean(cluster_stats_all_neg,length(size(cluster_stats_all_neg)));
                cluster_stats_summary_neg.std=std(cluster_stats_all_neg,0,length(size(cluster_stats_all_neg)));
            end
            % REMOVED for now since not used yet: get positive statistic values at every repetition
%                 % make sure positives are in same space as cluster-level stats
%                 size_cluster_stats_all=size(cluster_stats_all);
%                 n_dim__cluster_stats_all=length(size_cluster_stats_all); %note that matrices may have different sizes, so we summarize over the last dimension)
%                 if ~isequal(size(positives),size(cluster_stats_all))
%                     if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
%                         error('Something went wrong - this shouldn''t happen anymore, only in old summaries created by old script.')
%                     elseif numel(positives)==numel(cluster_stats_all)
%                         % reshape positives to matrix to match cluster_stats_all
%                         positives=reshape(positives,n_nodes,n_nodes,n_repetitions);
%                         positives_neg=reshape(positives_neg,n_nodes,n_nodes,n_repetitions);
%                     elseif strcmp(UI.statistic_type.ui,'Omnibus')
%                         warning('This is ONLY a temporary quick fix for the new omnibus.');
%                         cluster_stats_all=squeeze(cluster_stats_all(1,1,:))';
%                     else; error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
%                     end
%                 end 
% 
%                 % mask stats at every repetition by whether positivs 
%                 cluster_stats_sig_all=cluster_stats_all.*positives;
%                 cluster_stats_sig_summary.mean=mean(cluster_stats_sig_all,n_dim__cluster_stats_all);
%                 cluster_stats_sig_summary.std=std(cluster_stats_sig_all,0,n_dim__cluster_stats_all);
%                 
%                 cluster_stats_sig_all_neg=cluster_stats_all_neg.*positives_neg;
%                 cluster_stats_sig_summary_neg.mean=mean(cluster_stats_sig_all_neg,n_dim__cluster_stats_all);
%                 cluster_stats_sig_summary_neg.std=std(cluster_stats_sig_all_neg,0,n_dim__cluster_stats_all);
            
            % double check FWER calculation
            if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA') || strcmp(UI.statistic_type.ui,'Omnibus') || strcmp(UI.statistic_type.ui,'Parametric_FDR') || strcmp(UI.statistic_type.ui,'Parametric_Bonferroni')
                FWER_manual=sum(+any(positives))/n_repetitions;
                FWER_manual_neg=sum(+any(positives_neg))/n_repetitions;
            else
                positives_reshaped=reshape(positives,n_nodes^2,n_repetitions);
                positives_reshaped_neg=reshape(positives_neg,n_nodes^2,n_repetitions);
                FWER_manual=sum(+any(positives_reshaped))/n_repetitions;
                FWER_manual_neg=sum(+any(positives_reshaped_neg))/n_repetitions;
            end
            
            n_subs_subset=rep_params.n_subs_subset;
            n_perms=UI.perms.ui;
            if exist('run_time','var')
                run_time_h=run_time/(60*60);
            else
                run_time_h=NaN;
            end
            
            save(benchmarking_summary_filename,'edge_stats_summary','edge_stats_summary_neg','cluster_stats_summary','cluster_stats_summary_neg','positives','positives_neg','positives_total','positives_total_neg','FWER_manual','FWER_manual_neg','n_repetitions','n_subs_subset','run_time_h','n_perms','-v7.3');
        else
            load(benchmarking_summary_filename,'positives_total','positives_total_neg','n_repetitions','n_subs_subset','run_time_h','n_perms')
            [warnmsg,~] = lastwarn;
            if contains(warnmsg,'Variable ') && contains(warnmsg,'not found.')
                error(['Unable to load all necessary variables from ',benchmarking_summary_filename,'.\nPlease check these exist and try again. (Recently added grsize to the filename string--check that this is included.)']);
            end
            if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA') % need for summary in edge_groups
                load(results_filename,'UI');
            end
        end
        
        % calculate TPR
        
        ids_pos_vec=dcoeff>0;
        ids_neg_vec=dcoeff<0;
        
        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
            edge_groups_triu=UI.edge_groups.ui';
            edge_groups_vec=edge_groups_triu(ids_triu);
            ids_pos=edge_groups_vec(ids_pos_vec);
            ids_neg=edge_groups_vec(ids_neg_vec);
        elseif contains(stat_type,'Omnibus')
            ids_pos=1;
            ids_neg=1;
        else
            ids_pos=ids_triu(ids_pos_vec);
            ids_neg=ids_triu(ids_neg_vec);
        end
        
        true_positives=zeros(size(dcoeff));
        if strcmp(stat_type,'Parametric_FDR') || strcmp(stat_type,'Parametric_Bonferroni')
            % already upper triangle
            true_positives(ids_pos_vec)=positives_total(ids_pos_vec);
            true_positives(ids_neg_vec)=positives_total_neg(ids_neg_vec);
        else
            true_positives(ids_pos_vec)=positives_total(ids_pos);
            true_positives(ids_neg_vec)=positives_total_neg(ids_neg);
        end
        tpr=true_positives*100/n_repetitions;
        
        
        
        %% VISUALIZATION
        
        % Setup for figures and logs
        
        % - get network-level summary of dcoeff and tpr
        dcoeff_scaled{1}=dcoeff;
        tpr_scaled{1}=tpr;
        
        ids_pos_summat=dcoeff_net>0;
        ids_neg_summat=dcoeff_net<0;
        
        % true postive counts by networks
        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
            true_positives_summat(ids_pos_summat)=positives_total(ids_pos_summat);
            true_positives_summat(ids_neg_summat)=positives_total_neg(ids_neg_summat);
        else
            true_positives_summat=summarize_matrix_by_atlas(true_positives,'suppressimg',1,'do_std',1)';
            true_positives_summat=true_positives_summat(triu_msk_net)';
            % [true_positives_summat,true_positives_summat_std]=summarize_matrix_by_atlas(true_positives,'suppressimg',1,'do_std',1)';
            % true_positives_summat_std=true_positives_summat_std(triu_msk_summat)';
        end
        
        dcoeff_scaled{2}=dcoeff_net;
        tpr_scaled{2}=(true_positives_summat*100/n_repetitions)';
        
        % - set up data for log
        log_data.n_repetitions=n_repetitions;
        log_data.n_perms=n_perms;
        log_data.n_subs_subset=n_subs_subset;
        log_data.n_subs_total=n_subs_total;
        log_data.run_time_h=run_time_h;
        
        % Do figures and logs for both levels of scaling (edge and network)
        for scaling=1:2
            
            if make_figs || save_log
                % Fit curves for TPR v effect size to all available task data - needed for both visualization and log
                [tpr_fit{scaling},res_scaled{scaling},dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},~]=...
                fit_spline(dcoeff_scaled{scaling},tpr_scaled{scaling},pp.spline_smoothing{scaling},pp.window_sz{scaling});
                       
                % Make figs and log
                if make_figs
                    save_settings_for_all=visualize_tprs(dcoeff_scaled{scaling},tpr_scaled{scaling},dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},tpr_fit{scaling},res_scaled{scaling},triu_msk,summary_prefix,pp,scaling,save_figs,save_settings_for_all);
%                     save_settings_for_all=visualize_tprs(dcoeff_scaled{scaling},tpr_scaled{scaling},dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},dcoeff_scaled{scaling},tpr_fit{scaling},res_scaled{scaling},res_mat,summary_prefix,pp,scaling,save_figs,save_settings_for_all);
                end

                if save_log
                    save_settings_for_all=write_summary_to_log(dcoeff_windowed{scaling},tpr_windowed{scaling},log_data,summary_prefix,pp,scaling,save_log,save_settings_for_all);
                end
            end
            
>>>>>>> expanded_tests

stat_level_map.stat_gt_levels=zeros(1,length(stat_level_map.stat_levels_str));
stat_level_map.stat_gt_levels_str=cell(1,length(stat_level_map.stat_levels_str));
for s=1:size(statlevel_gtlevel_map,1)
    idx=contains(stat_level_map.stat_levels_str,statlevel_gtlevel_map{s,1});
    stat_level_map.stat_gt_levels(idx)=statlevel_gtlevel_map{s,3};
    stat_level_map.stat_gt_levels_str(idx)=statlevel_gtlevel_map(s,2);
    % TODO: use this for ground truth categories
end

% TODO: double-check saving for combined tasks


%% MAIN

switch summary_type
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE GROUND TRUTH EFFECT SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'dcoeff'
fprintf('** Calculating d-coefficients **\n');

for t=1:length(tasks)
    
    task=tasks{t};
    fprintf(['Doing: ',task,'\n'])
    
    % calculate and save d-coefficients
    set_datetimestr_and_files; % set datetimestr for specified task, stat_type, data_source
    save_settings = summary_tools.calculate_dcoefficients(ground_truth_filename,stat_level_map,ground_truth_dcoeff_filename,save_settings);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POSITIVES
% note: this is the only one where individual stat_types can be specificied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'positives'
fprintf('**** Calculating positives ****\n');

for t=1:length(tasks)
    for s=1:length(stat_types)
        
        task=tasks{t};
        stat_type=stat_types{s};
        fprintf(['Doing: ',task,'::',stat_type,'\n'])
        
<<<<<<< HEAD
        % calculate and save tpr (checking before saving is built in)
        set_datetimestr_and_files; % set datetimestr for specified task, stat_type, data_source
        if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        save_settings=summary_tools.calculate_positives(results_filename,benchmarking_summary_filename,save_settings);
=======
        % set results filenames
        date_time_str_now=datestr(now,'mmddyyyy');
        bench_results_basename_prefix=['results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_now];
%         bench_results_basename_prefix=['results__',task,'_',stat_type,omnibus_str,'_','grsize',num2str(grsize),'_',date_time_str_now];
        benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];
        
        % set summary prefixes for combined tasks
        summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
        summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_now];
        
        % make summary output dir
        if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        
        % save intermediate data
        [save_settings_for_all,save_combined] = check_whether_to_save(save_settings_for_all,save_combined,'combined','Combined summary data',benchmarking_summary_filename);
        if save_combined
            save(benchmarking_summary_filename,'dcoeff_scaled_all','tpr_scaled_all','log_data_combined','-v7.3'); % TODO: change here and above *_all -> *_combined_tasks to match naming needed for "compare_methods"
        end
        
        for scaling=1:length(dcoeff_windowed)            
            
%             dcoeff_scaled{scaling}=dcoeff_scaled_all{s}{scaling}(:);
%             tpr_scaled{scaling}=tpr_scaled_all{s}{scaling}(:);

            if make_figs || save_log
                % Fit curves for TPR v effect size to all available task data - needed for both visualization and log
                [tpr_fit{scaling},res_scaled{scaling},dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},~]=...
                fit_spline(dcoeff_scaled_all{s}{scaling}(:),tpr_scaled_all{s}{scaling}(:),pp.spline_smoothing{scaling},pp.window_sz{scaling});
            end
            
            % Make figs and log
             if make_figs || make_figs__only_combined
                save_settings_for_all=visualize_tprs(dcoeff_scaled_all{s}{scaling}(:),tpr_scaled_all{s}{scaling}(:),dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},tpr_fit{scaling},res_scaled{scaling},triu_msk,summary_prefix,pp,scaling,save_figs,save_settings_for_all);
             end
            
            if save_log
                save_settings_for_all=write_summary_to_log(dcoeff_windowed{scaling},tpr_windowed{scaling},log_data_combined,summary_prefix,pp,scaling,save_log,save_settings_for_all); 
            end
            
        end
>>>>>>> expanded_tests
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE AND AGGREGATE TRUE POSITIVES
% i.e., combine effect size and positives above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'calculate_tpr'
fprintf('**** Aggregating true positives ****\n');
if length(tasks)<length(all_tasks)
    error('Not enough tasks specified to combine all.')
end

for t=1:length(tasks)
    for s=1:length(stat_types)

        task=tasks{t};
        stat_type=stat_types{s};
        fprintf(['Doing: ',task,'::',stat_type,'\n'])

        % set files for specified task, stat_type, data_source and make output
        set_datetimestr_and_files;
        if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end

        % need original bench results for summary in edge_groups - TODO: save edge_groups into summary matfile
        if strcmp(stat_level_map.stat_gt_levels_str(s),'network')
            if exist(results_filename,'file')
                load(results_filename,'UI');
                edge_groups=UI.edge_groups.ui;
            else
                warning('Loading edge groups from edge groups file, not from original results.');
                load(edge_groups_filename,'edge_groups');
            end
        else
            edge_groups=[];
        end
        
        % calculate true positives
        [dcoeff.(stat_types{s})(:,t),tpr.(stat_types{s})(:,t),log_data.(stat_types{s}).(tasks{t})]=summary_tools.calculate_tpr(benchmarking_summary_filename,ground_truth_dcoeff_filename,stat_type,stat_level_map.stat_gt_levels_str{s},pp.remove_matrix_diag,edge_groups);

    end
end
    
% save combined-task data
save_settings = summary_tools.check_whether_to_save(save_settings,'save_summarized_data','Combined summary data',combined_summary_filename); % TODO
if save_settings.do.save_summarized_data
    if ~exist(combined_summary_dir,'dir'); mkdir(combined_summary_dir); end
    save(combined_summary_filename,'dcoeff','tpr','log_data','-v7.3');
    fprintf(['Saved combined data in ',combined_summary_filename,'.\n']);
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZE GROUND TRUTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'visualize_gt'
fprintf('**** Visualizing ground truth ****\n');
if length(tasks)<length(all_tasks)
    error('Not enough tasks specified to combine for visualizations.')
end
if length(stat_types)<length(all_stat_types)
    error('Not enough statistic types specified to combine for visualizations.')
end

set_datetimestr_and_files;

if ~exist(ground_truth_vis_dir,'dir'); mkdir(ground_truth_vis_dir); end

load(combined_summary_filename,'dcoeff');

% count dimensions and make upper triangular masks
pp.n_stat_types=length(stat_types);
[~,matching_stat_idx]=unique(stat_level_map.stat_gt_levels);
pp.n_gt_levels=length(matching_stat_idx);
pp.n_tasks=length(tasks);

for g=1:pp.n_gt_levels
    
    s_idx=matching_stat_idx(g);
    gt_level_str=stat_level_map.stat_gt_levels_str{s_idx};
    
    pp.n_features.(gt_level_str)=length(dcoeff.(stat_types{s_idx})(:,1));
    pp.n_nodes.(gt_level_str)=int16(roots([1 1 -2*pp.n_features.(gt_level_str)])); % assuming n_nets x n_nets, x = n*(n+1)/2 -> n^2 + n - 2x
    pp.n_nodes.(gt_level_str)=pp.n_nodes.(gt_level_str)(end) + pp.remove_matrix_diag.(gt_level_str);
    pp.triu_msk.(gt_level_str)=triu(true(pp.n_nodes.(gt_level_str)),pp.remove_matrix_diag.(gt_level_str));
    pp.ids_triu.(gt_level_str)=find(pp.triu_msk.(gt_level_str));
        
end

summary_tools.visualize_ground_truth(dcoeff,ground_truth_vis_filename_prefix,pp,stat_level_map,save_settings);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZE INDIVIDUAL TASKS AND COMBINED TASK COMPARISON (ACROSS STATS)
% depending on whether pp.do_combined=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'visualize_tpr'
fprintf('**** Visualizing tpr results ****\n');
if length(tasks)<length(all_tasks)
    error('Not enough tasks specified to combine for visualizations.')
end
if length(stat_types)<length(all_stat_types)
    error('Not enough statistic types specified to combine for visualizations.')
end

set_datetimestr_and_files;

if pp.do_combined
    filename_prefix=combined_filename_prefix;
else
    filename_prefix=combined_by_task_filename_prefix;
    if ~exist(combined_by_task_dir,'dir'); mkdir(combined_by_task_dir); end
    if ~exist(log_dir,'dir'); mkdir(log_dir); end
end

load(combined_summary_filename,'dcoeff','tpr','log_data');

for s=1:length(stat_types)
    
    if ~strcmp(stat_level_map.stat_gt_levels_str{s},'whole_brain')
        
        pp.level=s;
        
        if pp.do_combined
            % fit spline
            [tpr_fit.(stat_types{s}),res.(stat_types{s}),dcoeff_windowed.(stat_types{s}),tpr_windowed.(stat_types{s}),tpr_std.(stat_types{s}),~]=...
                fit_spline(dcoeff.(stat_types{s})(:),tpr.(stat_types{s})(:),pp.spline_smoothing{stat_level_map.stat_gt_levels(s)},pp.window_sz{stat_level_map.stat_gt_levels(s)});
        else
            for t=1:length(tasks)
                [tpr_fit.(stat_types{s})(:,t),res.(stat_types{s})(:,t),dcoeff_windowed.(stat_types{s}).(tasks{t}),tpr_windowed.(stat_types{s}).(tasks{t}),tpr_std.(stat_types{s}).(tasks{t}),~]=...
                    fit_spline(dcoeff.(stat_types{s})(:,t),tpr.(stat_types{s})(:,t),pp.spline_smoothing{stat_level_map.stat_gt_levels(s)},pp.window_sz{stat_level_map.stat_gt_levels(s)});
            
                %log - TODO: consider a single combined log
                this_log_filename_prefix=[log_filename_prefix,'_',tasks{t},'_',stat_types{s}];
                save_settings=summary_tools.write_summary_log(dcoeff_windowed.(stat_types{s}).(tasks{t}),tpr_windowed.(stat_types{s}).(tasks{t}),log_data.(stat_types{s}).(tasks{t}),stat_level_map,pp,save_settings,this_log_filename_prefix);
                
            end
        end
        
        % count dimensions and make upper triangular masks
        pp.n_stat_types=length(stat_types);
        [~,matching_stat_idx]=unique(stat_level_map.stat_gt_levels);
        pp.n_gt_levels=length(matching_stat_idx);
        pp.n_tasks=length(tasks);

        for g=1:pp.n_gt_levels

            s_idx=matching_stat_idx(g);
            gt_level_str=stat_level_map.stat_gt_levels_str{s_idx};

            pp.n_features.(gt_level_str)=length(dcoeff.(stat_types{s_idx})(:,1));
            pp.n_nodes.(gt_level_str)=int16(roots([1 1 -2*pp.n_features.(gt_level_str)])); % assuming n_nets x n_nets, x = n*(n+1)/2 -> n^2 + n - 2x
            pp.n_nodes.(gt_level_str)=pp.n_nodes.(gt_level_str)(end) + pp.remove_matrix_diag.(gt_level_str);
            pp.triu_msk.(gt_level_str)=triu(true(pp.n_nodes.(gt_level_str)),pp.remove_matrix_diag.(gt_level_str));
            pp.ids_triu.(gt_level_str)=find(pp.triu_msk.(gt_level_str));

        end
        
        % reshape residuals
        res.(stat_types{s})=reshape(res.(stat_types{s}),pp.n_features.(stat_level_map.stat_gt_levels_str{s}),pp.n_tasks);
       
    else % whole-brain
        for t=1:length(tasks)
            %log - TODO: consider a single combined log
            this_log_filename_prefix=[log_filename_prefix,'_',tasks{t},'_',stat_types{s}];
            save_settings=summary_tools.write_summary_log(dcoeff.(stat_types{s})(:,t),tpr.(stat_types{s})(:,t),log_data.(stat_types{s}).(all_tasks{t}),stat_level_map,pp,save_settings,this_log_filename_prefix);

        end
    end
end

pp.all_tasks=tasks; % TODO: need for individual tasks, but maybe there's a better way
summary_tools.visualize_tprs(dcoeff,tpr,tpr_fit,dcoeff_windowed,tpr_windowed,tpr_std,res,filename_prefix,pp,stat_level_map,save_settings);



otherwise % catch mis-specified summary_type
    error('Specified summary procedure doesn''t exist')
end
end



