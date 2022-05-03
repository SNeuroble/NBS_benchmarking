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
%       'calculate_tpr':    calculate true positive rates for each task/stat
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
% Requirements:
%   for 'dcoefficients': ground truth test statistics (see calculate_ground_truth.m)
%   for 'positives': benchmarking results
%   for 'tpr': dcoefficients and positives
%   for 'visualize: 'tpr' for all tasks and stat types
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
addpath(genpath(['../',current_path])); % need to add the complete script dir to get the config files
do_fpr=0;
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
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'do_combined',do_combined_default);
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
make_figs=p.Results.make_figs;
pp.do_combined=p.Results.do_combined;
if pp.do_combined; make_figs=0; else; make_figs=p.Results.make_figs; end
if class(stat_types)=='char'; stat_types={stat_types}; end


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

% set up a single mapping from statistic type to level of inference to ground truth level based on separate definitions in setparams
% stat_level_map: stat_types to stat_levels_str to stat_gt_levels to stat_gt_levels_str

% stat_level_map.do_overlay=1; %TODO: set as user-defined param

stat_level_map.stat_types=all_stat_types;
stat_level_map.stat_levels_str=all_stat_types;
for s=1:size(stats_levelstr_map,1)
%     stat_level_map.stat_levels_str{contains(stat_level_map.stat_levels_str,stats_levelstr_map{s,1})}=stats_levelstr_map{s,2};
    stat_level_map.stat_levels_str{strcmp(stat_level_map.stat_levels_str,stats_levelstr_map{s,1})}=stats_levelstr_map{s,2};
end

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
        
        % calculate and save tpr (checking before saving is built in)
        set_datetimestr_and_files; % set datetimestr for specified task, stat_type, data_source
        %if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        save_settings=summary_tools.calculate_positives(results_filename,benchmarking_summary_filename,save_settings);
       
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE AND AGGREGATE TRUE POSITIVES
% i.e., combine effect size and positives above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'calculate_tpr'
fprintf('**** Aggregating true positives ****\n');
if ~use_preaveraged_constrained
if length(tasks)<length(all_tasks)
    error('Not enough tasks specified to combine all.')
end
end

for t=1:length(tasks)
    for s=1:length(stat_types)

        task=tasks{t};
        stat_type=stat_types{s};
        this_stat_level_str=stat_level_map.stat_levels_str{strcmp(stat_level_map.stat_types,stat_type)};
        this_stat_gt_level_str=stat_level_map.stat_gt_levels_str{strcmp(stat_level_map.stat_types,stat_type)};
        fprintf(['Doing: ',task,'::',stat_type,'\n'])

        % set files for specified task, stat_type, data_source and make output
        set_datetimestr_and_files;
        %if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end

        if use_preaveraged_constrained
            this_stat_level_str='edge';
            this_stat_gt_level_str='edge';
        end

        % need original bench results for summary in edge_groups - TODO: save edge_groups into summary matfile
        if strcmp(this_stat_gt_level_str,'network')
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
        [dcoeff.(stat_types{s})(:,t),tpr.(stat_types{s})(:,t),fpr.(stat_types{s})(:,t),fwer_strong.(stat_types{s})(t),fdr.(stat_types{s})(:,t),localizing_power.(stat_types{s})(t),num_fp.(stat_types{s})(:,t),spatial_extent_fp.(stat_types{s})(:,t),log_data.(stat_types{s}).(tasks{t})]...
            =summary_tools.calculate_tpr(benchmarking_summary_filename,ground_truth_dcoeff_filename,this_stat_level_str,this_stat_gt_level_str,pp.remove_matrix_diag,edge_groups);
        % TODO: remove the fwer_strong here bc it's already passed to log_data
        
    end
end

% TODO: remove
% save(sprintf('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/fwer_fdr_lp_indvid_files/lp_fp_gr%d',grsize),'fdr','num_fp','spatial_extent_fp','-append')

% save combined-task data
save_settings = summary_tools.check_whether_to_save(save_settings,'save_summarized_data','Combined summary data',combined_summary_filename); % TODO
if save_settings.do.save_summarized_data
    if ~exist(combined_summary_dir,'dir'); mkdir(combined_summary_dir); end
    save(combined_summary_filename,'dcoeff','tpr','log_data','fpr','fwer_strong','fdr','localizing_power','num_fp','spatial_extent_fp','-v7.3');
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

load(combined_summary_filename,'dcoeff','tpr','log_data','fpr','fwer_strong','fdr','localizing_power','spatial_extent_fp');

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

% TODO: TEMPORARY FOR TESTING ONLY, APPEND TO EXISTING FILES AND REMOVE HERE
% load(sprintf('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/fwer_fdr_lp_indvid_files/lp_fp_gr%d.mat',grsize))

pp.all_tasks=tasks; % TODO: need for individual tasks, but maybe there's a better way
% summary_tools.visualize_tpr(dcoeff,tpr,tpr_fit,dcoeff_windowed,tpr_windowed,tpr_std,res,filename_prefix,pp,stat_level_map,save_settings);
% summary_tools.visualize_tpr(dcoeff,tpr,tpr_fit,dcoeff_windowed,tpr_windowed,tpr_std,res,fpr,fwer_strong,fdr,localizing_power,num_fp,spatial_extent_fp,filename_prefix,pp,stat_level_map,save_settings);
summary_tools.visualize_tpr(dcoeff,tpr,tpr_fit,dcoeff_windowed,tpr_windowed,tpr_std,res,fpr,fwer_strong,fdr,localizing_power,spatial_extent_fp,filename_prefix,pp,stat_level_map,save_settings);


otherwise % catch mis-specified summary_type
    error('Specified summary procedure doesn''t exist. Must be one of: ''dcoeff'', ''positives'', ''calculate_tpr'', ''visualize_gt'', or ''visualize_tpr''.')
end
end



