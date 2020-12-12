function summarize_tprs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script summarizes and visualizes true positive rates
%
% Usage: summarize_tprs('tasks',{'SOCIAL_v_REST'},'stat_types',{'Size_Extent'},'grsize',40,'make_figs',0);
%   Task choices: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION; MOTOR; GAMBLING
%
% Required: Ground truth summary (see calculate_ground_truth), benchmarking results
% 
% STEPS
% Summarization: fits spline to effect size vs. mean TPR
% Plot: d v. TPR spline, d v. TPR residual map
%
% Recommended to first run with local access to data to obtain intermediate 
% summary file--this intermediate file is slow to create but can then be
% used to recreate any summaries/visualizations. Note that this step
% doesn't rely on the ground truth
% (When using remote data, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARSE PARAMETERS

all_tasks={'EMOTION_v_REST','GAMBLING_v_REST','LANGUAGE_v_REST','MOTOR_v_REST','RELATIONAL_v_REST','SOCIAL_v_REST','WM_v_REST'};
%all_tasks={'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};
stat_types_default={'Size_Extent','TFCE','Constrained','Omnibus_Threshold_Both_Dir'};
% omnibus_types_default={''};
combine_all_tasks_default=0;
grsize_default=40;
make_figs_default=1;
make_figs__only_combined_default=0;
save_figs_default=1;
save_log_default=1;
save_combined_default=1;

p = inputParser;
% addRequired(p,'date_time_str_results',@ischar);
addOptional(p,'tasks',all_tasks);
addOptional(p,'stat_types',stat_types_default);
% addOptional(p,'omnibus_types',omnibus_types_default);
addOptional(p,'combine_all_tasks',combine_all_tasks_default);
addOptional(p,'grsize',grsize_default);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'make_figs__only_combined',make_figs__only_combined_default);
addOptional(p,'save_figs',save_figs_default);
addOptional(p,'save_log',save_log_default);
addOptional(p,'save_combined',save_combined_default);
parse(p,varargin{:});

make_figs__only_combined=p.Results.make_figs__only_combined;
if make_figs__only_combined; make_figs=0; else; make_figs=p.Results.make_figs; end
save_figs=p.Results.save_figs;
save_log=p.Results.save_log;
save_combined=p.Results.save_combined;
tasks=p.Results.tasks;
stat_types=p.Results.stat_types;
% omnibus_types=p.Results.omnibus_types;
combine_all_tasks=p.Results.combine_all_tasks;
% date_time_str_results=p.Results.date_time_str_results;
grsize=p.Results.grsize;


%% OTHER SETUP

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarking
addpath(genpath(current_path));
setpaths;
save_settings_for_all.asked.summarize=0;
save_settings_for_all.asked.figs=0;
save_settings_for_all.asked.log=0;
save_settings_for_all.asked.combined=0;

warning('off', 'STRUCTURE_DATA:ASSUME_LOWER_TRI');
warning('off', 'SPLINES:CHCKXYWP:NaNs');

% Check for curve fitting toolbox, required for figs or log
if make_figs || save_log || combine_all_tasks
    v_=ver;
    [installedToolboxes{1:length(v_)}] = deal(v_.Name);
    curve_toolbox_exists = all(ismember('Curve Fitting Toolbox',installedToolboxes));
    if ~curve_toolbox_exists
        error('Curve fitting toolbox required for making summary figures, logs, or combining tasks but not installed. Please re-run with the toolbox installed OR with options ''make_figs'' and ''save_log'' set to 0.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE GROUND TRUTH EFFECTS AND TRUE POSITIVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('* Summarizing true positive benchmarking results.\n');

for t=1:length(tasks)
    for s=1:length(stat_types)
%         for omn=1:length(omnibus_types)
            
            % task-/stat-specific setup
            task=tasks{t};
            stat_type=stat_types{s};
%             if contains(stat_type,'Omnibus')
% %                 omnibus_type=omnibus_types{omn};
%                 omnibus_str=['_',omnibus_type];
%             else
% %                 omnibus_type='';
%                 omnibus_str='';
%             end
            
            fprintf(['Summarizing TPRs - ',task,'::',stat_type,'\n'])
            setparams_summary;
            
            ground_truth_results_basename_prefix=['ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
            bench_results_basename_prefix=['results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
%             bench_results_basename_prefix=['results__',task,'_',stat_type,omnibus_str,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
            %ground_truth_results_basename_prefix=['nbs_ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
            %bench_results_basename_prefix=['nbs_benchmark_results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
            
            % set results filenames
            ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
            results_filename=[output_dir,bench_results_basename_prefix,'.mat'];
            benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];
            
            % set summary prefixes
            summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
            summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_',date_time_str_results.(task)];
%             summary_output_dir=[output_dir,task,'_',stat_type,omnibus_str,'_summary/'];
%             summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,omnibus_str,'_',date_time_str_results.(task)];
            summary_output_dir_gt=[output_dir,task,'_',stat_type_gt,'_summary/'];
            % ground_truth_summary_prefix=[summary_output_dir_gt,'ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
            % ground_truth_summary_prefix=[summary_output_dir_gt,'nbs_ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
            %summary_prefix=[summary_output_dir,'nbs_benchmark_results__',task,'_',stat_type,'_',date_time_str_results.(task)];
            
            % make summary output dir
            if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
            if ~exist(summary_output_dir_gt,'dir'); mkdir(summary_output_dir_gt); end
            
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
                size_cluster_stats_all=size(cluster_stats_all);
                n_repetitions=rep_params.n_repetitions;
                n_dim__cluster_stats_all=length(size_cluster_stats_all); %note that matrices may have different sizes, so we summarize over the last dimension)
                
                % summarize edge and cluster stats
                edge_stats_summary.mean=mean(edge_stats_all,length(size(edge_stats_all)));
                edge_stats_summary.std=std(edge_stats_all,0,length(size(edge_stats_all)));
                edge_stats_summary_neg.mean=mean(edge_stats_all_neg,length(size(edge_stats_all_neg)));
                edge_stats_summary_neg.std=std(edge_stats_all_neg,0,length(size(edge_stats_all_neg)));
                
                cluster_stats_summary.mean=mean(cluster_stats_all,length(size(cluster_stats_all)));
                cluster_stats_summary.std=std(cluster_stats_all,0,length(size(cluster_stats_all)));
                cluster_stats_summary_neg.mean=mean(cluster_stats_all_neg,length(size(cluster_stats_all_neg)));
                cluster_stats_summary_neg.std=std(cluster_stats_all_neg,0,length(size(cluster_stats_all_neg)));
                
                % get positives
                positives=+(pvals_all<str2double(UI.alpha.ui));
                positives_neg=+(pvals_all_neg<str2double(UI.alpha.ui));
                
                % removed this and cluster_stats_sig* calculation below since not used for now (why weight the positives by the effect size? don't we just care about the positives?)
%                 % before significance masking, make sure positives are in same space as cluster-level stats
%                 if ~isequal(size(positives),size(cluster_stats_all))
%                     if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
%                         error('Something went wrong - this shouldn''t happen anymore, only in old summaries created by old script.')
%                     elseif numel(positives)==numel(cluster_stats_all)
%                         
%                         % reshape positives to matrix to match cluster_stats_all
%                         positives=reshape(positives,n_nodes,n_nodes,n_repetitions);
%                         positives_neg=reshape(positives_neg,n_nodes,n_nodes,n_repetitions);
%                         
%                     elseif strcmp(UI.statistic_type.ui,'Omnibus')
%                         warning('This is ONLY a temporary quick fix for the new omnibus.');
%                         cluster_stats_all=squeeze(cluster_stats_all(1,1,:))';
%                     else
%                         error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
%                     end
%                 end
                
                % summarize positives, and mask with cluster_stats (all and significant-only)
                positives_total=sum(positives,length(size(positives)));
                positives_total_neg=sum(positives_neg,length(size(positives)));

                % again, removed bc not used for now (see above)
%                 cluster_stats_sig_all=cluster_stats_all.*positives;
%                 cluster_stats_sig_summary.mean=mean(cluster_stats_sig_all,n_dim__cluster_stats_all);
%                 cluster_stats_sig_summary.std=std(cluster_stats_sig_all,0,n_dim__cluster_stats_all);
%                 
%                 cluster_stats_sig_all_neg=cluster_stats_all_neg.*positives_neg;
%                 cluster_stats_sig_summary_neg.mean=mean(cluster_stats_sig_all_neg,n_dim__cluster_stats_all);
%                 cluster_stats_sig_summary_neg.std=std(cluster_stats_sig_all_neg,0,n_dim__cluster_stats_all);
                
                % double check FWER calculation
                if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA') || strcmp(UI.statistic_type.ui,'Omnibus')
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
                    error(['Unable to load all necessary variables from ',benchmarking_summary_filename,'.\nPlease check these exist and try again.']);
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
            elseif strcmp(stat_type,'Omnibus')
                ids_pos=1;
                ids_neg=1;
            else
                ids_pos=ids_triu(ids_pos_vec);
                ids_neg=ids_triu(ids_neg_vec);
            end
            
            true_positives=zeros(size(dcoeff));
            true_positives(ids_pos_vec)=positives_total(ids_pos);
            true_positives(ids_neg_vec)=positives_total_neg(ids_neg);
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
                

                % If combining tasks, append the current data
                if combine_all_tasks
                    if t==1
                        dcoeff_scaled_all{s}{scaling}=dcoeff_scaled{scaling};
                        tpr_scaled_all{s}{scaling}=tpr_scaled{scaling};
                        if scaling==1
                            log_data_combined.n_repetitions=n_repetitions;
                            log_data_combined.n_perms=n_perms;
                            log_data_combined.n_subs_subset=n_subs_subset;
                            log_data_combined.n_subs_total=n_subs_total;
                            log_data_combined.run_time_h=run_time_h;
                        end
                    else
                        dcoeff_scaled_all{s}{scaling}(:,t)=dcoeff_scaled{scaling};
                        tpr_scaled_all{s}{scaling}(:,t)=tpr_scaled{scaling};
                        if scaling==1
                            if log_data_combined.n_repetitions==n_repetitions; log_data_combined.n_repetitions=nan; end % same val for all or nan
                            if log_data_combined.n_perms==n_perms; log_data_combined.n_perms=nan; end % same val for all or nan
                            if log_data_combined.n_subs_subset==n_subs_subset; log_data_combined.n_subs_subset=nan; end % same val for all or nan
                            log_data_combined.n_subs_total=log_data_combined.n_subs_total+(n_subs_total-log_data_combined.n_subs_total)/t; % avg num subs
                            log_data_combined.run_time_h=log_data_combined.run_time_h+(run_time_h-log_data_combined.run_time_h)/t; % avg run time
                        end
                    end
                end
                
                
            
        end
    end
    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINE ALL TASKS
% note: not too much point running this section if
% not doing figs/log. We don't return any summary 
% info right now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if combine_all_tasks
    
    % re-ask about saving figures and log settings
    save_settings_for_all.asked.figs=0;
    save_log=1;
    save_settings_for_all.asked.log=0;
    
    if length(tasks)<length(all_tasks)
        error('Not enough tasks specified to combine all.')
    end
    
    for s=1:length(stat_types)
        
        % all tasks + stat-specific setup
        
        task='all_tasks';
        stat_type=stat_types{s};
        fprintf(['Summarizing TPRs - ',task,'::',stat_type,'\n'])
        setparams_summary;
        
        % redefine axis limits
        pp.ax_ymax_esz_hist{1}=pp.ax_ymax_esz_hist{1}*6;
        pp.ax_ymax_esz_hist{2}=pp.ax_ymax_esz_hist{2}*6;
        
        % set results filenames
        date_time_str_now=datestr(now,'mmddyyyy');
        bench_results_basename_prefix=['results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_now];
%         bench_results_basename_prefix=['results__',task,'_',stat_type,omnibus_str,'_','grsize',num2str(grsize),'_',date_time_str_now];
        benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];
        
        % set summary prefixes
        summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
        summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_',date_time_str_now];
%         summary_output_dir=[output_dir,task,'_',stat_type,omnibus_str,'_summary/'];
%         summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,omnibus_str,'_',date_time_str_now];
        
        % make summary output dir
        if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        
        % save intermediate data
        [save_settings_for_all,save_combined] = check_whether_to_save(save_settings_for_all,save_combined,'combined','Combined summary data',benchmarking_summary_filename);
        if save_combined
            save(benchmarking_summary_filename,'dcoeff_scaled_all','tpr_scaled_all','log_data_combined','-v7.3');
        end
        
        for scaling=1:length(dcoeff_windowed)            
            
%             dcoeff_scaled{scaling}=dcoeff_scaled_all{s}{scaling}(:);
%             tpr_scaled{scaling}=tpr_scaled_all{s}{scaling}(:);

            % Fit curves for TPR v effect size to all available task data - needed for both visualization and log
            [tpr_fit{scaling},res_scaled{scaling},dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},~]=...
            fit_spline(dcoeff_scaled_all{s}{scaling}(:),tpr_scaled_all{s}{scaling}(:),pp.spline_smoothing{scaling},pp.window_sz{scaling});
            
            % Make figs and log
             if make_figs || make_figs__only_combined
                save_settings_for_all=visualize_tprs(dcoeff_scaled_all{s}{scaling}(:),tpr_scaled_all{s}{scaling}(:),dcoeff_windowed{scaling},tpr_windowed{scaling},tpr_windowed_std{scaling},tpr_fit{scaling},res_scaled{scaling},triu_msk,summary_prefix,pp,scaling,save_figs,save_settings_for_all);
             end
            
            if save_log
                save_settings_for_all=write_summary_to_log(dcoeff_windowed{scaling},tpr_windowed{scaling},log_data_combined,summary_prefix,pp,scaling,save_log,save_settings_for_all); 
            end
            
        end
        
    end
end
end























function save_settings_for_all=visualize_tprs(dcoeff,tpr,dcoeff_windowed,tpr_windowed,tpr_std,tpr_fit,res,triu_msk,summary_prefix,pp,scaling,save_figs,save_settings_for_all)

%% Unpack plot params

% scaling descrip
scaling_str=pp.scaling_str{scaling};

% font
fontsz=pp.fontsz;

% spline parameters
window_sz=pp.window_sz;
spline_smoothing=pp.spline_smoothing;

% axis limits (used for histogram counting too)
ax_ymin=pp.ax_ymin;
ax_ymax_tp=pp.ax_ymax_tp;
ax_xmin_delta=pp.ax_xmin_delta; ax_xmax_delta=pp.ax_xmax_delta; ax_ymax_esz_hist_delta=pp.ax_ymax_esz_hist_delta;
ax_xmin=pp.ax_xmin{scaling}; ax_xmax=pp.ax_xmax{scaling};
ax_ymax_esz_hist=pp.ax_ymax_esz_hist{scaling};

% histograms params (keep an eye out for NAN/empty bins)
bin_width=pp.bin_width{scaling};
tpr_bin_width=pp.tpr_bin_width{scaling};
%     bin_width_at_summary_thresh=pp.bin_width_at_summary_thresh;

% effect size thresholds
thresh_small=pp.thresh_small; thresh_med=pp.thresh_med; thresh_large=pp.thresh_large;

% for visualizing residuals
n_std_residual_outlier=pp.n_std_residual_outlier;

% color limits
%     clim=pp.clim;
clim_res=pp.clim_res{scaling}; % for N=40

% feature size
n_features=pp.n_features{scaling};

% define dcoeff_fit for sanity's sake
dcoeff_fit=dcoeff;

% Check figs already created
esz_v_tpr_file=[summary_prefix,'_tpr_v_esz',scaling_str,'.png'];
if save_figs
    save_figs__results=1;
    [save_settings_for_all,save_figs__results] = check_whether_to_save(save_settings_for_all,save_figs__results,'figs','Results Figures',esz_v_tpr_file);
else
    save_figs__results=0;
end

% 0. Plot TPR histogram

nbins=ceil((ax_ymax_tp-ax_ymin)/tpr_bin_width);
bin_edges=linspace(ax_ymin,ax_ymax_tp,nbins+1);
tpr_histcounts=histcounts(tpr,bin_edges,'Normalization','probability');
bin_centers=bin_edges(1:end-1)+tpr_bin_width/2;
figure; plot(bin_centers,tpr_histcounts,'b-','LineWidth',2)

% h=histogram(tpr,bin_edges,'Normalization','probability');
% h=histfit(tpr,nbins,'kernel');
% figure; plot(h(2).XData,h(2).YData)
% figure; plot(bin_centers,h.Values,'b-','LineWidth',2)
% [y_smhist]=fit_spline(bin_centers,N,tpr_smhist_smoothing,tpr_smhist_window_sz);
% figure; plot(bin_centers,y_smhist,'b-','LineWidth',2)

if save_figs__results
    saveas(gcf,[summary_prefix,'_TPR_hist',scaling_str],'png')
end

% 1. Plot effect size vs. TPR

% first get that hist so can underlay in plot
tmp=figure;
nbins=ceil((ax_xmax-ax_xmin)/bin_width); % edges - good if about 75, 60 for larger axis; nets - should be about 1
bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
h=histogram(dcoeff_fit,bin_edges,'Normalization','probability');
% consider: 
% dcoeff_histcounts=histcounts(dcoeff_fit,bin_edges,'Normalization','probability');
% bin_centers=bin_edges(1:end-1)+bin_width/2;
% figure; plot(bin_centers,dcoeff_histcounts,'b-','LineWidth',2)

figure
hold on
yyaxis left
[~,ind]=sort(dcoeff_fit);
plot(dcoeff_fit(ind),tpr_fit(ind),'b-','LineWidth',2)
if scaling==1
    if exist('shadedErrorBar','file')
        shadedErrorBar(dcoeff_windowed,tpr_windowed,tpr_std,'noLine',1,'lineProps','-b')
    else
        warning('Can''t find shadedError function, so won''t draw shaded error bars.')
    end
else
    %                 if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
    scatter(dcoeff,tpr,1,'b.')
    %         errorbar(dcoeff_plt,tpr_plt,tpr_std_plt,'.')
    %                 else
    %                     shadedErrorBar(dcoeff_windowed,tpr_windowed,tpr_std,'noLine',1,'lineProps','-b')
    %                 end
end
hold off

% add stuff to TPR by esz
axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
set(gca,'fontsize',fontsz)
% add trace of previous hist
hold on
yyaxis right
axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz_hist])
plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/n_features,'--','LineWidth',2)
rectangle('Position',[-thresh_large,ax_ymin,2*thresh_large,ax_ymax_esz_hist],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
hold off

if save_figs__results
    saveas(gcf,esz_v_tpr_file,'png')
end
close(tmp);

% 2. Plot effect size vs. TPR residuals - diagnostics

figure
hold on;
scatter(dcoeff_fit,res,1,'b.')

std_thresh=n_std_residual_outlier*std(res);
idx=abs(res)>std_thresh;
scatter(dcoeff_fit(idx),res(idx),1,'.')
plot(dcoeff_fit,zeros(size(dcoeff_fit)),'k-','LineWidth',2) % plot zero residual line

hold off;

if save_figs__results
    saveas(gcf,[summary_prefix,'_esz_v_TPR__residuals',scaling_str],'png')
end


% 3. Plot effect size vs. TPR residuals - spatial distribution

if scaling==1
    
    if numel(res)==n_features
        % put edge-level residuals back into upper triangle
        res_mat=+triu_msk;
        res_mat(triu_msk)=res;
    else
        res_mat=+triu_msk;
        res_mat(triu_msk)=mean(reshape(res,n_features,numel(res)/n_features),2);
    end
    
    draw_atlas_boundaries(res_mat');
    caxis(clim_res); colormap(bipolar([],0.1));
    if save_figs__results; saveas(gcf,[summary_prefix,'_residuals',scaling_str],'png'); end
    
    summarize_matrix_by_atlas(res_mat');
    caxis(clim_res); colormap(bipolar([],0.1));
    if save_figs__results; saveas(gcf,[summary_prefix,'_residuals',scaling_str,'_in_networks'],'png'); end
    
end

end



function save_settings_for_all=write_summary_to_log(dcoeff_windowed,tpr_windowed,log_data,summary_prefix,pp,scaling,save_log,save_settings_for_all)


%% Mean TPR within effect size thresholds

% unpack plot params
scaling_str=pp.scaling_str{scaling};
thresholds=[pp.thresh_small, pp.thresh_med, pp.thresh_large];
bin_width_at_summary_thresh=pp.bin_width_at_summary_thresh{scaling};

for i=1:length(thresholds)
    
    % Get IDs of edges below/between d-thresh (e.g., edges < thresh_high; edges < thresh_high & > thresh_low)
    ids_lt_thr=abs(dcoeff_windowed) <= thresholds(i);
    if i~=1
        ids_btw_thr_and_thr_below=abs(dcoeff_windowed) <= thresholds(i) & abs(dcoeff_windowed) >= thresholds(i-1);
    end
    
    % -> Get TPR of edges below/between d-threshold
    tpr_lt_thr(i) = nanmean(tpr_windowed(ids_lt_thr));
    if i~=1
        tpr_btw_thr_and_thr_below(i-1) = nanmean(tpr_windowed(ids_btw_thr_and_thr_below));
    end
    
    % Get IDs of edges at (around) dcoeff (divide by 2 to get both halves of the bin)
    ids_at_pos_thr = abs(dcoeff_windowed-thresholds(i)) <= bin_width_at_summary_thresh / 2;
    ids_at_neg_thr = abs(dcoeff_windowed+thresholds(i)) <= bin_width_at_summary_thresh / 2;
    
    % -> Get TPR of edges at (around) +thresh, -thresh, and mean
    tpr_at_thr(3*(i-1) + 1) = nanmean(tpr_windowed(ids_at_pos_thr));
    tpr_at_thr(3*(i-1) + 2) = nanmean(tpr_windowed(ids_at_neg_thr));
    tpr_at_thr(3*(i-1) + 3) = mean(tpr_at_thr((3*i-2):(3*i-1)));
    
end

%% Log percent esz and TP at thresholds

logfile=[summary_prefix,'_log',scaling_str,'.txt'];
[save_settings_for_all,save_log] = check_whether_to_save(save_settings_for_all,save_log,'log','Log file',logfile);

if save_log
    fprintf('Saving log in %s.\n',logfile);
    
    fid=fopen(logfile,'w');
    fprintf(fid,'Mean TPR between d=+/-%1.1f: %1.3f\n',[thresholds; tpr_lt_thr]);
    fprintf(fid,'Mean TPR between d=%1.1f and %1.1f: %1.3f\n',[thresholds(2:end); thresholds(1:end-1); tpr_btw_thr_and_thr_below]);
    fprintf(fid,'Mean TPR at d=+/-%1.1f: %f (+), %f (-), %f (mean)\n',[thresholds; reshape(tpr_at_thr,3,length(thresholds))]);
    fprintf(fid,'%d total repetitions',log_data.n_repetitions);
    fprintf(fid,'\n%s total permutations',log_data.n_perms);
    fprintf(fid,'\n%d subjects sampled out of %d total subjects',log_data.n_subs_subset,round(log_data.n_subs_total));
    fprintf(fid,'\nRun time: %1.2f hours',log_data.run_time_h); % toc is in sec
    fclose(fid);
end
end



function [save_settings_for_all,do_thing] = check_whether_to_save(save_settings_for_all,do_thing,step_name,descriptive_string,filename)

if do_thing
    if exist(filename, 'file') == 2
        if ~save_settings_for_all.asked.(step_name) || ~save_settings_for_all.use_same.(step_name)
            
            user_response=input(sprintf([descriptive_string,' already exists. Overwrite? [yes/no]\n> ']),'s');
            if strcmp(user_response,'yes')
                fprintf(['Replacing previous ',descriptive_string,'.\n']);
%                 warning('DONT WANT TO SAVE FOR NOW - NOT REPLACING'); do_thing=0; % TODO: remove
            else
                fprintf(['Keeping existing ',descriptive_string,'.\n']);
                do_thing=0;
            end
            
            if ~save_settings_for_all.asked.(step_name)
                user_response=input(sprintf('Repeat for all? [yes/no]\n> '),'s');
                if strcmp(user_response,'yes')
                    fprintf('Using this setting for all.\n');
                    save_settings_for_all.use_same.(step_name)=1;
                    save_settings_for_all.(step_name)=do_thing;
                else
                    fprintf('Okay, will ask each time.\n');
                    save_settings_for_all.use_same.(step_name)=0;
                end
                save_settings_for_all.asked.(step_name)=1;
            end
            
        else
            do_thing=save_settings_for_all.(step_name);
        end
        
    end
end

end
