function summarize_fprs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting locally, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/
% This script summarizes false positives from the "fake" task
% This script is mainly used for FWER estimation; the rest of the features are under construction
% Summarization: FWER calculation, comparison with "real" task effect size
% Plot: spatial distribution of false positives
% Usage: summarize_fprs('stat_types','Size_Extent','save_summarized_data',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables
stat_types_default={'Parametric_Bonferroni','Parametric_FDR','Size_Extent','TFCE','Constrained_FWER','Constrained'};
do_fpr=1;
summary_type='summarize_fprs';
tasks={'REST_v_REST2'};
setparams_summary;
compare_effect_size_default=0; % TODO: move

p = inputParser;
% addRequired(p,'date_time_str_results',@ischar);
addOptional(p,'stat_types',stat_types_default);
addOptional(p,'grsize',grsize_default);
addOptional(p,'save_summarized_data',save_settings.defaults.save_summarized_data);
addOptional(p,'make_figs',make_figs_default); % usage of make_figs is under construction
addOptional(p,'save_figs',save_settings.defaults.save_figs);
addOptional(p,'save_logs',save_settings.defaults.save_logs);
addOptional(p,'compare_effect_size',compare_effect_size_default);

parse(p,varargin{:});

save_settings.do.save_summarized_data=p.Results.save_summarized_data;
make_figs=p.Results.make_figs;
save_settings.do.save_figs=p.Results.save_figs;
save_settings.do.save_logs=p.Results.save_logs;
stat_types=p.Results.stat_types;
grsize=p.Results.grsize;
compare_effect_size=p.Results.compare_effect_size;

%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;
task=tasks{1};

if compare_effect_size
    if ~exist('dcoeff','var')
        summary_type_orig=summary_type; % going to temporarily change summary type to get dcoeff data
        summary_type='visualize_gt';
        set_datetimestr_and_files;
        load(combined_summary_filename,'dcoeff');
        summary_type=summary_type_orig;
    end
end

fprintf('* Summarizing false positive benchmarking results.\n');
figure();

for s=1:length(stat_types)
    
    % task- & stat-specific setup
    stat_type=stat_types{s};
    fprintf(['Summarizing FPRs - ',task,'::',stat_type,'\n'])
    set_datetimestr_and_files;
    summary_tools; % summary functions

    % setup paths
    bench_results_basename_prefix=['results__',task,fpr_str,stat_type,'_grsize',num2str(grsize),'_',date_time_str_results.(task)];

    % set results files
    results_filename=[output_dir,bench_results_basename_prefix,'.mat'];
    fpr_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];

%     % set summary prefixes
%     summary_output_dir=[output_dir,task,fpr_str,'_',stat_type,'_summary/'];
%     summary_prefix=[summary_output_dir,'results__',task,fpr_str,stat_type,'_grsize',num2str(grsize),'_',date_time_str_results.(task)];

%     if make_figs==1
    %     % define a few output files to save for testing already created
    %     fpr_vis_filename=[summary_prefix,'_fpr_by_edges.png'];
    %     fpr_vis_filename=[combined_filename_prefix,'_',task,fpr_str,stat_type,'_fp_spatial'];
        fpr_vis_filename=[combined_filename_prefix,'_',task,fpr_str,'_fp_spatial'];
    
        % setup summary output dir
        % TODO: check removal is okay - what else uses summary_output_dir
        %if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
        if ~exist(combined_summary_dir,'dir'); mkdir(combined_summary_dir); end
%     end
    
    %% Check for needed data and existing summaries
    
    save_settings = summary_tools.check_whether_to_save(save_settings,'save_summarized_data','Combined summary data',fpr_summary_filename);
    save_settings = summary_tools.check_whether_to_save(save_settings,'save_figs','Save figs',fpr_vis_filename);
 
    %% Summarize benchmarking results: 'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual'

    % Option 1: Calculate and save positives and simple summaries
    if save_settings.do.save_summarized_data
        load(results_filename);
        
        %% Additional setup

        % count stuff
        n_subs=rep_params.n_subs_subset;
        %n_subs_total=rep_params.n_subs;
        n_edges=size(edge_stats_all,1); %TODO: find right dim
        n_nodes=size(cluster_stats_all,1);

        % re-create upper triangular mask
        triu_msk=triu(true(n_nodes),1);
        ids_triu=find(triu_msk);
        
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

        % before significance masking, make sure positives are in same space as cluster-level stats
        if ~isequal(size(positives),size(cluster_stats_all)) && ~(contains(UI.statistic_type.ui,'FDR') || contains(UI.statistic_type.ui,'Bonferroni'))
            if numel(positives)==numel(cluster_stats_all)
                % reshape positives to matrix to match cluster_stats_all
                positives=reshape(positives,n_nodes,n_nodes,n_repetitions);
                positives_neg=reshape(positives_neg,n_nodes,n_nodes,n_repetitions);
            else
                error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
            end
        end

        % summarize positives, and mask with cluster_stats (all and significant-only)
        positives_total=sum(positives,length(size(positives)));
        positives_total_neg=sum(positives_neg,length(size(positives)));
        
        % double check FWER calculation
        if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA') || strcmp(UI.statistic_type.ui,'Parametric_FDR') || strcmp(UI.statistic_type.ui,'Parametric_Bonferroni')
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
    
        save(fpr_summary_filename,'edge_stats_summary','edge_stats_summary_neg','cluster_stats_summary','cluster_stats_summary_neg','positives','positives_neg','positives_total','positives_total_neg','FWER_manual','FWER_manual_neg','n_repetitions','n_subs_subset','run_time_h','n_perms','-v7.3');
        fprintf(['Saved summary data in ',fpr_summary_filename,'.\n']);
    

    else % Option 2: Use pre-calculated positives to summarize FPR, FWER, etc.

        load(fpr_summary_filename,'positives_total','positives_total_neg','n_repetitions','n_subs_subset','run_time_h','n_perms','FWER_manual','FWER_manual_neg')
        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA') % need for summary in edge_groups
            load(results_filename,'UI');
        end
        
        
        
        
        % set up naming/indexing vars based on levels of inference and ground truth

        stat_level_map.stat_types=all_stat_types;
        stat_level_map.stat_levels_str=all_stat_types;
        for i=1:size(stats_levelstr_map,1)
            stat_level_map.stat_levels_str{strcmp(stat_level_map.stat_levels_str,stats_levelstr_map{i,1})}=stats_levelstr_map{i,2};
        end

        stat_level_map.stat_gt_levels=zeros(1,length(stat_level_map.stat_levels_str));
        stat_level_map.stat_gt_levels_str=cell(1,length(stat_level_map.stat_levels_str));
        for i=1:size(statlevel_gtlevel_map,1)
            idx=contains(stat_level_map.stat_levels_str,statlevel_gtlevel_map{i,1});
            stat_level_map.stat_gt_levels(idx)=statlevel_gtlevel_map{i,3};
            stat_level_map.stat_gt_levels_str(idx)=statlevel_gtlevel_map(i,2);
            % TODO: use this for ground truth categories
        end
        
        % count dimensions and make upper triangular masks
        pp.n_stat_types=length(stat_types);
        [~,matching_stat_idx]=unique(stat_level_map.stat_gt_levels);
        pp.n_gt_levels=length(matching_stat_idx);
        pp.n_tasks=1;

         if contains(stat_type,'Parametric') 
             g=1;
         elseif contains(stat_type,'Constrained')
             g=2;
         elseif contains(stat_type,'Size') || contains(stat_type,'TFCE')
             g=3; % TODO: just a workaround
         end
             
%         for g=1:pp.n_gt_levels

            s_idx=matching_stat_idx(g);
            gt_level_str=stat_level_map.stat_gt_levels_str{s_idx};

%             pp.n_features.(gt_level_str)=length(dcoeff.(stat_types{s_idx})(:,1));
            pp.n_features.(gt_level_str)=length(positives_total);
            pp.n_nodes.(gt_level_str)=int16(roots([1 1 -2*pp.n_features.(gt_level_str)])); % assuming n_nets x n_nets, x = n*(n+1)/2 -> n^2 + n - 2x
            pp.n_nodes.(gt_level_str)=pp.n_nodes.(gt_level_str)(end) + pp.remove_matrix_diag.(gt_level_str);
            pp.triu_msk.(gt_level_str)=triu(true(pp.n_nodes.(gt_level_str)),pp.remove_matrix_diag.(gt_level_str));
            pp.ids_triu.(gt_level_str)=find(pp.triu_msk.(gt_level_str));

%         end

        % Plot
        
        if contains(stat_type,'Parametric')
            p=+pp.triu_msk.edge;
            p(pp.triu_msk.edge)=positives_total;
            p_neg=+pp.triu_msk.edge;
            p_neg(pp.triu_msk.edge)=positives_total_neg;

            % pos
            subplot_tight(4,pp.n_stat_types,s);
            draw_atlas_boundaries(p');
            caxis(pp.clim_fpr.edge);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types+s);
            summarize_matrix_by_atlas(p');
            caxis(pp.clim_fpr.edge_bynet);

            % neg
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*2+s);
            draw_atlas_boundaries(p_neg');
            caxis(pp.clim_fpr.edge);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*3+s);
            summarize_matrix_by_atlas(p_neg');
            caxis(pp.clim_fpr.edge_bynet);
        
        elseif contains(stat_type,'Size') || contains(stat_type,'TFCE')

            %pos
            subplot_tight(4,pp.n_stat_types,s);
            draw_atlas_boundaries(positives_total');
            caxis(pp.clim_fpr.edge);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types+s);
            summarize_matrix_by_atlas(positives_total');
            caxis(pp.clim_fpr.edge_bynet);

            %neg
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*2+s);
            draw_atlas_boundaries(positives_total_neg');
            caxis(pp.clim_fpr.edge);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*3+s);
            summarize_matrix_by_atlas(positives_total_neg');
            caxis(pp.clim_fpr.edge_bynet);
            
        elseif contains(stat_type,'Constrained')
            
            % pos
            p_mat_tmp=+pp.triu_msk.network;
            p_mat_tmp(pp.triu_msk.network)=positives_total; % post
            p_mat_tmp2=p_mat_tmp';
            
            if ~isfield(pp.n_nodes,'edge')
                % have to step in temporarily
                pp.n_nodes.edge=268;
            end
            map=load_atlas_mapping(pp.n_nodes.edge,'subnetwork');
            p_net2edge=summary_to_full_matrix(p_mat_tmp2,map);
            p_net2edge=p_net2edge{1};

            subplot_tight(4,pp.n_stat_types,s);
            summarize_matrix_by_atlas(p_net2edge);
            caxis(pp.clim_fpr.network);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types+s);
            summarize_matrix_by_atlas(p_net2edge);
            caxis(pp.clim_fpr.network);
            
            %neg
            p_mat_tmp=+pp.triu_msk.network;
            p_mat_tmp(pp.triu_msk.network)=positives_total_neg; % neg
            p_mat_tmp2=p_mat_tmp';
            
            map=load_atlas_mapping(pp.n_nodes.edge,'subnetwork');
            p_net2edge=summary_to_full_matrix(p_mat_tmp2,map);
            p_net2edge=p_net2edge{1};

            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*2+s);
            summarize_matrix_by_atlas(p_net2edge);
            caxis(pp.clim_fpr.network);
            
            subplot_tight(4,pp.n_stat_types,pp.n_stat_types*3+s);
            summarize_matrix_by_atlas(p_net2edge);
            caxis(pp.clim_fpr.network);
            
        end
        
        if compare_effect_size
%             if ~exist('dcoeff','var')
%                 summary_type_orig=summary_type; % going to temporarily change summary type to get dcoeff data
%                 summary_type='visualize_gt';
%                 set_datetimestr_and_files;
%                 load(combined_summary_filename,'dcoeff');
%                 summary_type=summary_type_orig;
%             end
            
            if contains(stat_type,'Constrained')
                tmp=tril(true(10));
                d=structure_data(mean(dcoeff.Constrained,2),'mask',tmp);
                d=d(pp.triu_msk.network);
                pos=positives_total;
                pos_neg=positives_total_neg;
            else
                if contains(stat_type,'FDR')
                    pos=positives_total;
                    pos_neg=positives_total_neg;
                else
                    pos=positives_total(pp.triu_msk.edge);
                    pos_neg=positives_total_neg(pp.triu_msk.edge);
                end
                d=mean(dcoeff.Size_Extent,2);
            end
            
            [r__pos_v_esz,p__pos_v_esz]=corr(pos,d);
            [r__pos_v_esz__neg,p__pos_v_esz__neg]=corr(pos_neg,d);
            
            % check association w subnetwork size - TODO: separate out
            load('/Users/steph/Documents/data/mnt/project/benchmarking_results/edge_groups.mat');
            for i=unique(edge_groups(edge_groups~=0))'
                
                msk=edge_groups'==i;
                if ~contains(stat_type,'Constrained')
                    n_edges_in_group(i)=sum(sum(+msk));
                    tmp=triu(true(size(edge_groups)),1);
                    if contains(stat_type,'FDR')
                        pos=structure_data(positives_total,'mask',tmp);
                        pos_neg=structure_data(positives_total_neg,'mask',tmp);
                    else
                        pos=positives_total;
                        pos_neg=positives_total_neg;
                    end
                    pos_in_group__sum(i)=sum(pos(msk));
                    pos_in_group__mean(i)=mean(pos(msk));
                    pos_in_group__sum__neg(i)=sum(pos_neg(msk));
                    pos_in_group__mean__neg(i)=mean(pos_neg(msk));
                else
                    pos=positives_total;
                    pos_neg=positives_total_neg;
                end
                
            end
            
            [r__sum,p__sum]=corr(n_edges_in_group',pos_in_group__sum');
            [r__sum__neg,p__sum__neg]=corr(n_edges_in_group',pos_in_group__sum__neg');
            [r__mean,p__mean]=corr(n_edges_in_group',pos_in_group__mean');
            [r__mean__neg,p__mean__neg]=corr(n_edges_in_group',pos_in_group__mean__neg');
            
        end
       
        

    %     colormap(bipolar([],0.1));

        %% Log percent esz and TP at thresholds
        
    %     this_log_filename_prefix=[log_filename_prefix,'_',task,'_',stat_type];
    %     logfile=[log_filename_prefix,task,fpr_str,stat_type,'_grsize',num2str(grsize),'_',date_time_str_results.(task),'_log.txt'];

        logfile=[log_filename_prefix,'_',task,fpr_str,'_',stat_type,'_fp_log','.txt'];    
        save_settings = summary_tools.check_whether_to_save(save_settings,'save_logs','Save log',logfile);
        
        if save_settings.do.save_logs
            fprintf('Saving log in %s.\n',logfile);
            fid=fopen(logfile,'w');
            fprintf(fid,'Manual FWER: %1.4f (FWER neg: %1.4f)\n',[FWER_manual, FWER_manual_neg]);
            fprintf(fid,'%d total repetitions',n_repetitions);
            fprintf(fid,'\n%s total permutations',n_perms);
            fprintf(fid,'\n%d subjects sampled out of %d total subjects',n_subs_subset);
            fprintf(fid,'\nRun time: %1.2f hours',run_time_h); % toc is in sec
            if compare_effect_size
                fprintf(fid,'\nCorr with mean effect size - pos tail: r=%1.4f , p=%1.4f; neg tail: r=%1.4f , p=%1.4f',r__pos_v_esz,p__pos_v_esz,r__pos_v_esz__neg,p__pos_v_esz__neg);
                fprintf(fid,'\nCorr btw # pos and network size - pos tail: r=%1.4f , p=%1.4f; neg tail: r=%1.4f , p=%1.4f',r__sum,p__sum,r__sum__neg,p__sum__neg);
                fprintf(fid,'\nCorr btw mean pos and network size - pos tail: r=%1.4f , p=%1.4f; neg tail: r=%1.4f , p=%1.4f',r__mean,p__mean,r__mean__neg,p__mean__neg);
            end
            fclose(fid);
        end
        


        colormap(pp.cmap_fpr);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, pp.fig_width_1plt_norm_big*1.5, pp.fig_height_1plt_norm_big*1.7])
        if save_settings.do.save_figs
            print(gcf,fpr_vis_filename,'-dpng','-r300'); 
        end

    end % end of "if save_settings.do.save_summarized_data"
end


