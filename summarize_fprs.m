function summarize_fprs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting locally, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/
% This script summarizes and visualizes true positive rates
% Summarization: fits spline to effect size vs. mean TPR
% Plot: d v. TPR spline, d v. TPR residual map
% Usage: summarize_tprs('LANGUAGE','Size_Extent','02102020_1759',40);
%   Task choices: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION; MOTOR; GAMBLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables
stat_types_default={'Size_Extent','TFCE','Constrained'};
grsize_default=40;
make_figs_default=1;
save_figs_default=1;
save_log_default=1;

p = inputParser;
% addRequired(p,'date_time_str_results',@ischar);
addOptional(p,'stat_types',stat_types_default);
addOptional(p,'grsize',grsize_default);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'save_figs',save_figs_default);
addOptional(p,'save_log',save_log_default);
parse(p,varargin{:});

make_figs=p.Results.make_figs;
save_figs=p.Results.save_figs;
save_log=p.Results.save_log;
stat_types=p.Results.stat_types;
% date_time_str_results=p.Results.date_time_str_results;
grsize=p.Results.grsize;

%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;
save_settings_for_all.asked.summarize=0;
save_settings_for_all.asked.figs=0;
save_settings_for_all.asked.log=0;

fprintf('* Summarizing false positive benchmarking results.\n');

for s=1:length(stat_types)
    
    % stat-specific setup
    task='REST';
    stat_type=stat_types{s};
    fprintf(['Summarizing FPRS - ',task,'::',stat_type,'\n'])
    setparams_summary;

    fpr_results_basename_prefix=['nbs_benchmark_results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];

    % set results files
    results_filename=[output_dir,bench_results_basename_prefix,'.mat'];
    benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];

    % set summary prefixes
    summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
    summary_output_dir_gt=[output_dir,task,'_',stat_type_gt,'_summary/'];
    summary_prefix=[summary_output_dir,'nbs_benchmark_results__',task,'_',stat_type,'_',date_time_str_results.(task)];

    % define a few output files to save for testing already created
    fpr_by_edges_file=[summary_prefix,'_fpr_by_edges.png'];
    logfile=[summary_prefix,'_log.txt'];

    % setup summary output dir
    if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end
    
    
    %% Check for needed data and existing summaries

    summarize_benchmarking=1;
    if exist(benchmarking_summary_filename, 'file') == 2
        if ~save_settings_for_all.asked.summarize || ~save_settings_for_all.use_same.summarize
         
            user_response=input(sprintf('Summary data already exists. Overwrite? [yes/no]\n> '),'s');
            if strcmp(user_response,'yes')
                fprintf('Replacing previous summary.\n');
            else
                fprintf('Keeping existing summary.\n');
                summarize_benchmarking=0;
            end
            
            if ~save_settings_for_all.asked.summarize
                user_response=input(sprintf('Repeat for all? [yes/no]\n> '),'s');
                if strcmp(user_response,'yes')
                    fprintf('Using this setting for all.\n');
                    save_settings_for_all.use_same.summarize=1;
                    save_settings_for_all.summarize=summarize_benchmarking;
                else
                    fprintf('Okay, will ask each time.\n');
                    save_settings_for_all.use_same.summarize=0;
                end
                save_settings_for_all.asked.summarize=1;
            end
            
        else
            summarize_benchmarking=save_settings_for_all.summarize;
        end
        
    end

    if make_figs && save_figs
        save_figs__results=1;
        if exist(fpr_by_edges_file,'file')
            if ~save_settings_for_all.asked.figs || ~save_settings_for_all.use_same.figs
                resp=input(sprintf('Results figures already exist in %s. \nOverwrite? (Otherwise will plot without saving.) [y/n]\n> ',fpr_by_edges_file),'s');
                if strcmp(resp,'y')
                    fprintf('Replacing results figures.\n');
                else
                    save_figs__results=0;
                    fprintf('Okay, won''t overwrite.\n');
                end

                if ~save_settings_for_all.asked.figs
                    user_response=input(sprintf('Repeat for all? [yes/no]\n> '),'s');
                    if strcmp(user_response,'yes')
                        fprintf('Using this setting for all.\n');
                        save_settings_for_all.use_same.figs=1;
                        save_settings_for_all.figs=save_figs__results;
                    else
                        fprintf('Okay, will ask each time.\n');
                        save_settings_for_all.use_same.figs=0;
                    end
                    save_settings_for_all.asked.figs=1;
                end

            else
                save_figs__results=save_settings_for_all.figs;
            end
        end
    else
        save_figs__results=0;
    end

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


    %% Summarize benchmarking results: 'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual'
    
    if summarize_benchmarking

        load(results_filename);
        
        %% Additional setup

        % count stuff
        n_subs=rep_params.n_subs_subset;
        n_subs_total=rep_params.n_subs;
        n_edges=length(edge_stats);
        n_nodes=size(cluster_stats,1);

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
        if ~isequal(size(positives),size(cluster_stats_all))
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
        if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
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
        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA') % need for summary in edge_groups
            load(results_filename,'UI');
        end
    end

    %% Calculate TPR

    ids_pos_vec=dcoeff>0;
    ids_neg_vec=dcoeff<0;

    if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
        edge_groups_triu=UI.edge_groups.ui';
        edge_groups_vec=edge_groups_triu(ids_triu);
        ids_pos=edge_groups_vec(ids_pos_vec);
        ids_neg=edge_groups_vec(ids_neg_vec);
    else
        ids_pos=ids_triu(ids_pos_vec);
        ids_neg=ids_triu(ids_neg_vec);
    end

    true_positives=zeros(size(dcoeff));
    true_positives(ids_pos_vec)=positives_total(ids_pos);
    true_positives(ids_neg_vec)=positives_total_neg(ids_neg);
    tpr=true_positives*100/n_repetitions;
    
    
    %% Mean TPR within effect size thresholds

    thresholds=[thresh_small, thresh_med, thresh_large];

    for i=1:length(thresholds)

        % Get IDs of edges below/between d-thresh (e.g., edges < thresh_high; edges < thresh_high & > thresh_low)
        ids_lt_thr=abs(dcoeff) <= thresholds(i);
        if i~=1
            ids_btw_thr_and_thr_below=abs(dcoeff) <= thresholds(i) & abs(dcoeff) >= thresholds(i-1);
        end

        % -> Get TPR of edges below/between d-threshold
        tpr_lt_thr(i) = sum(tpr(ids_lt_thr)) / sum(+ids_lt_thr);
        if i~=1
            tpr_btw_thr_and_thr_below(i-1) = sum(tpr(ids_btw_thr_and_thr_below)) / sum(+ids_btw_thr_and_thr_below) ;
        end

        % Get IDs of edges at (around) dcoeff (divide by 2 to get both halves of the bin)
        ids_at_pos_thr = abs(dcoeff-thresholds(i)) <= bin_width_at_summary_thresh / 2;
        ids_at_neg_thr = abs(dcoeff+thresholds(i)) <= bin_width_at_summary_thresh / 2;

        % -> Get TPR of edges at (around) +thresh, -thresh, and mean
        tpr_at_thr(3*(i-1) + 1) = sum(tpr(ids_at_pos_thr)) / sum(+ids_at_pos_thr);
        tpr_at_thr(3*(i-1) + 2) = sum(tpr(ids_at_neg_thr)) / sum(+ids_at_neg_thr);
        tpr_at_thr(3*(i-1) + 3) = mean(tpr_at_thr((3*i-2):(3*i-1)));

    end
    
    
    %% Fit TPR v effect size
    % curve fitting toolbox required - check - thanks https://www.mathworks.com/matlabcentral/fileexchange/51794-istoolboxavailable

    v_=ver;
    [installedToolboxes{1:length(v_)}] = deal(v_.Name);
    curve_toolbox_exists = all(ismember('Curve Fitting Toolbox',installedToolboxes));
    if curve_toolbox_exists
        tpr_fit=zeros(n_edges,1);
        res=zeros(n_edges,1);
        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
            dcoeff_mat=zeros(n_nodes);
            dcoeff_mat(triu_msk)=dcoeff;

            dcoeff_summat=summarize_matrix_by_atlas(dcoeff_mat')';
            close
            triu_msk_summat=logical(triu(ones(size(dcoeff_summat))));
            dcoeff_summat=dcoeff_summat(triu_msk_summat);

            ids_pos_summat=dcoeff_summat>0;
            ids_neg_summat=dcoeff_summat<0;

            %TBD
            true_positives_summat(ids_pos_summat)=positives_total(ids_pos_summat);
            true_positives_summat(ids_neg_summat)=positives_total_neg(ids_neg_summat);

            tpr_summat=true_positives_summat*100/n_repetitions;

            spline_smoothing_set=0.99999;
            [tpr_fit_summat,res_summat,dcoeff_windowed_summat,tpr_windowed_summat,tpr_windowed_std_summat,~]=fit_spline(dcoeff_summat,tpr_summat,spline_smoothing_set,[summary_prefix,'_esz_v_TPR_summat_pos']);
            [tpr_fit_by_edge,res_by_edge,dcoeff_windowed_by_edge,tpr_windowed_by_edge,tpr_windowed_std_by_edge,~]=fit_spline(dcoeff,tpr,spline_smoothing,[summary_prefix,'_esz_v_TPR_pos']);
        else
            [tpr_fit,res,dcoeff_windowed,tpr_windowed,tpr_windowed_std,~]=fit_spline(dcoeff,tpr,spline_smoothing,[summary_prefix,'_esz_v_TPR_pos']);
        %         [tpr_fit(ids_neg),res(ids_neg),~]=fit_spline(dcoeff(ids_neg),true_positives(ids_neg,spline_smoothing,strcat(out_prefix,'_esz_v_TPR_neg'));
        end
    else
        warning('Curve fitting toolbox required for fitting spline but not installed - you won''t be able to plot residuals.');
    end



    %% VISUALIZATION

    if make_figs && curve_toolbox_exists

        % Plot results from fitting spline

        % first get that hist so can underlay in plot
        tmp=figure;
        bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
        h=histogram(dcoeff,bin_edges,'Normalization','probability');

        % 1. Plot effect size vs. TPR

        if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
    %         dcoeff_plt=dcoeff_summat;
    %         tpr_plt=tpr_summat;
            dcoeff_plt=dcoeff_windowed_by_edge;
            tpr_plt=tpr_windowed_by_edge;
            tpr_std_plt=tpr_windowed_std_by_edge;

            dcoeff_fit_plt=dcoeff;
            tpr_fit_plt=tpr_fit_by_edge;
            
            res_plt=res_by_edge;
        else
    %         dcoeff_plt=dcoeff;
    %         tpr_plt=tpr;
            dcoeff_plt=dcoeff_windowed;
            tpr_plt=tpr_windowed;
            tpr_std_plt=tpr_windowed_std;

            dcoeff_fit_plt=dcoeff;
            tpr_fit_plt=tpr_fit;
            
            res_plt=res;
        end

        figure
        hold on
        yyaxis left
        [~,ind]=sort(dcoeff_fit_plt);
        plot(dcoeff_fit_plt(ind),tpr_fit_plt(ind),'k-','LineWidth',2)
        if exist('shadedErrorBar','file')
            shadedErrorBar(dcoeff_plt,tpr_plt,tpr_std_plt,'noLine',1)
    %         scatter(dcoeff_plt,tpr_plt,1,'b.')
    %         errorbar(dcoeff_plt,tpr_plt,tpr_std_plt,'.')
        else
            warning('Must add function shadedError to the path to draw shaded error bars.')
        end
        hold off

        % add stuff to TPR by esz
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
        set(gca,'fontsize',fontsz)
        % add trace of previous hist
        hold on
        yyaxis right
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
        plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/n_edges,'--','LineWidth',2)
        rectangle('Position',[-thresh_large,ax_ymin,2*thresh_large,ax_ymax_tp],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
        hold off

        if save_figs__results
            saveas(gcf,fpr_by_edges_file,'png')
        end
        close(tmp);

        % 2. Plot effect size vs. TPR residuals - diagnostics

        figure
        hold on;
        scatter(dcoeff_fit_plt,res_plt,1,'b.')

        std_thresh=n_std_residual_outlier*std(res_plt);
        idx=abs(res_plt)>std_thresh;
        scatter(dcoeff_fit_plt(idx),res_plt(idx),1,'.')
        plot(dcoeff_fit_plt,zeros(size(dcoeff_fit_plt)),'k-','LineWidth',2) % plot zero residual line

        hold off;

        if save_figs__results
            % save plot
            saveas(gcf,[summary_prefix,'_esz_v_TPR__residuals'],'png')
        end


        % 3. Plot effect size vs. TPR residuals - spatial distribution
        % edge-level results

        % put stuff back into upper triangle
        res_mat=zeros(n_nodes);
        res_mat(triu_msk)=res_plt;

        draw_atlas_boundaries(res_mat');
        colormap(bipolar([],0.1));
        caxis(clim_res_detailed);

        if save_figs__results
            saveas(gcf,[summary_prefix,'_residuals_by_edges'],'png')
        end

        % network-level results
        summarize_matrix_by_atlas(res_mat');
        colormap(bipolar([],0.1));
        caxis(clim_res);

        if save_figs__results
            saveas(gcf,[summary_prefix,'_residuals_by_networks'],'png')
        end



    end

    %% Log percent esz and TP at thresholds

    if save_log
        fprintf('Saving log in %s.\n',logfile);

        fid=fopen(logfile,'w');
        fprintf(fid,'Mean TPR between d=+/-%1.1f: %1.3f\n',[thresholds; tpr_lt_thr]);
        fprintf(fid,'Mean TPR between d=%1.1f and %1.1f: %1.3f\n',[thresholds(2:end); thresholds(1:end-1); tpr_btw_thr_and_thr_below]);
        fprintf(fid,'Mean TPR at d=+/-%1.1f: %f (+), %f (-), %f (mean)\n',[thresholds; reshape(tpr_at_thr,3,length(thresholds))]);
        fprintf(fid,'%d total repetitions',n_repetitions);
        fprintf(fid,'\n%s total permutations',n_perms);
        fprintf(fid,'\n%d subjects sampled out of %d total subjects',n_subs_subset,n_subs_total);
        fprintf(fid,'\nRun time: %1.2f hours',run_time_h); % toc is in sec
        fclose(fid);
    end
    
end



