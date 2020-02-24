function summarize_TPRs(task_type,stat_type,date_time_str_ground_truth,date_time_str_results,grsize,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/
% e.g., summarize_TPRs('LANGUAGE','Size_Extent','02102020_1759');
% Task can be: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION
% This script summarizes and visualizes ground truth effect sizes % TBD: and TPR data
% TBD: Summarization: bins edges by effect size, gets mean TPR within bins,
% and fits spline to effect size vs. mean TPR
% TBD: Plot: binned effect size, binned TPR,  d v. TPR spline, d v. TPR residual map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables

make_figs_default=1;
save_figs_default=1;
save_log_default=1;
% calc_TPR_default=1;

p = inputParser;
addRequired(p,'task_type',@ischar);
addRequired(p,'stat_type',@ischar);
addRequired(p,'date_time_str_ground_truth',@ischar);
addRequired(p,'date_time_str_results',@ischar);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'save_figs',save_figs_default);
addOptional(p,'save_log',save_log_default);
% addOptional(p,'calc_TPR',calc_TPR_default);
parse(p,task_type,stat_type,date_time_str_ground_truth,date_time_str_results,varargin{:});

task_type=p.Results.task_type;
stat_type=p.Results.stat_type;
date_time_str_ground_truth=p.Results.date_time_str_ground_truth;
date_time_str_results=p.Results.date_time_str_results;
make_figs=p.Results.make_figs;
save_figs=p.Results.save_figs;
save_log=p.Results.save_log;
% calc_TPR=p.Results.calc_TPR;


%% Plot params

ax_xmin=-2.5; ax_xmax=2.5; ax_ymin=0; ax_ymax_esz=0.15; ax_ymax_tp=100;
fontsz=25;
spline_smoothing=0.995;

% histograms params (keep an eye out for NAN/empty bins)
nbins=75;
half_bin_width=0.05; % ad hoc bin size

% effect size thresholds
thresh_small=0.2; thresh_med=0.5; thresh_high=0.8;

% cmap threshold
clim=[-thresh_med, thresh_med];
%clim_res=[-0.001,0.001]; % for N=20
if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
    clim_res=[-10,10]; % for N=40
    clim_res_detailed=[-60,60]; % for N=40
else
    clim_res=[-0.5,0.5]; % for N=40
    clim_res_detailed=[-3,3]; % for N=40
end



%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;

ground_truth_results_basename_prefix=['nbs_ground_truth__',task_type,'_Size_Extent_',date_time_str_ground_truth];
bench_results_basename_prefix=['nbs_benchmark_results__',task_type,'_',stat_type,'_','_grsize',grsize,'_',date_time_str_results];

% set results files
ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
results_filename=[output_dir,bench_results_basename_prefix,'.mat'];
benchmarking_summary_filename=[output_dir,bench_results_basename_prefix,'_summary.mat'];

% set summary prefixes
summary_output_dir=[output_dir,task_type,'_',stat_type,'_summary/'];
ground_truth_summary_prefix=[summary_output_dir,'nbs_ground_truth__',task_type,'_',stat_type,'_',date_time_str_ground_truth];
summary_prefix=[summary_output_dir,'nbs_benchmark_results__',task_type,'_',stat_type,'_',date_time_str_results];

% define a few output summary files to test existence
esz_hist_file=[ground_truth_summary_prefix,'_esz_hist.png'];
esz_v_tpr_file=[summary_prefix,'_tpr_v_esz.png'];
logfile=[summary_prefix,'_log.txt'];

% setup summary output dir
mkdir(summary_output_dir);

%% Check for existing summaries

summarize_benchmarking=1;
if exist(benchmarking_summary_filename, 'file') == 2
    user_response=input(sprintf('Summary data already exists. Overwrite? [yes/no]\n> '),'s');
    if strcmp(user_response,'yes')
        fprintf('Replacing previous summary.\n');
    else
        fprintf('Keeping existing summary.\n');
        summarize_benchmarking=0;
    end
end

if make_figs && save_figs
    save_figs__gt=1;
    save_figs__results=1;
    if exist(esz_hist_file,'file')
        resp=input(sprintf('Ground truth figures already exist in %s. \nOverwrite? (Otherwise will plot without saving.) [y/n]\n> ',esz_hist_file),'s');
        if strcmp(resp,'y')
            fprintf('Replacing ground truth figures.\n');
        else
            save_figs__gt=0;
            fprintf('Okay, won''t overwrite.\n');
        end
    end
    if exist(esz_v_tpr_file,'file')
        resp=input(sprintf('Results figures already exist in %s. \nOverwrite? (Otherwise will plot without saving.) [y/n]\n> ',esz_v_tpr_file),'s');
        if strcmp(resp,'y')
            fprintf('Replacing results figures.\n');
        else
            save_figs__results=0;
            fprintf('Okay, won''t overwrite.\n');
        end
    end
end

if save_log
    if exist(logfile,'file')
        resp=input('Log file already exists. \nOverwrite? [y/n]\n> ','s');
        if strcmp(resp,'y')
            fprintf('Replacing log.\n');
        else
            save_log=0;
            fprintf('Okay, won''t overwrite.\n');
        end
    end
end


%% Load ground truth and estimate d-coefficients

load(ground_truth_filename);

% t-stat -> d-coefficient - transpose because need for fitting spline
n_subs=size(UI.design.ui,2)-1;
dcoeff=(edge_stats/sqrt(n_subs))';

% other
n_nodes=size(cluster_stats,1); % TODO: this is also set above
n_edges=length(dcoeff);

% re-create upper triangular mask
triu_msk=triu(true(n_nodes),1);
ids_triu=find(triu_msk);

% to visualize residual outliers, highlight points less or greater than 2 std
n_std_residual_outlier=2;


%% Load and summarize benchmarking results: 'edge_stats_summary','cluster_stats_summary','positives','positives_total','FWER_manual'

if summarize_benchmarking
    %note that these matrices may have different sizes, so we summarize over the last dimension)
    
    
    load(results_filename);
    size_cluster_stats_all=size(cluster_stats_all);
    n_repetitions=size_cluster_stats_all(end);
    n_dim__cluster_stats_all=length(size_cluster_stats_all);
    % TODO: this includes positive TPs only, not negative - include negative
    
    
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
        
        if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
            
%             % the old cNBS saves cluster_stats_all as matrix - convert to
%             % vector to match positives
%             edge_groups_msk=UI.edge_groups.ui;
%             % Force mask to be upper triangular only if not already
%             tmp=tril(edge_groups_msk,-1); if any(tmp(:)); edge_groups_msk=triu(edge_groups_msk'); end
%             groups=unique(edge_groups_msk); groups=groups(2:end);
%             n_groups=length(groups);
%             cluster_stats_summary_mat=cluster_stats_summary;
%             cluster_stats_summary.mean=zeros(n_groups,1);
%             cluster_stats_summary.std=zeros(n_groups,1);
%             n_repetitions=size(positives,2);
%             cluster_stats_all_mat=cluster_stats_all;
%             cluster_stats_all=zeros(n_groups,n_repetitions);
%             for i=1:n_groups
%                 [mat_idx_i,mat_idx_j]=find(edge_groups_msk==groups(i),1);
%                 cluster_stats_summary.mean(i)=cluster_stats_summary_mat.mean(mat_idx_i,mat_idx_j);
%                 cluster_stats_summary.std(i)=cluster_stats_summary_mat.std(mat_idx_i,mat_idx_j);
%                 for j=1:n_repetitions
%                     cluster_stats_all(i,j)=cluster_stats_all_mat(mat_idx_i,mat_idx_j,j);
%                 end
%             end
            error('Something went wrong - this shouldn''t happen anymore, only in old data.')
            
        elseif numel(positives)==numel(cluster_stats_all)
            
            % reshape positives to matrix to match cluster_stats_all
            n_nodes=size(cluster_stats_all,1);
            n_repetitions=size(cluster_stats_all,3);
            positives=reshape(positives,n_nodes,n_nodes,n_repetitions);
             positives_neg=reshape(positives_neg,n_nodes,n_nodes,n_repetitions);
            
        else
            error('Cluster stats and p-value dimensions don''t match. We can only fix this in two ways and they must have failed.')
        end
        
    end
    
    % summarize positives, and mask with cluster_stats (all and
    % significant-only)
    positives_total=sum(positives,length(size(positives)));
     positives_total_neg=sum(positives_neg,length(size(positives)));
    cluster_stats_sig_all=cluster_stats_all.*positives; % why weight the positives by the effect size? don't we just care about the positives?
    cluster_stats_sig_summary.mean=mean(cluster_stats_sig_all,n_dim__cluster_stats_all);
    cluster_stats_sig_summary.std=std(cluster_stats_sig_all,0,n_dim__cluster_stats_all);
     cluster_stats_sig_all_neg=cluster_stats_all_neg.*positives_neg;
     cluster_stats_sig_summary_neg.mean=mean(cluster_stats_sig_all_neg,n_dim__cluster_stats_all);
     cluster_stats_sig_summary_neg.std=std(cluster_stats_sig_all_neg,0,n_dim__cluster_stats_all);
    
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
    
    n_repetitions=size(positives,3);

    save(benchmarking_summary_filename,'edge_stats_summary','edge_stats_summary_neg','cluster_stats_summary','cluster_stats_summary_neg','positives','positives_neg','positives_total','positives_total_neg','FWER_manual','FWER_manual_neg','n_repetitions','-v7.3');
else
    load(benchmarking_summary_filename)
    size_positives=size(positives);
    n_repetitions=size_positives(end);
end


%% Calc percent edges within thresholds

% get voxels within thresholds (vox < thresh_high; vox < thresh_high & > thresh_low)
ids_lt_thr_med=abs(dcoeff) <= thresh_med;
ids_btw_thr_med_and_small=abs(dcoeff) <= thresh_med & abs(dcoeff) >= thresh_small;

% calc percent edges within dcoeff thresholds
perc_edges_lt_thr_med=sum(+ids_lt_thr_med) * 100 / n_edges;
perc_edges_btw_thr_med_and_small=sum(+ids_btw_thr_med_and_small) * 100 / n_edges;

%% Calculate TPR from positives
    
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
true_positives(ids_neg_vec)=positives_total_neg(ids_neg); % TODO: all the neg stuff
tpr=true_positives*100/n_repetitions;

% mean TPR within thresholds
% TODO: double check this - not sure this is TPR - should be TPR, not just TP
tpr_lt_thr_med = sum(tpr(ids_lt_thr_med)) / sum(+ids_lt_thr_med);
tpr_btw_thr_med_and_small = sum(tpr(ids_btw_thr_med_and_small)) / sum(+ids_btw_thr_med_and_small) ;

% mean TPR "at" (around) thresholds
t=[thresh_med, thresh_small];
for i=1:length(t)
    ids_at_pos_thr = abs(dcoeff-t(i)) <= half_bin_width;
    ids_at_neg_thr = abs(dcoeff+t(i)) <= half_bin_width;
    tpr_at_thr(3*(i-1) + 1) = sum(tpr(ids_at_pos_thr)) / sum(+ids_at_pos_thr);
    tpr_at_thr(3*(i-1) + 2) = sum(tpr(ids_at_neg_thr)) / sum(+ids_at_neg_thr);
    tpr_at_thr(3*(i-1) + 3) = mean(tpr_at_thr((3*i-2):(3*i-1)));
end

% fit TPR v effect size
% curve fitting toolbox required - check - thanks https://www.mathworks.com/matlabcentral/fileexchange/51794-istoolboxavailable
v_=ver;
[installedToolboxes{1:length(v_)}] = deal(v_.Name);
curve_toolbox_exists = all(ismember('Curve Fitting Toolbox',installedToolboxes));
if curve_toolbox_exists
    tpr_fit=zeros(n_edges,1);
    res=zeros(n_edges,1);
    [tpr_fit,res,~]=fit_spline(dcoeff,tpr,spline_smoothing,[summary_prefix,'_esz_v_TPR_pos']);
    %         [tpr_fit(ids_neg),res(ids_neg),~]=fit_spline(dcoeff(ids_neg),true_positives(ids_neg,spline_smoothing,strcat(out_prefix,'_esz_v_TPR_neg'));
else
    warning('Curve fitting toolbox required for fitting spline but not installed - you won''t be able to get residuals,');
end

%% Visualize

if make_figs 
    
    % 1. Plot effect size histograms

    bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
    h=histogram(dcoeff,bin_edges,'Normalization','probability');
    hold on;
    plot(h.BinEdges(1:end-1) + h.BinWidth/2, h.BinCounts/length(dcoeff))
    hold off;

    % add stuff to hist
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
    set(gca,'fontsize',fontsz)
    % highlight
    hold on
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off

    if save_figs__gt
        saveas(gcf,esz_hist_file,'png')
    end
    

    % 2. Plot effect size spatial distributions

    % put stuff back into upper triangle
    dcoeff_mat=zeros(n_nodes);
    dcoeff_mat(triu_msk)=dcoeff;
    
    % edge-level results
    draw_atlas_boundaries(dcoeff_mat');
    colormap(bipolar([],0.1));
    caxis(clim);
    if save_figs__gt
        saveas(gcf,[ground_truth_summary_prefix,'_esz_by_edges'],'png')
    end
    
    % network-level results
    summarize_matrix_by_atlas(dcoeff_mat');
    colormap(bipolar([],0.1));
    caxis(clim);

    if save_figs__gt
        saveas(gcf,[ground_truth_summary_prefix,'_esz_by_networks'],'png')
    end
    

    % 3. Plot effect size vs. TPR
     
    figure
    hold on
    [~,ind]=sort(dcoeff);
    yyaxis left
    scatter(dcoeff,tpr,1,'b.')
    plot(dcoeff(ind),tpr_fit(ind),'k-','LineWidth',2)
    hold off
    
    % add stuff to TPR by esz
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
    set(gca,'fontsize',fontsz)
    % add trace of previous hist
    hold on
    yyaxis right
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
    plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/n_edges,'--','LineWidth',2)
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_tp],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off
   
    if save_figs__results
        saveas(gcf,esz_v_tpr_file,'png')
    end
    
    
    % 4. Plot effect size vs. TPR residuals - diagnostics

    figure;
    hold on;
    scatter(dcoeff,res,1,'b.')

    std_thresh=n_std_residual_outlier*std(res);
    idx=abs(res)>std_thresh;
    scatter(dcoeff(idx),res(idx),1,'.')
    plot(dcoeff,zeros(size(dcoeff)),'k-','LineWidth',2) % plot zero residual line

    hold off;
        
    if save_figs__results
        % save plot
        saveas(gcf,[summary_prefix,'_esz_v_TPR__residuals'],'png')
    end
    

    % 5. Plot effect size vs. TPR residuals - spatial distribution
    % edge-level results

    % put stuff back into upper triangle
    res_mat=zeros(n_nodes);
    res_mat(triu_msk)=res;

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
    fprintf(fid,'Percent less than d=%1.1f: %f; ALSO greater than %1.1f: %f',thresh_med,perc_edges_lt_thr_med,thresh_small,perc_edges_btw_thr_med_and_small);
        fprintf(fid,'\nAvg percent detected between d=+/-%1.1f: %f; ALSO greater than %1.1f: %f',thresh_med,tpr_lt_thr_med,thresh_small,tpr_btw_thr_med_and_small);
        fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_med,tpr_at_thr(1:3));
        fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_small,tpr_at_thr(4:6));
        fprintf(fid,'\n(%1.0f total permutations)',n_repetitions);
    fclose(fid);
end

