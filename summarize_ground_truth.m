function summarize_ground_truth(task_type,stat_type,date_time_str,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ f
% e.g., summarize_ground_truth('LANGUAGE','Size_Extent','02102020_1759');
% Task can be: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION
% This script summarizes and visualizes ground truth effect sizes % TBD: and TPR data
% TBD: Summarization: bins edges by effect size, gets mean TPR within bins,
% and fits spline to effect size vs. mean TPR
% TBD: Plot: binned effect size, binned TPR,  d v. TPR spline, d v. TPR residual map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables

make_figs_default=0;
do_log_default=1;

p = inputParser;
addRequired(p,'task_type',@ischar);
addRequired(p,'stat_type',@ischar);
addRequired(p,'date_time_str',@ischar);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'do_log',do_log_default);
parse(p,task_type,stat_type,date_time_str,varargin{:});

task_type=p.Results.task_type;
stat_type=p.Results.stat_type;
date_time_str=p.Results.date_time_str;
make_figs=p.Results.make_figs;
do_log=p.Results.do_log;

%% User-defined plot params
ax_xmin=-2.5; ax_xmax=2.5; ax_ymin=0; ax_ymax_esz=0.15; ax_ymax_tp=100;
fontsz=25;
new_palette=[[45,170,180]/256; [255,166,0]/256];

% for histograms (keep an eye out for NAN/empty bins)
nbins=75;
bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);

% effect size thresholds
thresh_small=0.2; thresh_med=0.5; thresh_high=0.8;

% cmap threshold
clim=[-thresh_med, thresh_med];

%% Setup

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarkin
addpath(genpath(current_path));
setpaths;

% results file full path
results_prefix=[output_dir,'nbs_ground_truth__',task_type,'_',stat_type,'_',date_time_str];
results_filename=[results_prefix,'.mat'];

%% Load results and convert to d-coefficient
load(results_filename);
n_subs=size(UI.design.ui,2)-1;
dcoeff=edge_stats/sqrt(n_subs);

%% Calc percent edges within thresholds
n_edges=length(dcoeff);
perc_edges_lt_thr_med=sum(+abs(dcoeff)<thresh_med) * 100 / n_edges;
perc_edges_btw_thr_med_and_small=sum(+(abs(dcoeff)<thresh_med & abs(dcoeff)>thresh_small)) * 100 / n_edges;

%% Summarize for interpretation
n_nodes=size(cluster_stats,1);
if make_figs
    
    % plot
    h=histogram(dcoeff,bin_edges,'Normalization','probability');
    hold on;
    plot(h.BinEdges(1:end-1) + h.BinWidth/2, h.BinCounts/length(dcoeff))
    hold off;
    
    % add stuff to esz hist
    adorn_plot=1;
    if adorn_plot
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
        set(gca,'fontsize',fontsz)
        
        % highlight
        hold on
        rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
        hold off
    end
    saveas(gcf,[results_prefix,'_esz_hist_test'],'png')

    % put stuff back into upper triangle 
    msk=triu(logical(ones(n_nodes)),1);
    dcoeff_mat=zeros(n_nodes);
    dcoeff_mat(msk)=dcoeff;
    
    % network-level results
    draw_atlas_boundaries(dcoeff_mat');
    colormap(bipolar([],0.1));
    caxis(clim);
    saveas(gcf,[results_prefix,'_esz_edges_test'],'png')
    
    summarize_matrix_by_atlas(dcoeff_mat');
    colormap(bipolar([],0.1));
    caxis(clim);
    saveas(gcf,[results_prefix,'_esz_by_networks_test'],'png')

end

%% Log Summaries: perc esz and TP at thresholds
if do_log
    logfile=[results_prefix,'_log'];
    fid=fopen(logfile,'w');
    fprintf(fid,'Percent less than d=%1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,perc_edges_lt_thr_med,thresh_small,perc_edges_btw_thr_med_and_small);
%     fprintf(fid,'\nAvg percent detected between d=+/-%1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,tpr_lt_thr_high,thresh_low,tpr_btw_thr_high_and_low);
%     fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_high,tpr_at_thr(1:3));
%     fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_low,tpr_at_thr(4:6));
%     fprintf(fid,'\n(%1.0f total permutations)',rep);
    fclose(fid);
end

