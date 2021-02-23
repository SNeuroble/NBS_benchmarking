% set params for summarization - comparing all methods, with each method
% representing a combination of all tasks


%% Defaults

%default summary type
% summary_type_default='compare_all';

% default tasks, stats, and grsize
all_tasks={'EMOTION_v_REST','GAMBLING_v_REST','LANGUAGE_v_REST','MOTOR_v_REST','RELATIONAL_v_REST','SOCIAL_v_REST','WM_v_REST'};
all_stat_types={'Parametric_FDR','Size_Extent','TFCE','Constrained','Omnibus_Multidimensional_cNBS'};
grsize_default=80;

% saving defaults
save_settings.defaults.save_benchmarking_summary=1;
save_settings.defaults.save_dcoeff=1; % ground truth-specific
save_settings.defaults.save_figs=1;
save_settings.defaults.save_logs=1;
save_settings.defaults.save_summarized_data=1; % combined-specific % TODO: may have to partition for saving combined summary v combined visualization

% plotting defaults
% combine_all_tasks_default=0;
make_figs_default=1;
do_combined_default=0; % combined-specific

% combined-specific strings
task_combined='all_tasks'; % TODO: gotta rename for combined
stat_type_combined='all_stats'; % TODO: gotta rename

%% Mapping from stats to level of inference to ground truth level of pooling

% For combined: stat-to-level mapping for plots (first column can be partial strings - using pattern matching for first column)
stats_levelstr_map =    {'FDR','edge';
                        'Size','cluster (size)';
                        'TFCE','cluster (tfce)';
                        'Constrained','network';
                        'Omnibus','whole_brain'};
                    
% For combined: stat-to-ground truth level (first column can be partial strings - using pattern matching for first column)
statlevel_gtlevel_map =  {'edge', 'edge', 1;
                        'cluster', 'edge', 1;
                        'network', 'network', 2;
                        'whole_brain', 'whole_brain', 3};

% ground truth "stat_type" for filename (although it actually contains edge, network, and wb results)
stat_type_gt='Size_Extent';


%% Plot parameters (pp)
% indices correspond with ground truth levels (e.g., 1 = edge, 2=network, 3=whole_brain, often undefined)

% whether to combine all tasks to compare across stats
pp.do_combined=0;

% whether to use diag in making triu mask for spatial plots
pp.remove_matrix_diag.edge=1;
pp.remove_matrix_diag.network=0;
pp.remove_matrix_diag.whole_brain=0;

% bins for ground truth histograms
pp.dcoeff_hist_nbins.edge=40;
pp.dcoeff_hist_nbins.network=20;
pp.dcoeff_hist_nbins.whole_brain=2;

% for spline
pp.window_sz{1}=0.01;
pp.window_sz{2}=0.08; % pp.window_sz{2}=0.17; % 0.5 % 0.2
pp.spline_smoothing{1}=0.9999; %0.995;
pp.spline_smoothing{2}=0.99999; % 0.99995

% axis limits (used for histogram counting too)
pp.ax_ylim_tpr=[0,100];
pp.ax_xlim_dcoeff_v_tpr=[-1.5,1.5]; % for dcoeff v tpr, also works for edge-level dcoeff hist
pp.ax_xlim_dcoeff_hist=[-3.5,3.5];
pp.ax_ylim_dcoeff_hist=[0,0.15];
pp.ax_ylim_revcumhist=[0,1];
pp.ax_ylim_res=[-15,15];

% std thresholds visualizing residuals
pp.n_std_residual_outlier=2;

% color limits for spatial maps
% pp.clim=[-pp.thresh_med, pp.thresh_med];
pp.clim_dcoeff.edge=[-0.9,0.9];
pp.clim_dcoeff.network=[-2,2];
pp.clim_dcoeff.whole_brain=[-2,2];
% pp.clim_res=[-6,6];
pp.clim_res.edge=[-.5,.5];
pp.clim_res.network=[-.2,.2];
pp.clim_res.whole_brain=[-.2,.2];

% Color order for subsequent statistic types
pp.stat_color_order =   [ 0       0.4470    0.7410
                        0.8500    0.3250    0.0980
                        0.9290    0.6940    0.1250
                        0.4940    0.1840    0.5560
                        0.4660    0.6740    0.1880
                        0.3010    0.7450    0.9330
                        0.6350    0.0780    0.1840];



% atlas line width for plotting single tasks
pp.atlas_line_width_single_tasks=0.4;
                    
% font
pp.fontsz_sm=4;
pp.fontsz_lg=20;
pp.fontsz=pp.fontsz_lg; % set default to the larger
% pp.fontsz_single_task_scaling=0.2;

% figure sizes for saving (single tasks),lnih9
pp.fig_width_single_tasks=30;
pp.fig_height_single_tasks=3;
pp.fig_width_single_tasks_stats=20;
pp.fig_height_single_tasks_stats=30;
pp.fig_width_combined=30;
pp.fig_height_combined=20;
% pp.res=300; %300 dpi

% log params: binning and thresholding
pp.bin_width_at_summary_thresh{1}=0.1;
pp.bin_width_at_summary_thresh{2}=1;
pp.bin_width_at_summary_thresh{3}=1.5;
pp.tpr_bin_width.edge=0.2;
pp.tpr_bin_width.network=3; % TODO: 2
pp.tpr_bin_width.whole_brain=3;
pp.thresh_small=0.2; pp.thresh_med=0.5; pp.thresh_large=0.8; % effect size thresholds for summarizing tpr



% OLD

% pp.clim_single_task_scaling=0.5;
% pp.clim_res{1}=[-0.9,0.9];
% pp.clim_res{2}=[-2,2];
%         pp.clim_res{1}=[-60,60]; % for N=40
%         pp.clim_res{2}=[-10,10]; % for N=40
%         %clim_res{2}=[-0.001,0.001]; % for N=20
%         pp.clim_res{1}=[-3,3]; % for N=40 net
%         pp.clim_res{2}=[-0.5,0.5]; % for N=40 net


% % descrip for levels of summary
% pp.level_str{1}='_by_edges';
% pp.level_str{2}='_by_networks';
% pp.level_str{3}='_by_whole_brain';

% pp.ax_xmin_dcoeff{2}=-4; pp.ax_xmax_dcoeff{2}=4; % works for network-level dcoeff hist
% pp.ax_ymax_esz_hist{1}=0.1;
% pp.ax_ymax_esz_hist{2}=0.5;

% OLD histograms params (keep an eye out for NAN/empty bins)
% pp.bin_width{1}=0.07;
% pp.bin_width{2}=0.9;
% pp.bin_width{3}=1.5;

% pp.ax_xmin_delta=-1.5; pp.ax_xmax_delta=1.5;pp.ax_ymax_esz_hist_delta=0.25; % special for delta