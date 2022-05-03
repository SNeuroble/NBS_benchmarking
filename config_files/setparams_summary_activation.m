%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% User-defined parameters for summarization & comparing all methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%% Defaults %%%%%%%%%

% Default: tasks, stats, and grsize
all_tasks={'EMOTION_v_REST','GAMBLING_v_REST','LANGUAGE_v_REST','MOTOR_v_REST','RELATIONAL_v_REST','SOCIAL_v_REST','WM_v_REST'};
all_stat_types={'Parametric_Bonferroni','Parametric_FDR','Size_Extent','TFCE','Constrained_FWER','Constrained','Omnibus_Multidimensional_cNBS'};
grsize_default=80;
use_preaveraged_constrained=1;

% Default: saving settings
save_settings.defaults.save_benchmarking_summary=1;
save_settings.defaults.save_dcoeff=1; % ground truth-specific
save_settings.defaults.save_figs=1;
save_settings.defaults.save_logs=1;
save_settings.defaults.save_summarized_data=1; % combined-specific 

% Default: plotting params
make_figs_default=1;
do_combined_default=0; % combined versus task-specific summaries

% Default: combined summary strings
task_combined='all_tasks'; % TODO: may rename for combined
stat_type_combined='all_stats'; % TODO: may rename

% Default: FPR-specific strings
if do_fpr
    fpr_str='_shuffled_for_FPR_';
else
    fpr_str='';
end




%%%%%% Mapping: stats -> level of inference -> ground truth level %%%%%%%%

% For combined: stat -> level mapping for plots 
% first column can be partial strings - using pattern matching for that col
stats_levelstr_map =    {'Parametric_Bonferroni','edge';
                        'Parametric_FDR','edge (fdr)';
                        'Size_Extent','cluster size';
                        'TFCE','cluster tfce';
                        'Constrained_FWER','network';
                        'Constrained','network (fdr)';
                        'Omnibus_Multidimensional_cNBS','whole brain'};
                    
% For combined: stat -> ground truth level
% first column can be partial strings - using pattern matching for that col
statlevel_gtlevel_map =  {'edge', 'edge', 1;
                        'cluster', 'edge', 1;
                        'network', 'network', 2;
                        'whole', 'whole_brain', 3};

% Ground truth "stat_type" for filename
% although it actually contains edge, network, and wb results
stat_type_gt='Size_Extent';





%%%%%%%%%%%%%%% Plot parameters (pp) %%%%%%%%%%%%%%%

% indices correspond with ground truth levels (e.g., 1 = edge, 2=network, 3=whole_brain, often undef)

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
pp.window_sz{2}=0.08; 
pp.spline_smoothing{1}=0.9999; 
pp.spline_smoothing{2}=0.99999; 

% axis limits (used for histogram counting too)
pp.ax_ylim_tpr=[0,100];
pp.ax_xlim_dcoeff_v_tpr=[-1.1,1.1]; % for dcoeff v tpr, also works for edge-level dcoeff hist
pp.ax_xlim_dcoeff_hist=[-3.5,3.5];
pp.ax_ylim_dcoeff_hist=[0,0.15];
pp.ax_ylim_revcumhist=[0,1];
pp.ax_ylim_res.edge=[-7,7];
pp.ax_ylim_res.network=[-15,15];
pp.ax_ylim_lp=[0,100];
pp.ax_ylim_FWERstrong=[0,20];
pp.ax_ylim_FWERweak=[0,7];
pp.ax_ylim_FDR=[0,5];
pp.ax_ylim_num_fp=[0,100];
pp.ax_ylim_spatial_extent_fp=[0,5];

% std thresholds visualizing residuals
pp.n_std_residual_outlier=2;

% color limits for spatial maps
% pp.clim=[-pp.thresh_med, pp.thresh_med];
pp.clim_dcoeff.edge=[-0.5,0.5];
pp.clim_dcoeff.network=[-1.5,1.5];
pp.clim_dcoeff.whole_brain=[-1.5,1.5];
% pp.clim_res=[-6,6];
pp.clim_res.edge=[-.6,.6];
pp.clim_res.edge_bynet=[-.15,.15];
pp.clim_res.network=[-.6,.6];
pp.clim_res.whole_brain=[-.15,.15];


% Color palette for each statistic type
% colorful nested
pp.stat_color_order =   [235, 215, 28
                        254, 167, 0
                        31, 192, 198
                        26, 99, 102
                        254, 88, 118
                        188, 78, 144
                        30, 60, 110];
                    

                   
% Alternative color palettes
%{
pp.stat_color_order =   [ 0       0.4470    0.7410
                         0.8500    0.3250    0.0980
                         0.9290    0.6940    0.1250
                         0.4940    0.1840    0.5560
                         0.4660    0.6740    0.1880
                         0.3010    0.7450    0.9330
                         0.6350    0.0780    0.1840];


% colorful nested
pp.stat_color_order =   [235, 215, 28
                        254, 167, 0
                        31, 192, 198
                        26, 99, 102
                        254, 88, 118
                        188, 78, 144
                        30, 60, 110];
                    


% colorful nested
pp.stat_color_order =   [235, 215, 28
                        254, 167, 0
                        31, 192, 198
                        26, 99, 102
                        254, 88, 118
                        188, 78, 144
                        30, 60, 110];
                    

% orig
pp.stat_color_order =   [235, 215, 28
                        255, 166, 0
                        26, 192, 198
                        188, 80, 144
                        255, 88, 77
                        111, 77, 168
                        0, 63, 92];

%muted
pp.stat_color_order =   [252, 234, 165
                        251, 198, 40
                        186, 245, 245
                        31, 192, 198
                        236, 205, 228
                        188, 78, 146
                        20, 53, 77];
% modern             
pp.stat_color_order =   [185, 169, 233
                        106, 79, 163
                        186, 245, 245
                        31, 192, 198
                        254, 190, 186
                        254, 88, 77
                        20, 53, 77];

% orig ordered
pp.stat_color_order =   [188, 80, 144
                        255, 88, 77
                        255, 166, 0
                        235, 215, 28
                        26, 192, 198
                        111, 77, 168
                        0, 63, 92];

% soft
pp.stat_color_order =   [186, 245, 245
                        11, 174, 190
                        246, 201, 219
                        222, 97, 119
                        252, 223, 188
                        242, 170, 60
                        116, 84, 145];
%}

pp.stat_color_order=pp.stat_color_order/255;

% FPR colors/color limits
pp.clim_fpr.edge=[0,1];
pp.clim_fpr.edge_bynet=[0,0.04];
pp.clim_fpr.network=[0,2];
pp.cmap_fpr =   [0 0 0
                1 0 0];
pp.cmap_fpr = interp1([1,size(pp.cmap_fpr,1)],pp.cmap_fpr, 1:0.01:2,'linear');

% line type for FDR
pp.line_type_fdr='--';

% atlas line width for plotting single tasks
pp.atlas_line_width_single_tasks=0.4;
                    
% font
pp.fontsz_sm=5;
pp.fontsz_lg=18;
pp.fontsz=pp.fontsz_lg; % set default to the larger

% figure sizes for saving (single tasks)
pp.margins=[0.015,0.015];
pp.margins_bigx=[0.04,0.04];
pp.margins_dcoeff=[0.02,0.02];
pp.margins_bigx_dcoeff=[0.09,0.09];
pp.fig_width_single_tasks_dcoeff=20;
pp.fig_height_single_tasks_dcoeff=7;
pp.fig_width_combined_dcoeff=12;
pp.fig_height_combined_dcoeff=10;

pp.fig_sz_1plot=7;
pp.fig_sz_3plots=20;
pp.fig_sz_7plots=20;

pp.fig_height_1plt_norm=0.21;
pp.fig_width_1plt_norm=0.14;

pp.fig_height_3plts_norm=pp.fig_height_1plt_norm*3;
pp.fig_height_4plts_norm=pp.fig_height_1plt_norm*4;
pp.fig_width_3plts_norm=pp.fig_width_1plt_norm*3;
pp.fig_width_4plts_norm=pp.fig_width_1plt_norm*4;
pp.fig_width_6plts_norm=pp.fig_width_1plt_norm*6;
pp.fig_width_7plts_norm=pp.fig_width_1plt_norm*7;

pp.fig_height_1plt_norm_big=pp.fig_height_1plt_norm*2.5;
pp.fig_width_1plt_norm_big=pp.fig_width_1plt_norm*3;


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
pp.tpr_bin_width.network=3;
pp.tpr_bin_width.whole_brain=3;
pp.thresh_small=0.2; pp.thresh_med=0.5; pp.thresh_large=0.8; % effect size thresholds for summarizing tpr

