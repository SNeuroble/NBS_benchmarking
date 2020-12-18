% set params for meta-summary - comparing all methods, with each method
% representing a combination of all tasks

%% Task params

% stat_type_gt='Size_Extent'; % TBD

% date_time_str_ground_truth_combinedtasks.FDR='08062020_0932'; % edges
date_time_str_ground_truth_combined.Size_Extent='12022020'; % edges
date_time_str_ground_truth_combined.TFCE='12022020'; % edges
date_time_str_ground_truth_combined.Constrained='12022020'; % networks
% date_time_str_ground_truth_combinedtasks.Omnibus='08062020_0932'; % omnibus A


switch grsize
    case 40
%         date_time_str_combined.FDR='08062020_0932';
        date_time_str_combined.Size_Extent='12022020';
        date_time_str_combined.TFCE='12022020';
        date_time_str_combined.Constrained='12022020';
%         date_time_str_combined.Omnibus='08062020_0932';
    case 80
%         date_time_str_combinedtasks.FDR='08062020_0932';
%         date_time_str_combinedtasks.Size_Extent='03032020_1816';
%         date_time_str_combinedtasks.TFCE='08062020_0932';
%         date_time_str_combinedtasks.Constrained='08062020_0932';
%         date_time_str_combinedtasks.Omnibus='08062020_0932';
    case 120
%         date_time_str_combinedtasks.FDR='08062020_0932';
%         date_time_str_combinedtasks.Size_Extent='03032020_1816';
%         date_time_str_combinedtasks.TFCE='08062020_0932';
%         date_time_str_combinedtasks.Constrained='08062020_0932';
%         date_time_str_combinedtasks.Omnibus='08062020_0932';
    otherwise
        error('Group size not defined');
end

%% Plot params

% descrip for levels of summary
% pp.scaling_str{1}='_by_edges';
% pp.scaling_str{2}='_by_networks';

% font
pp.fontsz=20;

% spline parameters
pp.window_sz{1}=0.01;
pp.spline_smoothing{1}=0.995;
pp.window_sz{2}=0.17;
pp.spline_smoothing{2}=0.999995;

% axis limits (used for histogram counting too)
pp.ax_ymin=0;
pp.ax_ymax_tp=100;
% pp.ax_xmin_delta=-1.5; pp.ax_xmax_delta=1.5; pp.ax_ymax_esz_hist_delta=0.25; % special for delta
% %   - edges
pp.ax_xmin=-1.5; pp.ax_xmax=1.5;
% pp.ax_xmin{1}=-2; pp.ax_xmax{1}=2;
% pp.ax_ymax_esz_hist{1}=0.1;
% % ax_xmin=-2.5; ax_xmax=2.5;
% % ax_ymin=0; ax_ymax_esz_hist=0.15; ax_ymax_tp=100; 
% %   - nets
% pp.ax_xmin{2}=-4; pp.ax_xmax{2}=4;
% pp.ax_ymax_esz_hist{2}=0.5;
% 
% if exist('stat_type') % stat_type reserved for benchmarking and stat_type_gt for ground truth - sometimes both are used in one script
%     if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
%         pp.spline_smoothing{2}=0.99995;
%         pp.window_sz{2}=0.2;
%     end
% end
% 
% % histograms params (keep an eye out for NAN/empty bins)
% %   - edges
pp.tpr_bin_width=2;
% % pp.bin_width{1}=0.07;
% % pp.bin_width_at_summary_thresh{1}=0.1;
% % pp.tpr_bin_width{1}=0.2;
% %   - nets
% % pp.bin_width{2}=0.9;
% % pp.bin_width_at_summary_thresh{2}=1;
% % pp.tpr_bin_width{2}=3;
% 
% 
% % effect size thresholds
% pp.thresh_small=0.2; pp.thresh_med=0.5; pp.thresh_large=0.8;
% 
% % for visualizing residuals
% pp.n_std_residual_outlier=2;
% 
% % color limits
% pp.clim=[-pp.thresh_med, pp.thresh_med];
% if exist('stat_type') % stat_type reserved for benchmarking and stat_type_gt for ground truth - sometimes both are used in one script
%     if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
%         pp.clim_res{1}=[-60,60]; % for N=40
%         pp.clim_res{2}=[-10,10]; % for N=40
%     else
%         %clim_res{2}=[-0.001,0.001]; % for N=20
%         pp.clim_res{1}=[-3,3]; % for N=40
%         pp.clim_res{2}=[-0.5,0.5]; % for N=40
%     end
% end

% colors
stat_color_order=[ 0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

