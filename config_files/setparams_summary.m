% set params for summary

%% Task params

stat_type_gt='Size_Extent';

date_time_str_ground_truth.EMOTION='03012020_1722';
date_time_str_ground_truth.GAMBLING='03012020_1652';
date_time_str_ground_truth.LANGUAGE='03012020_1704';
date_time_str_ground_truth.MOTOR='03012020_1717';
date_time_str_ground_truth.RELATIONAL='03012020_1736';
date_time_str_ground_truth.SOCIAL_v_REST='08052020_1304';
%date_time_str_ground_truth.SOCIAL='03012020_1733';
date_time_str_ground_truth.WM='03012020_1709';
date_time_str_ground_truth.REST_v_REST2='03012020_1709';

switch stat_type
    case 'Size_Extent'
        date_time_str_results.EMOTION='03032020_1816';
        date_time_str_results.GAMBLING='03022020_1739';
        date_time_str_results.LANGUAGE='02232020_0317';
        date_time_str_results.MOTOR='02292020_1853';
        date_time_str_results.RELATIONAL='02242020_1715';
        date_time_str_results.SOCIAL_v_REST='08032020_1807';
        %date_time_str_results.SOCIAL='03012020_1942';
        date_time_str_results.WM='02252020_1931';
        date_time_str_results.REST_v_REST2='08062020_0932';
    case 'TFCE'
        date_time_str_results.EMOTION='03042020_0355';
        date_time_str_results.GAMBLING='03032020_0007';
        date_time_str_results.LANGUAGE='02242020_1327';
        date_time_str_results.MOTOR='03012020_0101';
        date_time_str_results.RELATIONAL='02252020_0829';
        date_time_str_results.SOCIAL_v_REST='08042020_0119';
        %date_time_str_results.SOCIAL='03022020_0203';
        date_time_str_results.WM='02262020_0139';
        date_time_str_results.REST_v_REST2='08062020_0519';
    case 'Constrained'
        date_time_str_results.EMOTION='03042020_0724';
        date_time_str_results.GAMBLING='03032020_0332';
        date_time_str_results.LANGUAGE='02242020_0355';
        date_time_str_results.MOTOR='03012020_0417';
        date_time_str_results.RELATIONAL='02252020_1519';
        date_time_str_results.SOCIAL_v_REST='08042020_0504';
        %date_time_str_results.SOCIAL='03022020_0531';
        date_time_str_results.WM='02262020_0457';
        date_time_str_results.REST_v_REST2='08052020_2143';
    case 'NA'
        % okay, trusting that won't need to set a benchmarking stat_type, e.g., not needed for running ground truth
    otherwise
        error('Stat type not defined');
end

%% Plot params

% x and y axis limits for esz and tpr plots
ax_xmin=-2; ax_xmax=2;
ax_ymin=0; ax_ymax_esz=0.1; ax_ymax_tp=100;
ax_xmin_delta=-1.5; ax_xmax_delta=1.5; ax_ymax_esz_delta=0.25; % special for delta
% ax_xmin=-2.5; ax_xmax=2.5;
% ax_ymin=0; ax_ymax_esz=0.15; ax_ymax_tp=100; 

% font
fontsz=25;

% smoothing factor for spline
spline_smoothing=0.995;
spline_smoothing_set=0.99995;

% histograms params (keep an eye out for NAN/empty bins)
bin_width=0.07;
nbins=ceil((ax_xmax-ax_xmin)/bin_width);
% nbins=75; % this fits for the narrower x-axis, 60 for larger axis
bin_width_at_summary_thresh=0.1;
bin_width_at_summary_thresh__network=0.1;
% half_bin_width=bin_width/2; % ad hoc bin size for 

% effect size thresholds
thresh_small=0.2; thresh_med=0.5; thresh_large=0.8;

% for visualizing residuals
n_std_residual_outlier=2;

% color limits

clim=[-thresh_med, thresh_med];

if exist('stat_type') % stat_type reserved for benchmarking and stat_type_gt for ground truth - sometimes both are used in one script
    if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
        clim_res=[-10,10]; % for N=40
        clim_res_detailed=[-60,60]; % for N=40
    else
        %clim_res=[-0.001,0.001]; % for N=20
        clim_res=[-0.5,0.5]; % for N=40
        clim_res_detailed=[-3,3]; % for N=40
    end
end

clim_all_tasks=[-0.3,0.3]; % for N=40
clim_all_tasks_detailed=[-0.5,0.5]; % for N=40
