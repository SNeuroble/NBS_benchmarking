%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Task params

stat_type_gt='Size_Extent';

date_time_str_ground_truth.EMOTION_v_REST='10192020_1537'; % with new networks
date_time_str_ground_truth.GAMBLING_v_REST='10192020_1631';
date_time_str_ground_truth.LANGUAGE_v_REST='10192020_1642';
date_time_str_ground_truth.MOTOR_v_REST='10192020_1649';
date_time_str_ground_truth.RELATIONAL_v_REST='11112020_1345';
date_time_str_ground_truth.SOCIAL_v_REST='10192020_1656';
date_time_str_ground_truth.WM_v_REST='10192020_1714';
date_time_str_ground_truth.REST_v_REST2='03012020_1709';

% these are old ground truth without updated net/all summaries
% %date_time_str_ground_truth.EMOTION_v_REST='03012020_1722';
% date_time_str_ground_truth.EMOTION_v_REST='09222020_2337'; % with new networks
% date_time_str_ground_truth.GAMBLING_v_REST='03012020_1652';
% date_time_str_ground_truth.LANGUAGE_v_REST='03012020_1704';
% date_time_str_ground_truth.MOTOR_v_REST='03012020_1717';
% date_time_str_ground_truth.RELATIONAL_v_REST='03012020_1736';
% date_time_str_ground_truth.RELATIONAL_v_REST='09232020_0036';
% date_time_str_ground_truth.SOCIAL_v_REST='10192020_1656';
% % date_time_str_ground_truth.SOCIAL_v_REST='08052020_1304';
% %date_time_str_ground_truth.SOCIAL='03012020_1733';
% date_time_str_ground_truth.WM_v_REST='03012020_1709';
% date_time_str_ground_truth.REST_v_REST2='03012020_1709';

switch data_origin
case 'mrrc'
    switch stat_type
        case 'FDR'
            if grsize==80
                date_time_str_results.LANGUAGE_v_REST='12102020_1655';
                date_time_str_results.EMOTION_v_REST='03032020_1816'; % differences in edges too, plus sub nums
                date_time_str_results.GAMBLING_v_REST='03022020_1739';
            end
%             date_time_str_results.LANGUAGE_v_REST='02232020_0317';
%             date_time_str_results.MOTOR_v_REST='02292020_1853';
%             date_time_str_results.RELATIONAL_v_REST='02242020_1715';
%             date_time_str_results.SOCIAL_v_REST='03012020_1942';
%             date_time_str_results.WM_v_REST='02252020_1931';
    
        case 'Size_Extent'

            % IMPORTANT: these are all just copies of the originally named
            % summaries, just with _v_REST appended (e.g., EMOTION ->
            % EMOTION_v_REST)
            if grsize==40
                date_time_str_results.EMOTION_v_REST='03032020_1816'; % differences in edges too, plus sub nums
                date_time_str_results.GAMBLING_v_REST='03022020_1739';
                date_time_str_results.LANGUAGE_v_REST='02232020_0317';
                date_time_str_results.MOTOR_v_REST='02292020_1853';
                date_time_str_results.RELATIONAL_v_REST='02242020_1715';
                date_time_str_results.SOCIAL_v_REST='03012020_1942';
                date_time_str_results.WM_v_REST='02252020_1931';
                date_time_str_results.REST_v_REST2='08062020_0932';
            elseif grsize==80
                date_time_str_results.SOCIAL_v_REST='08032020_1807'; % no summary yet
            end

        case 'TFCE'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='03042020_0355';
                date_time_str_results.GAMBLING_v_REST='03032020_0007';
                date_time_str_results.LANGUAGE_v_REST='02242020_1327';
                date_time_str_results.MOTOR_v_REST='03012020_0101';
                date_time_str_results.RELATIONAL_v_REST='02252020_0829';
                date_time_str_results.SOCIAL_v_REST='03022020_0203';
                date_time_str_results.WM_v_REST='02262020_0139';
%                 date_time_str_results.REST_v_REST2='08062020_0519'; % I don't think this exists...
        
            elseif grsize==80
                date_time_str_results.SOCIAL_v_REST='08042020_0119';
            end

        case 'Constrained'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='03042020_0724';
                date_time_str_results.GAMBLING_v_REST='03032020_0332';
                date_time_str_results.LANGUAGE_v_REST='02242020_0355';
                date_time_str_results.MOTOR_v_REST='03012020_0417';
                date_time_str_results.RELATIONAL_v_REST='02252020_1519';
                date_time_str_results.SOCIAL_v_REST='03022020_0531';
                date_time_str_results.WM_v_REST='02262020_0457';
                date_time_str_results.REST_v_REST2='08052020_2143';
            elseif grsize==80
                date_time_str_results.SOCIAL_v_REST='08042020_0504';
            end

       case 'Omnibus'
           switch omnibus_type
               case 'Multidimensional_cNBS'
                  if grsize==80
                    date_time_str_results.EMOTION_v_REST='09112020_1539';
%                     date_time_str_results.GAMBLING='03032020_0332';
%                     date_time_str_results.LANGUAGE='02242020_0355';
%                     date_time_str_results.MOTOR='03012020_0417';
%                     date_time_str_results.RELATIONAL='02252020_1519';
%                     date_time_str_results.SOCIAL_v_REST='08042020_0504';
%                     %date_time_str_results.SOCIAL='03022020_0531';
%                     date_time_str_results.WM='02262020_0457';
%                     date_time_str_results.REST_v_REST2='08052020_2143';
                  end
                  
               case 'Multidimensional_all_edges'
                   if grsize==80
                       date_time_str_results.EMOTION_v_REST='09132020_0021';
                   end
               case 'Threshold_Both_Dir'
                   if grsize==80
                       date_time_str_results.EMOTION_v_REST='09102020_0244';
                   end
               otherwise
                    error('Omnibus type not defined');
           end
        case 'NA'
            % okay, trusting that won't need to set a benchmarking stat_type, e.g., not needed for running ground truth
        otherwise
            error('Stat type not defined');
    end
    
case 'farnam'
    switch stat_type
        case 'FDR'
                date_time_str_results.EMOTION_v_REST='testing_12232020_2201';
        case 'Size_Extent'
            if grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_0752';
                date_time_str_results.GAMBLING_v_REST='12182020_0617';
                date_time_str_results.LANGUAGE_v_REST='12202020_0453';
                date_time_str_results.MOTOR_v_REST='12202020_0348';
                date_time_str_results.RELATIONAL_v_REST='12202020_0356';
                date_time_str_results.SOCIAL_v_REST='12202020_0436';
                date_time_str_results.WM_v_REST='12202020_0436';
            end
        case 'TFCE'
            if grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_1731';
                date_time_str_results.GAMBLING_v_REST='12182020_1501';
		date_time_str_results.LANGUAGE_v_REST='12202020_1501';
                date_time_str_results.MOTOR_v_REST='12202020_1215';
                date_time_str_results.RELATIONAL_v_REST='12202020_1217';
                date_time_str_results.SOCIAL_v_REST='12202020_1308';
                date_time_str_results.WM_v_REST='12202020_1259';
            end
        case 'Constrained'
            if grsize==80
                date_time_str_results.EMOTION_v_REST='12252020_0030';
                date_time_str_results.GAMBLING_v_REST='12182020_2202';
            	date_time_str_results.LANGUAGE_v_REST='12202020_2342';
                date_time_str_results.MOTOR_v_REST='12202020_1848';
                date_time_str_results.RELATIONAL_v_REST='12202020_1845';
                date_time_str_results.SOCIAL_v_REST='12202020_1940';
                date_time_str_results.WM_v_REST='12202020_1928';
	    end
        case 'Omnibus'
           switch omnibus_type
               case 'Multidimensional_cNBS'
		    if grsize==80
			date_time_str_results.EMOTION_v_REST='12242020_2038';
                	date_time_str_results.GAMBLING_v_REST='12242020_2301';
			date_time_str_results.LANGUAGE_v_REST='12252020_0030';
                	date_time_str_results.MOTOR_v_REST='12242020_2323';
                	date_time_str_results.RELATIONAL_v_REST='12242020_2354';
                	date_time_str_results.SOCIAL_v_REST='12242020_2347';
                	date_time_str_results.WM_v_REST='12252020_0007';
		    end
	   end
    end
                    
end

%% Plot params

% descrip for levels of summary
pp.scaling_str{1}='_by_edges';
pp.scaling_str{2}='_by_networks';

% font
pp.fontsz=25;

% spline parameters
pp.window_sz{1}=0.01;
pp.spline_smoothing{1}=0.995;
pp.window_sz{2}=0.5;
pp.spline_smoothing{2}=0.999995;

% axis limits (used for histogram counting too)
pp.ax_ymin=0;
pp.ax_ymax_tp=100;
pp.ax_xmin_delta=-1.5; pp.ax_xmax_delta=1.5; pp.ax_ymax_esz_hist_delta=0.25; % special for delta
%   - edges
pp.ax_xmin{1}=-2; pp.ax_xmax{1}=2;
pp.ax_ymax_esz_hist{1}=0.1;
% ax_xmin=-2.5; ax_xmax=2.5;
% ax_ymin=0; ax_ymax_esz_hist=0.15; ax_ymax_tp=100; 
%   - nets
pp.ax_xmin{2}=-4; pp.ax_xmax{2}=4;
pp.ax_ymax_esz_hist{2}=0.5;

if exist('stat_type') % stat_type reserved for benchmarking and stat_type_gt for ground truth - sometimes both are used in one script
    if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
        pp.spline_smoothing{2}=0.99995;
        pp.window_sz{2}=0.2;
    end
end

% histograms params (keep an eye out for NAN/empty bins)
%   - edges
pp.bin_width{1}=0.07;
pp.bin_width_at_summary_thresh{1}=0.1;
pp.tpr_bin_width{1}=0.2;
%   - nets
pp.bin_width{2}=0.9;
pp.bin_width_at_summary_thresh{2}=1;
pp.tpr_bin_width{2}=3;


% effect size thresholds
pp.thresh_small=0.2; pp.thresh_med=0.5; pp.thresh_large=0.8;

% for visualizing residuals
pp.n_std_residual_outlier=2;

% color limits
pp.clim=[-pp.thresh_med, pp.thresh_med];
if exist('stat_type') % stat_type reserved for benchmarking and stat_type_gt for ground truth - sometimes both are used in one script
    if strcmp(stat_type,'Constrained') || strcmp(stat_type,'SEA')
        pp.clim_res{1}=[-60,60]; % for N=40
        pp.clim_res{2}=[-10,10]; % for N=40
    else
        %clim_res{2}=[-0.001,0.001]; % for N=20
        pp.clim_res{1}=[-3,3]; % for N=40
        pp.clim_res{2}=[-0.5,0.5]; % for N=40
    end
end

