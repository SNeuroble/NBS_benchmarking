% set date_time_strs based on updated list

%% Now (for saving comparison across stats)

date_time_str_now=datestr(now,'mmddyyyy');

%% Ground truth

date_time_str_ground_truth.EMOTION_v_REST='02222021_1312';
date_time_str_ground_truth.GAMBLING_v_REST='02222021_1645';
date_time_str_ground_truth.LANGUAGE_v_REST='02222021_1335';
date_time_str_ground_truth.MOTOR_v_REST='02222021_1705';
date_time_str_ground_truth.RELATIONAL_v_REST='02222021_1728';
date_time_str_ground_truth.SOCIAL_v_REST='02222021_1751';
date_time_str_ground_truth.WM_v_REST='02222021_1820';
% date_time_str_ground_truth.EMOTION_v_REST='10192020_1537'; % with new networks
% date_time_str_ground_truth.GAMBLING_v_REST='10192020_1631';
% date_time_str_ground_truth.LANGUAGE_v_REST='10192020_1642';
% date_time_str_ground_truth.MOTOR_v_REST='10192020_1649';
% date_time_str_ground_truth.RELATIONAL_v_REST='11112020_1345';
% date_time_str_ground_truth.SOCIAL_v_REST='10192020_1656';
% date_time_str_ground_truth.WM_v_REST='10192020_1714';
% date_time_str_ground_truth.REST_v_REST2='03012020_1709';

%% Single task summary filename info



switch data_origin
case 'farnam'
    
switch summary_type
case {'visualize_tpr','visualize_gt','dcoeff'}
    
    switch grsize
        case 40
%             date_time_str_combined.Parametric_FDR='01312021';
%             date_time_str_combined.Size_Extent='12022020';
%             date_time_str_combined.TFCE='12022020';
%             date_time_str_combined.Constrained='12022020';
%             date_time_str_combined.Omnibus_Multidimensional_cNBS='01192021';
%             date_time_str_combined='02182021';
            date_time_str_combined='02222021';
        case 80
%             date_time_str_combined.FDR='01152021';
%             date_time_str_combined.Parametric_FDR='01312021';
%             date_time_str_combined.Size_Extent='01072021';
%             date_time_str_combined.TFCE='01072021';
%             date_time_str_combined.Constrained='01072021';
%             date_time_str_combined.Omnibus_Multidimensional_cNBS='01182021';
%             date_time_str_combined='02102021'; %02112021
            date_time_str_combined='02222021'; %'02172021'; % '02152021';'02102021'; %02112021
        case 120
%             date_time_str_combined.Parametric_FDR='01312021';
%             date_time_str_combined.Size_Extent='01312021';
%             date_time_str_combined.TFCE='01312021';
%             date_time_str_combined.Constrained='01312021';
%             date_time_str_combined.Omnibus_Multidimensional_cNBS='01312021';
            date_time_str_combined='02222021'; % '02182021';
        otherwise
            error('Group size not defined');
    end
    
    
    
otherwise
    switch stat_type
        case 'Parametric_FDR'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='01312021_0053';
                date_time_str_results.GAMBLING_v_REST='01312021_0100';
                date_time_str_results.LANGUAGE_v_REST='01312021_0103';
                date_time_str_results.MOTOR_v_REST='01312021_0113';
                date_time_str_results.RELATIONAL_v_REST='01312021_0116';
                date_time_str_results.SOCIAL_v_REST='01312021_0117';
                date_time_str_results.WM_v_REST='01312021_0120';
                date_time_str_results.REST_v_REST2='01312021_0047';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='01312021_0138';
                date_time_str_results.GAMBLING_v_REST='01312021_0141';
                date_time_str_results.LANGUAGE_v_REST='01312021_0145';
                date_time_str_results.MOTOR_v_REST='01312021_0146';
                date_time_str_results.RELATIONAL_v_REST='01312021_0152';
                date_time_str_results.SOCIAL_v_REST='01312021_0151';
                date_time_str_results.WM_v_REST='01312021_0153';
                date_time_str_results.REST_v_REST2='01312021_0057';
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01312021_0205';
                date_time_str_results.GAMBLING_v_REST='01312021_0216';
                date_time_str_results.LANGUAGE_v_REST='01312021_0240';
                date_time_str_results.MOTOR_v_REST='01312021_0239';
                date_time_str_results.RELATIONAL_v_REST='01312021_0237';
                date_time_str_results.SOCIAL_v_REST='01312021_0242';
                date_time_str_results.WM_v_REST='01312021_0240';
                date_time_str_results.REST_v_REST2='01312021_0120';
            end
        case 'FDR' % not using anymore bc invalid FWER control
            if grsize==40
                date_time_str_results.EMOTION_v_REST='01082021_1335';
                date_time_str_results.GAMBLING_v_REST='01062021_1402';
                date_time_str_results.LANGUAGE_v_REST='01082021_1641';
                date_time_str_results.MOTOR_v_REST='01082021_1519';
                date_time_str_results.RELATIONAL_v_REST='01082021_1722';
                date_time_str_results.SOCIAL_v_REST='01082021_1257';
                date_time_str_results.WM_v_REST='01082021_1655';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12062020_0521_corrected';
                date_time_str_results.GAMBLING_v_REST='01142021_0100';
                date_time_str_results.LANGUAGE_v_REST='12102020_1655_corrected';
                date_time_str_results.MOTOR_v_REST='01132021_1953';
                date_time_str_results.RELATIONAL_v_REST='01092021_0610';
                date_time_str_results.SOCIAL_v_REST='01092021_0151';
                date_time_str_results.WM_v_REST='01132021_2131';
            end
        case 'Size_Extent'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='02182021_1726';
%                 date_time_str_results.EMOTION_v_REST='03032020_1816'; % differences in edges too, plus sub nums
                date_time_str_results.GAMBLING_v_REST='03022020_1739';
                date_time_str_results.LANGUAGE_v_REST='02182021_1752';
%                 date_time_str_results.LANGUAGE_v_REST='02232020_0317';
                date_time_str_results.MOTOR_v_REST='02292020_1853';
                date_time_str_results.RELATIONAL_v_REST='02242020_1715';
                date_time_str_results.SOCIAL_v_REST='03012020_1942';
                date_time_str_results.WM_v_REST='02252020_1931';
                date_time_str_results.REST_v_REST2='08062020_0932';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_0752';
                date_time_str_results.GAMBLING_v_REST='12182020_0617';
                date_time_str_results.LANGUAGE_v_REST='12202020_0453';
                date_time_str_results.MOTOR_v_REST='12202020_0348';
                date_time_str_results.RELATIONAL_v_REST='12202020_0356';
                date_time_str_results.SOCIAL_v_REST='12202020_0436';
                date_time_str_results.WM_v_REST='12202020_0436';
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01262021_0007';
                date_time_str_results.GAMBLING_v_REST='01282021_0025';
                date_time_str_results.LANGUAGE_v_REST='01272021_2300';
                date_time_str_results.MOTOR_v_REST='01272021_2304';
                date_time_str_results.RELATIONAL_v_REST='01272021_2256';
                date_time_str_results.SOCIAL_v_REST='01282021_0204';
                date_time_str_results.WM_v_REST='01282021_0209';
                date_time_str_results.REST_v_REST2='01282021_0826';
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
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_1731';
                date_time_str_results.GAMBLING_v_REST='12182020_1501';
                date_time_str_results.LANGUAGE_v_REST='12202020_1501';
                date_time_str_results.MOTOR_v_REST='12202020_1215';
                date_time_str_results.RELATIONAL_v_REST='12202020_1217';
                date_time_str_results.SOCIAL_v_REST='12202020_1308';
                date_time_str_results.WM_v_REST='12202020_1259';
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01272021_1155';
                date_time_str_results.GAMBLING_v_REST='01272021_1153';
                date_time_str_results.LANGUAGE_v_REST='01272021_1134';
                date_time_str_results.MOTOR_v_REST='01292021_1204';
                date_time_str_results.RELATIONAL_v_REST='01272021_1207';
                date_time_str_results.SOCIAL_v_REST='01272021_1234';
                date_time_str_results.WM_v_REST='01272021_1216';
                date_time_str_results.REST_v_REST2='01282021_2336';
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
                date_time_str_results.EMOTION_v_REST='12252020_0030';
                date_time_str_results.GAMBLING_v_REST='12182020_2202';
            	date_time_str_results.LANGUAGE_v_REST='12202020_2342';
                date_time_str_results.MOTOR_v_REST='12202020_1848';
                date_time_str_results.RELATIONAL_v_REST='12202020_1845';
                date_time_str_results.SOCIAL_v_REST='12202020_1940';
                date_time_str_results.WM_v_REST='12202020_1928';
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01262021_0053';
                date_time_str_results.GAMBLING_v_REST='01272021_2332';
                date_time_str_results.LANGUAGE_v_REST='01272021_2357';
                date_time_str_results.MOTOR_v_REST='01282021_0023';
                date_time_str_results.RELATIONAL_v_REST='01282021_0055';
                date_time_str_results.SOCIAL_v_REST='01272021_2244';
                date_time_str_results.WM_v_REST='01272021_2310';
                date_time_str_results.REST_v_REST2='i01282021_2135';
            end
        case 'Omnibus_Multidimensional_cNBS'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='01192021_0119';
                date_time_str_results.GAMBLING_v_REST='01192021_0158';
                date_time_str_results.LANGUAGE_v_REST='01192021_0215';
                date_time_str_results.MOTOR_v_REST='01192021_0157';
                date_time_str_results.RELATIONAL_v_REST='01192021_0228';
                date_time_str_results.SOCIAL_v_REST='01192021_0213';
                date_time_str_results.WM_v_REST='01192021_0246';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12242020_2038';
                date_time_str_results.GAMBLING_v_REST='12242020_2301';
                date_time_str_results.LANGUAGE_v_REST='12252020_0030';
                date_time_str_results.MOTOR_v_REST='12242020_2323';
                date_time_str_results.RELATIONAL_v_REST='12242020_2354';
                date_time_str_results.SOCIAL_v_REST='12242020_2347';
                date_time_str_results.WM_v_REST='12252020_0007';
             elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01272021_1104';
                date_time_str_results.GAMBLING_v_REST='01272021_1008';
                date_time_str_results.LANGUAGE_v_REST='01272021_1028';
                date_time_str_results.MOTOR_v_REST='01272021_1008';
                date_time_str_results.RELATIONAL_v_REST='01272021_1036';
                date_time_str_results.SOCIAL_v_REST='01272021_1156';
                date_time_str_results.WM_v_REST='01272021_2126';
                date_time_str_results.REST_v_REST2='01282021_1956';
            end
        otherwise; error('Statistic type does not exist or is not specified.')
    end
    end
otherwise; error('Data source does not exist or is not specified. Right now, only one source is permitted (farnam).')
end

%% Set io filenames

if ~contains(summary_type,'visualize') % visualizations rely on combined summary, not these individual task files
    
    % ground truth variables
    ground_truth_results_basename_prefix=['ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
    ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
    ground_truth_dcoeff_filename=[output_dir,ground_truth_results_basename_prefix,'_dcoeff.mat'];

    if strcmp(summary_type,'calculate_tpr')
        % original (un-summarized) benchmarking variables
        benchmarking_results_basename_prefix=['results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
        results_filename=[output_dir,benchmarking_results_basename_prefix,'.mat'];

        % summarized benchmarking variables (output 1)
        benchmarking_summary_filename=[output_dir,benchmarking_results_basename_prefix,'_summary.mat'];

        % summary figs/logs and combined summary vars
        summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];

        % create new timestamp for a new file
        date_time_str_combined=date_time_str_now;
    end

% else
    % combined_summary_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_combined.(task)];
    % method_comparison_basename_prefix=[summary_output_dir,'results__',task,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_combined.(stat_type)];

    % TODO: rename as follows - lots of reps here
    % summary_basename_prefix=summary_prefix;
    % combined_basename_prefix=method_comparison_basename_prefix; 
    % combined_summary_filename=[combined_basename_prefix,'_summary.mat'];
    % TODO: consider stat_type='full_comparison' so can reuse the summary_prefix above
    
    % use timestamp from already-created file
end

% combined summary filename
combined_summary_dir=[output_dir,'combined_summary/'];
combined_basename_prefix=['combined_grsize',num2str(grsize),'_',date_time_str_combined];
combined_filename_prefix=[combined_summary_dir,combined_basename_prefix];
combined_summary_filename=[combined_filename_prefix,'_summary.mat'];

% ground truth filename
% name with combined prefix bc using combined data
ground_truth_vis_dir=[combined_summary_dir,'ground_truth/'];
ground_truth_vis_filename_prefix=[ground_truth_vis_dir,'ground_truth_',date_time_str_combined];

% individual task summaries filename
% name with combined prefix bc using combined data
combined_by_task_dir=[combined_summary_dir,'individual_task_summary/'];
% combined_by_task_basename_prefix=['tasks_grsize',num2str(grsize),'_',date_time_str_combined];
combined_by_task_filename_prefix=[combined_by_task_dir,combined_basename_prefix];

% log filenames
% name with combined prefix bc using combined data
log_dir=[combined_by_task_dir,'logs/'];
log_filename_prefix=[log_dir,combined_basename_prefix];

% set edge groups file in case can't get from original data
edge_groups_filename=[output_dir,'edge_groups.mat'];


