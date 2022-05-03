%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Set date/time strings and resulting filenames corresponding with 
% input/output data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%% Define now (for saving summaries) %%%%%%%%

date_time_str_now=datestr(now,'mmddyyyy');




%%%%%%%% Ground truth %%%%%%%%

% Task versus rest
date_time_str_ground_truth.EMOTION_v_REST='02222021_1312';
date_time_str_ground_truth.GAMBLING_v_REST='02222021_1645';
date_time_str_ground_truth.LANGUAGE_v_REST='02222021_1335';
date_time_str_ground_truth.MOTOR_v_REST='02222021_1705';
date_time_str_ground_truth.RELATIONAL_v_REST='02222021_1728';
date_time_str_ground_truth.SOCIAL_v_REST='02222021_1751';
date_time_str_ground_truth.WM_v_REST='02222021_1820';
date_time_str_ground_truth.REST_v_REST2='03012020_1709';

% Task versus rest, using trimmed resting runs (supplementary)
date_time_str_ground_truth.EMOTION_v_REST_176frames='04012021_0131';
date_time_str_ground_truth.GAMBLING_v_REST_253frames='04012021_1353';
date_time_str_ground_truth.LANGUAGE_v_REST_316frames='04012021_1412';
date_time_str_ground_truth.MOTOR_v_REST_284frames='04012021_1426';
date_time_str_ground_truth.RELATIONAL_v_REST_232frames='04012021_1442';
date_time_str_ground_truth.SOCIAL_v_REST_274frames='04012021_1503';
date_time_str_ground_truth.WM_v_REST_405frames='04012021_1531';

% Task singles (no contrast) (supplementary)
date_time_str_ground_truth.REST_176frames='04012021_1819';
date_time_str_ground_truth.REST_232frames='04012021_1711';
date_time_str_ground_truth.REST='03192021_0130';
date_time_str_ground_truth.EMOTION='03182021_2251';
date_time_str_ground_truth.GAMBLING='03192021_1359';
date_time_str_ground_truth.LANGUAGE='03192021_1407';
date_time_str_ground_truth.MOTOR='03192021_1415';
date_time_str_ground_truth.RELATIONAL='03222021_1201';
date_time_str_ground_truth.SOCIAL='03222021_1212';
date_time_str_ground_truth.WM='03222021_1219';
% Activation
warning("Altered for activation - see comments SMN: ALTERED")
date_time_str_ground_truth.SOCIAL='04122022_2317'; %SMN: ALTERED

% Task versus task (supplementary) - 05/03/22
date_time_str_ground_truth.EMOTION_v_GAMBLING='05032022_1425';


%%%%%%%% Single task summary filename info %%%%%%%%

switch data_origin
case 'farnam'
    
switch summary_type
case {'visualize_tpr','visualize_gt','dcoeff'}
    
    switch grsize
        case 40
            date_time_str_combined='12142021';
        case 80
            date_time_str_combined='12142021';
        case 120
            date_time_str_combined='12152021';
        otherwise
            error('Group size not defined');
    end
    
    
    
otherwise
    switch stat_type
        case 'Parametric_Bonferroni'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='05052021_2345';
                date_time_str_results.GAMBLING_v_REST='05092021_2248';
                date_time_str_results.LANGUAGE_v_REST='05092021_2253';
                date_time_str_results.MOTOR_v_REST='05092021_2301';
                date_time_str_results.RELATIONAL_v_REST='05092021_2306';
                date_time_str_results.SOCIAL_v_REST='05092021_2309';
                date_time_str_results.WM_v_REST='05092021_2348';
                date_time_str_results.REST_v_REST2='05062021_2130';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='05062021_1540';
                date_time_str_results.GAMBLING_v_REST='05082021_0915';
                date_time_str_results.LANGUAGE_v_REST='05082021_0921';
                date_time_str_results.MOTOR_v_REST='05082021_0925';
                date_time_str_results.RELATIONAL_v_REST='05082021_1030';
                date_time_str_results.SOCIAL_v_REST='05082021_1100';
                date_time_str_results.WM_v_REST='05082021_1123';
                date_time_str_results.REST_v_REST2='05062021_2157';
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='05062021_0041';
                date_time_str_results.GAMBLING_v_REST='05062021_1614';
                date_time_str_results.LANGUAGE_v_REST='05062021_2059';
                date_time_str_results.MOTOR_v_REST='05062021_2117';
                date_time_str_results.RELATIONAL_v_REST='05062021_2305';
                date_time_str_results.SOCIAL_v_REST='05062021_2331';
                date_time_str_results.WM_v_REST='05062021_2333';
                date_time_str_results.REST_v_REST2='05062021_2131';
            end
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
                date_time_str_results.EMOTION_v_GAMBLING='05022022_1822'; % 050322
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
        case 'FDR' % nonparametric FDR procedure - not using due to feasibility in completing req'd perms
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
                date_time_str_results.GAMBLING_v_REST='03022020_1739';
                date_time_str_results.LANGUAGE_v_REST='02182021_1752';
                date_time_str_results.MOTOR_v_REST='02292020_1853';
                date_time_str_results.RELATIONAL_v_REST='02242020_1715';
                date_time_str_results.SOCIAL_v_REST='03012020_1942';
                date_time_str_results.WM_v_REST='02252020_1931';
                date_time_str_results.REST_v_REST2='01232021_1500';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_0752';
                date_time_str_results.GAMBLING_v_REST='12182020_0617';
                date_time_str_results.LANGUAGE_v_REST='12202020_0453';
                date_time_str_results.MOTOR_v_REST='12202020_0348';
                date_time_str_results.RELATIONAL_v_REST='12202020_0356';
                date_time_str_results.SOCIAL_v_REST='12202020_0436';
                date_time_str_results.WM_v_REST='12202020_0436';
                date_time_str_results.REST_v_REST2='01292021_0726';
                date_time_str_results.EMOTION_v_GAMBLING='05022022_2059'; % 050322
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
                date_time_str_results.REST_v_REST2='01232021_1834';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12172020_1731';
                date_time_str_results.GAMBLING_v_REST='12182020_1501';
                date_time_str_results.LANGUAGE_v_REST='12202020_1501';
                date_time_str_results.MOTOR_v_REST='12202020_1215';
                date_time_str_results.RELATIONAL_v_REST='12202020_1217';
                date_time_str_results.SOCIAL_v_REST='12202020_1308';
                date_time_str_results.WM_v_REST='12202020_1259';
                date_time_str_results.REST_v_REST2='01292021_0739';
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
         case 'Constrained_FWER'
            if grsize==40
                date_time_str_results.EMOTION_v_REST='05032021_2121';
                date_time_str_results.GAMBLING_v_REST='05042021_0044';
                date_time_str_results.LANGUAGE_v_REST='05042021_0141';
                date_time_str_results.MOTOR_v_REST='05042021_0129';
                date_time_str_results.RELATIONAL_v_REST='05042021_0217';
                date_time_str_results.SOCIAL_v_REST='05042021_0224';
                date_time_str_results.WM_v_REST='05042021_0251';
                date_time_str_results.REST_v_REST2='05042021_1805';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='05042021_1814';
                date_time_str_results.GAMBLING_v_REST='05042021_1801';
                date_time_str_results.LANGUAGE_v_REST='05042021_1740';
                date_time_str_results.MOTOR_v_REST='05042021_1916';
                date_time_str_results.RELATIONAL_v_REST='05042021_1915';
                date_time_str_results.SOCIAL_v_REST='05042021_2007';
                date_time_str_results.WM_v_REST='05042021_2106';
                date_time_str_results.REST_v_REST2='05042021_1929';
                % Activation
                    date_time_str_results.SOCIAL='04122022_2317'; %SMN: ALTERED
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='05042021_2334';
                date_time_str_results.GAMBLING_v_REST='05052021_0013';
                date_time_str_results.LANGUAGE_v_REST='05052021_0252';
                date_time_str_results.MOTOR_v_REST='05052021_0233';
                date_time_str_results.RELATIONAL_v_REST='05052021_0634';
                date_time_str_results.SOCIAL_v_REST='05052021_0615';
                date_time_str_results.WM_v_REST='05052021_0702';
                date_time_str_results.REST_v_REST2='05052021_0950';
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
                date_time_str_results.REST_v_REST2='01232021_1559';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12252020_0030';
                date_time_str_results.GAMBLING_v_REST='12182020_2202';
            	date_time_str_results.LANGUAGE_v_REST='12202020_2342';
                date_time_str_results.MOTOR_v_REST='12202020_1848';
                date_time_str_results.RELATIONAL_v_REST='12202020_1845';
                date_time_str_results.SOCIAL_v_REST='12202020_1940';
                date_time_str_results.WM_v_REST='12202020_1928';
                date_time_str_results.REST_v_REST2='01292021_0609';
                date_time_str_results.EMOTION_v_GAMBLING='05032022_0110'; % 050322
                % Activation
                    date_time_str_results.SOCIAL='04122022_2100'; %SMN: ALTERED
            elseif grsize==120
                date_time_str_results.EMOTION_v_REST='01262021_0053';
                date_time_str_results.GAMBLING_v_REST='01272021_2332';
                date_time_str_results.LANGUAGE_v_REST='01272021_2357';
                date_time_str_results.MOTOR_v_REST='01282021_0023';
                date_time_str_results.RELATIONAL_v_REST='01282021_0055';
                date_time_str_results.SOCIAL_v_REST='01272021_2244';
                date_time_str_results.WM_v_REST='01272021_2310';
                date_time_str_results.REST_v_REST2='01282021_2135';
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
                date_time_str_results.REST_v_REST2='01232021_1506';
            elseif grsize==80
                date_time_str_results.EMOTION_v_REST='12242020_2038';
                date_time_str_results.GAMBLING_v_REST='12242020_2301';
                date_time_str_results.LANGUAGE_v_REST='12252020_0030';
                date_time_str_results.MOTOR_v_REST='12242020_2323';
                date_time_str_results.RELATIONAL_v_REST='12242020_2354';
                date_time_str_results.SOCIAL_v_REST='12242020_2347';
                date_time_str_results.WM_v_REST='12252020_0007';
                date_time_str_results.REST_v_REST2='01292021_0653';
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





%%%%%%%% Set io filenames %%%%%%%%

if ~contains(summary_type,'visualize') % visualizations rely on combined summary, not these individual task files

    % ground truth variables
    ground_truth_results_basename_prefix=['ground_truth__',task,'_',stat_type_gt,'_',date_time_str_ground_truth.(task)];
    ground_truth_filename=[output_dir,ground_truth_results_basename_prefix,'.mat'];
    ground_truth_dcoeff_filename=[output_dir,ground_truth_results_basename_prefix,'_dcoeff.mat'];

    if any(strcmp(summary_type,{'calculate_tpr','positives'})) || strcmp(summary_type,'summarize_fprs')
        
        % original (un-summarized) benchmarking variables
        benchmarking_results_basename_prefix=['results__',task,fpr_str,'_',stat_type,'_','grsize',num2str(grsize),'_',date_time_str_results.(task)];
        results_filename=[output_dir,benchmarking_results_basename_prefix,'.mat'];

        % summarized benchmarking variables (output 1)
        benchmarking_summary_filename=[output_dir,benchmarking_results_basename_prefix,'_summary.mat'];

        % create new timestamp for a new file
        date_time_str_combined=date_time_str_now;
    end

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
combined_by_task_filename_prefix=[combined_by_task_dir,combined_basename_prefix];

% log filenames
% name with combined prefix since using combined data
log_dir=[combined_by_task_dir,'logs/'];
log_filename_prefix=[log_dir,combined_basename_prefix];

% set Shen edge groups file (default for benchmarking) in case can't get from original data
edge_groups_filename=[output_dir,'edge_groups.mat'];


