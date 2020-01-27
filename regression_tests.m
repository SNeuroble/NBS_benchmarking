%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NBS-based method regression tests
% Simple test of NBS calculation of cluster statistic under different 
% conditions compared with gold standards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

config_file='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/cNBS/config_files/config_nbs.m';
run(config_file);
benchmark_setup;

%% Run NBS & compare with ground truth OR Create a new ground truth

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end

if create_new_standard
    
    % First ground truth was built on 11/20/2019
    % Algorithms were verified manually on 11/12/2019
    
    fprintf('NOT running regression test - creating new standard instead.\n')
    
    nbs_gt=NBSrun_smn(UI);
    mkdir(ground_truth_dir);
    gt_filename=[ground_truth_dir,'nbs_ground_truth__',UI.statistic_type.ui,size_str,'_',datestr(now,'mmddyyyy')];
    save(gt_filename,'nbs_gt','UI_gt');
    return;

else % calculate new and compare
    
    % load 'nbs_gt' which has the UI above as a field
    % thus, can replicate the results with: nbs_gt=NBSrun_smn(nbs_gt.UI);

    gt_filename=[ground_truth_dir,'nbs_ground_truth__',statistic_type,size_str,'_',ground_truth_date];
    load(gt_filename); % options: 'nbs_ground_truth__TFCE' | 'nbs_ground_truth__size_extent' | 'nbs_ground_truth__constrained' 
    nbs_new=NBSrun_smn(UI);
    % TODO: 
    
end

%% Compare results

% edge and cluster stats should be identical
ok.edge_stats=isequal(nbs_new.NBS.edge_stats,nbs_gt.NBS.edge_stats);
ok.cluster_stats=isequal(nbs_new.NBS.cluster_stats,nbs_gt.NBS.cluster_stats);

% pvals based on nonparametric test so will differ each time
% experimentally, this ranges from 0.2-0.85 (sometimes getting NAN)
ok.pvals_corr=corr(nbs_new.NBS.pval(:),nbs_gt.NBS.pval(:));

disp(ok)

