%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs NBS from command line
% Directions:
% 1. fill out setparams_cl
% 2. run this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup
setparams_cl;
addpath(genpath(nbs_dir));
addpath(genpath(nbs_addon_dir));

% Developer parameter changes
if testing 
    n_perms=test_n_perms;
end

% Move parameters into UI variable (comes from the NBS GUI)
UI.method.ui=nbs_method;
UI.design.ui=design_matrix_file;
UI.contrast.ui=contrast;
UI.test.ui=nbs_test_stat;
UI.perms.ui=n_perms;
UI.thresh.ui=tthresh_first_level;
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type;
UI.size.ui=cluster_size_type;
UI.edge_groups.ui=edge_groups_file;
UI.exchange.ui=exchange;
UI.matrices.ui=data_file;

% run command line NBS (TODO: maybe rename NBSrun_cl)
nbs=NBSrun_smn(UI);

% results for cNBS
if ischar(pthresh_second_level); pthresh_second_level=str2num(pthresh_second_level); end
sig_results=nbs.NBS.pval<pthresh_second_level;

