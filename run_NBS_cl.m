%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs NBS from command line
% Directions:
% 1. fill out setparams
% 2. run this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup
[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
setparams;
addpath(genpath(nbs_dir));

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

% Run command line NBS (TODO: maybe rename NBSrun_cl)
nbs=NBSrun_smn(UI);

% Significant results
if ischar(pthresh_second_level); pthresh_second_level=str2num(pthresh_second_level); end
sig_results=nbs.NBS.pval<pthresh_second_level;

% Map significant cNBS results to the edge level
% Note: Requires edge_groups_file to be a numerical matrix workspace 
% variable. cNBS results will be 1 x n subnetworks, where the first entry 
% corresponds to edges in edge group 1, second to edges in edge group 2, etc

if strcmp(cluster_stat_type,'Constrained')
    if ischar(edge_groups_file)
    else
        sig_edge_results=edge_groups;
        it=1;
        for i=unique(nonzeros(edge_groups))'
            sig_edge_results(edge_groups==i)=sig_results(it);
            it=it+1;
        end
    end
end

