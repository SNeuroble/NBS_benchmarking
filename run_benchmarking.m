% Do NBS-based method benchmarking (cNBS, TFCE, etc)
%
% main outputs:
% edge_stats_all: mean and sd of edge_stats_all
% cluster_stats_all: mean and sd of cluster_stats_all
% pvals_all: total # positives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function run_benchmarking

%% Setup
% make sure config files in NBS_benchmarking are correct

[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming NBS_benchmarking is current folder
addpath(genpath(current_path));
%global m; % make global var so parallel jobs don't need to copy
setup_benchmarking;
pkg load ndpar % must have installed ndpar

%% Setup for repetitions

%resume_from_previous=0;
%if exist('this_repetition') % resuming from another process
%    response=input(sprintf('%0.0d perms already completed. Resume? ',this_repetition),'s');
%    if strcmp(response,'y')
%        reps_completed_previously=this_repetition;
%        resume_from_previous=1;
%    end
%end
%if resume_from_previous==0 % starting from the beginning
%    reps_completed_previously=0;
%    any_sig_result=0;
%    any_sig_result_neg=0;
%    edge_stats_all=zeros(n_nodes*(n_nodes-1)/2,rep_params.n_repetitions);
%    edge_stats_all_neg=zeros(n_nodes*(n_nodes-1)/2,rep_params.n_repetitions);
%    if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
%        cluster_stats_all=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
%        cluster_stats_all_neg=zeros(length(unique(edge_groups))-1,1,rep_params.n_repetitions); % minus 1 to not count "zero"
%        pvals_all=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
%        pvals_all_neg=zeros(length(unique(UI.edge_groups.ui))-1,rep_params.n_repetitions); % minus 1 to not count "zero"
%    else
%        cluster_stats_all=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
%        cluster_stats_all_neg=zeros(n_nodes,n_nodes,rep_params.n_repetitions); 
%        pvals_all=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
%        pvals_all_neg=zeros(n_nodes*n_nodes,rep_params.n_repetitions);
%    end
%end

% randomly subsample subject IDs into groups 
for i=1:rep_params.n_repetitions
	ids=randperm(n_subs,n_subs_subset)';
	if do_TPR
		ids=[ids;ids+n_subs];
	end
	ids_sampled(:,i)=ids;
end

if n_workers>=nproc
	fprintf('Specified %d workers but only %d available - and don''t want to max out CPUs. Setting to max available minus 1.\n',n_workers,nproc);
	n_workers=nproc-1;
end

%% Do subsampling and NBS



if rep_params.testing; fprintf('*** TESTING MODE ***\n'); end
fprintf('Starting benchmarking repetitions across %d workers.\n',n_workers);

if testing
	response=input('About to begin parallel. Continue?'); % TODO: remove
end

[edge_stats_all,edge_stats_all_neg,cluster_stats_all,cluster_stats_all_neg,...
pvals_all,pvals_all_neg,any_sig_result,any_sig_result_neg] ...
= ndpar_arrayfun(n_workers,@(this_repetition) subsample_and_NBS(...
this_repetition,m(:,:,ids_sampled(:,this_repetition)),UI,n_subs,...
rep_params.n_subs_subset,do_TPR,nbs_contrast_neg),...
[1:rep_params.n_repetitions],"CatDimensions",[1,1,2,2,2,2,1,1]);
% this_repetition,m,UI,n_subs,rep_params.n_subs_subset,do_TPR,nbs_contrast_neg),...
%[1:rep_params.n_repetitions],"CatDimensions",[1,1,2,2,2,2,1,1]);
 
if testing
	response=input('About to begin parallel. Continue?'); % TODO: remove
end

if strcmp(UI.statistic_type.ui,'Constrained') || strcmp(UI.statistic_type.ui,'SEA')
   cluster_stats_all=squeeze(cluster_stats_all);
   cluster_stats_all_neg=squeeze(cluster_stats_all_neg);
end


%% Save results

mkdir(output_dir)

if strcmp(UI.statistic_type.ui,'Size'); size_str=['_',UI.size.ui];
else; size_str='';
end
if testing; test_str='_testing'; else test_str=''; end
if do_TPR; condition_str=['_',rep_params.task_condition]; else condition_str=['_',rep_params.non_task_condition]; end

output_filename=[output_dir,'nbs_benchmark_results__',condition_str,'_',UI.statistic_type.ui,size_str,test_str,'_',datestr(now,'mmddyyyy_HHMM'),'.mat'];

fprintf('Saving results in %s\n',output_filename)
save(output_filename,'edge_stats_all','cluster_stats_all','pvals_all','any_sig_result','edge_stats_all_neg','cluster_stats_all_neg','pvals_all_neg','any_sig_result_neg','UI','rep_params');

% show that results are available in the workspace
previous_results_filename__already_loaded=output_filename;

end


