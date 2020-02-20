% tiny script

%% Set paths, params, packages

% set paths and params
setpaths;
setparams;
addpath(genpath(nbs_dir));
addpath(genpath(other_scripts_dir));
subIDs_suffix='_subIDs.txt'; % TODO: move elsewhere

% octave stuff
page_screen_output (0); % for octave - also switched fprintf->printf
% pkg load struct
% pkg load parallel
% pkg load ndpar

%% Make data
n_subs=40;
n_nodes=268;
m=rand(n_nodes,n_nodes,n_subs*2);
t=triu(ones(size(m(:,:,1))),1);
t2=logical(repmat(t,1,1,size(m,3)));
m=m(t2);
m=reshape(m,n_nodes*(n_nodes-1)/2,n_subs*2);

%% Make design matrix and contrast

dmat=zeros(n_subs_subset*2,n_subs_subset+1);
dmat(1:(n_subs_subset),1)=1;
dmat((n_subs_subset+1):end,1)=-1;
for i=1:n_subs_subset
    dmat(i,i+1)=1;
    dmat(n_subs_subset+i,i+1)=1;
end

% set up contrasts - positive and negative
nbs_contrast=zeros(1,n_subs_subset+1);
nbs_contrast(1)=1;

nbs_contrast_neg=nbs_contrast;
nbs_contrast_neg(1)=-1;

% set up exchange
nbs_exchange=[1:n_subs_subset, 1:n_subs_subset];


%% Assign params to structures

% assign repetition parameters to rep_params
rep_params.data_dir=data_dir;
rep_params.testing=testing;
rep_params.mapping_category=mapping_category;
rep_params.n_repetitions=n_repetitions;
rep_params.n_subs_subset=n_subs_subset;
rep_params.do_TPR=do_TPR;
if do_TPR; rep_params.task_condition=task_condition;
end
rep_params.non_task_condition=non_task_condition;

% assign NBS parameters to UI (see NBS.m)
UI.method.ui=nbs_method; % TODO: revise to include vanilla FDR
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms; % previously: '5000'
UI.thresh.ui=tthresh_first_level; % p=0.01
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; % 'Size' | 'TFCE' | 'Constrained' | 'SEA'
UI.size.ui=cluster_size_type; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
if strcmp(cluster_stat_type,'cNBS') || strcmp(cluster_stat_type,'SEA')
    UI.edge_groups.ui=edge_groups; % smn
end
UI.exchange.ui=nbs_exchange;

%% Run
subsample_and_NBS(1,m,UI,n_subs,rep_params.n_subs_subset,do_TPR,nbs_contrast_neg);
