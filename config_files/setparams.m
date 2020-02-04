%% User-defined parameters for running NBS benchmarking

% Resampling parameters
testing=0; % developers only - speeds up analyses for troubleshooting but inappropriate for inference
n_workers=4; % num parallel workers for parfor, best to use # workers = # cores
do_simulated_effect=0;
networks_with_effects=[1,5]; % networks to add simulated effects into - only relevant if adding effect
mapping_category='subnetwork'; % for cNBS
n_repetitions=1000;
n_subs_subset=40; % size of subset is full group size (N=n*2 for two sample t-test or N=n for one-sample)

% NBS parameters
% TODO: right now dmat/contrast only designed for t-test, and edge_groups can only be Shen atlas
nbs_method='Run NBS'; % TODO: revise to include vanilla FDR
nbs_contrast='[1,-1]'; % Do not change - right now this is the only option
nbs_test_stat='t-test'; % alternatives are one-sample and F-test
n_perms='2000'; % previously: '5000'
zthresh_first_level='3.1'; % p=0.01
pthresh_second_level='0.05';
cluster_stat_type='Size'; % 'Size' | 'TFCE' | 'Constrained' | 'SEA' - smn
cluster_size_type='Intensity'; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
nbs_exchange='';

% DEVELOPERS ONLY - Use a small subset of reps and perms to speed up troubleshooting - inappropriate for inference
if testing 
    n_repetitions=70;   
    n_perms='20';
end
