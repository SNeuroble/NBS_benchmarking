% Set input, output, and script directories depending on computer

switch getenv('USER');

case 'smn33'
    
    % data directory - input data is assumed to be n_nodes x n_nodes x n_subjects typically HCP toy data)
    data_dir='/data15/mri_group/smn33_data/hcp_1200/matrices/';
    %data_path='/mnt/store1/mridata2/mri_group/smn33_data/hcp/data_01ffd_v7_3.mat';
    
    % results directory
    output_dir='/data15/mri_group/smn33_data/NBS_benchmarking_results/';

    % NBS toolbox
    nbs_dir='/mridata2/home2/smn33/scripts/NBS1.2';

    % misc scripts used for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas
    other_scripts_dir='/mridata2/home2/smn33/scripts/matlab/myscripts/general_mri_new/general_mri'; 

case 'steph'
    
    % data directory (mount!) - input data is assumed to be n_nodes x n_nodes x n_subjects
    data_dir='/mnt/d3_smn33/hcp_1200/matrices/'; % mounted
    %data_path='/mnt/smn33_data/hcp/data_01ffd_v7_3.mat'; % mounted
    %data_path='/Users/steph/Steph-Lab/Misc/ConstableLab/MRRC Neuroinformatics/resources/scripts/cpm_hackathon/test_data/HCP900_rest_n50.mat'; % HCP toy data
    
    % results directory
    output_dir='/Users/steph/Steph-Lab/NBS_benchmarking/results_benchmarking/';
    
    % NBS toolbox
    nbs_dir='/Users/steph/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
    
    % misc scripts used for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas
    other_scripts_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/general_mri';


case 'ubuntu'

    % data directory - input data is assumed to be n_nodes x n_nodes x n_subjects typically HCP toy data)
    %data_dir='/data15/mri_group/smn33_data/hcp_1200/matrices/';
    data_dir='/home/ubuntu/data/matrices/';
    %data_path='/mnt/store1/mridata2/mri_group/smn33_data/hcp/data_01ffd_v7_3.mat';
    
    % results directory
    output_dir='/home/ubuntu/data/NBS_benchmarking_results/';

    % NBS toolbox
    nbs_dir='/home/ubuntu/scripts/NBS1.2/';

    % misc scripts used for cNBS and for summarization - structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas
    other_scripts_dir='/home/ubuntu/scripts/NBS_benchmarking/atlas_scripts/'; 



end


