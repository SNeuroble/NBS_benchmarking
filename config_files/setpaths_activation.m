%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set input, output, and script directories depending on current workstation
% Will only run if system_dependent_paths=1 in setparams_bench.m
%
% directories
% - data_dir: data directory containing input data files each n_nodes x n_nodes x n_subjects
%   may be mounted, e.g., farnam results post-benchmarking onto local computer
% - output_dir: results directory
% - nbs_dir: NBS toolbox
% - other_scripts_dir: misc scripts used for summarization (structure_data,
% draw_atlas_boundaries, summarize_matrix_by_atlas)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_origin='farnam'; % 'farnam' or 'mrrc'

switch getenv('USER')
    
    case 'smn33' % farnam or mrrc
        
        switch data_origin
                
            case 'farnam'
                data_dir='/home/smn33/project/HCP_S1200_activation/by_nets/';
                output_dir='/home/smn33/project/benchmarking_results_activation/';
                nbs_dir='/home/smn33/scripts/NBS1.2';
                other_scripts_dir='/home/smn33/scripts/NBS_benchmarking/NBS_benchmarking/support_scripts/';
                
            case 'mrrc'
                data_dir='/data15/mri_group/smn33_data/hcp_1200/matrices/'; % symlink: '/mnt/dustin/data/S1200/matrices/S1200/matrices/'
                output_dir='/data15/mri_group/smn33_data/NBS_benchmarking_results/';
                nbs_dir='/mridata2/home2/smn33/scripts/NBS1.2';
                other_scripts_dir='/mridata2/home2/smn33/scripts/matlab/myscripts/general_mri_new/general_mri';
                
        end

    case 'steph' % input data dir is mounted as specified
        
        switch data_origin
            case 'farnam'
                data_dir='/Users/steph/Documents/data/mnt/project/HCP_S1200_activation/by_nets/'; 
                output_dir='/Users/steph/Documents/data/mnt/project/benchmarking_results_activation/';
            case 'mrrc'
                data_dir='/Users/steph/Documents/data/mnt/hcp_1200/matrices/';
                output_dir='/Users/steph/Documents/data/mnt/NBS_benchmarking_results/';
        end
        
        nbs_dir='/Users/steph/Steph-Lab/Misc/Software/scripts/Matlab/fmri/NBS1.2';
        other_scripts_dir='/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/NBS_benchmarking/support_scripts/';
        
    case 'ubuntu' % AWS
        
        data_dir='/home/ubuntu/data/matrices/';
        output_dir='/home/ubuntu/data/NBS_benchmarking_results/';
        nbs_dir='/home/ubuntu/scripts/NBS1.2/';
        other_scripts_dir='/home/ubuntu/scripts/NBS_benchmarking/support_scripts/';
        
    otherwise
        
        error('Specified to use system-dependent paths but paths undefined.');
        
end


