%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Setup for running NBS benchmarking
% This will load data and set up the parameters needed to run NBS benchmarking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% set paths and params
setpaths;
setparams;
addpath(genpath(nbs_dir));
addpath(genpath(other_scripts_dir));

if do_TPR

    % super hacky stuff

    data_dir='/data15/mri_group/smn33_data/hcp_1200/matrices/';
    data1_name='LANGUAGE_RL';
    data2_name='REST_RL';
    subnames1_file=[data_dir,data1_name,'_subnames_short.txt'];
    subnames2_file=[data_dir,data2_name,'_subnames_short.txt'];

    % compare subnames (thanks https://www.mathworks.com/matlabcentral/answers/358722-how-to-compare-words-from-two-text-files-and-get-output-as-number-of-matching-words)
    % Read the files into strings:
    subnames1_str = fileread(subnames1_file);
    subnames2_str = fileread(subnames2_file);
    fprintf([subnames1_file subnames2_file]);
    % Split the strings to words:
    subnames1 = strsplit(subnames1_str, newline);
    subnames2 = strsplit(subnames2_str, newline);
    % Get the common words:
    [subnames_common, iA, iB] = intersect(subnames1, subnames2);
    subnames_common=subnames_common(2:end); % bc the first find is empty - TODO add a check here first
    n_subs=length(subnames_common); % TODO: redefined below

    % read data
    load_data='y';
    if exist('m','var');
        load_data=input('Some data is already loaded in workspace. Replace? (y/n)','s');
    end
    if strcmp(load_data,'y')
        
        fprintf('Loading data: %s\n',data_path);
        fprintf('Reading subject ');
        template_file=[data_dir,data1_name,'/',subnames_common{1},'_',data1_name,'_GSR_matrix.txt'];
        fprintf(template_file)
        template=importdata(template_file);
        m=zeros(size(template,1),size(template,2),n_subs*2);

        for i = 1:length(subnames_common)
            fprintf([subnames_common{i},' ']);
            %'d3_smn33/hcp_1200/archives/home/ec2-user/data/matrices/LANGUAGE_RL/100206_LANGUAGE_RL_GSR_matrix.txt'
            this_file_gr1 = [data_dir,data1_name,'/',subnames_common{i},'_',data1_name,'_GSR_matrix.txt'];
            m(:,:,i) = importdata(this_file_gr1);
            this_file_gr2 = [data_dir,data2_name,'/',subnames_common{i},'_',data2_name,'_GSR_matrix.txt'];
            m(:,:,n_subs*2+i) = importdata(this_file_gr2);
        end
    elseif strcmp(load_data,'n')
        fprintf('Using previously loaded data and assuming already reordered.\n');
    else; error('Input must be y or n.');
    end

else
    % load data
    load_data='y';
    if exist('m','var');
        load_data=input('Some data is already loaded in workspace. Replace? (y/n)','s');
    end    

    if strcmp(load_data,'y')
        sprintf('Loading data: %s\n',data_path);
        variableInfo = who('-file',data_path);
        if ismember('data', variableInfo) % returns true
            m=struct2array(load(data_path,'data'));
        else
            error('Could not find variable ''data'' in data file.');
        end
        m=reorder_matrix_by_atlas(m,mapping_category);
    elseif strcmp(load_data,'n')
        fprintf('Using previously loaded data and assuming already reordered.\n');
    else; error('Input must be y or n.');
    end

    % get n subs from data
    n_subs=size(m,3);

    % assume data is rest, although it could technically also be within task
    data1_name='REST';

end

% get node size from data 
n_nodes=size(m,1); % assuming square


% make design matrix for 2-sample t-test (equal group size)
if do_TPR
    % set up design matrix for one-sample t-test
    % data should be organized: s1_gr1,s2_gr1, ... , sn-1_group2, sn_group2
    dmat=zeros(n_subs_subset*2,n_subs_subset+1);
    dmat(1:(n_subs_subset),1)=1;
    dmat((n_subs_subset+1):end,1)=-1;
    for i=1:n_subs_subset
        dmat(i,i+1)=1;
        dmat(n_subs_subset+i,i+1)=1;
    end

    nbs_contrast=zeros(1,n_subs_subset+1);
    nbs_contrast(1)=1;

    nbs_exchange=[1:n_subs_subset, 1:n_subs_subset];
else
    % set up design matrix for two-sample t-test
    % data should be organized: s1_gr1, ... sn_gr1, sn+1_gr2, ... s2*n_gr2
    dmat=zeros(n_subs_subset,2);
    dmat(1:(n_subs_subset/2),1)=1;
    dmat((n_subs_subset/2+1):end,2)=1;
    
    nbs_contrast=[1,-1];

    nbs_exchange='';
end

% make edge groupings (for cNBS)
edge_groups=load_atlas_edge_groups(n_nodes,mapping_category);
edge_groups=tril(edge_groups,-1);
% TODO: in NBS function, should we require zero diag? Automatically clear diag? Something else?

% assign params to structures
% should be able to run config file, load rep_params and UI from reference, and replicate reference results

% assign repetition parameters to rep_params
rep_params.data_path=data_path;
rep_params.testing=testing;
rep_params.do_simulated_effect=do_simulated_effect;
rep_params.networks_with_effects=networks_with_effects;
rep_params.mapping_category=mapping_category;
rep_params.n_repetitions=n_repetitions;
rep_params.n_subs_subset=n_subs_subset;
rep_params.do_TPR=do_TPR;
rep_params.task=data1_name;

% assign NBS parameters to UI (see NBS.m)
UI.method.ui=nbs_method; % TODO: revise to include vanilla FDR
UI.design.ui=dmat;
UI.contrast.ui=nbs_contrast;
UI.test.ui=nbs_test_stat; % alternatives are one-sample and F-test
UI.perms.ui=n_perms; % previously: '5000'
UI.thresh.ui=zthresh_first_level; % p=0.01
UI.alpha.ui=pthresh_second_level;
UI.statistic_type.ui=cluster_stat_type; % 'Size' | 'TFCE' | 'Constrained' | 'SEA'
UI.size.ui=cluster_size_type; % 'Intensity' | 'Extent' - only relevant if stat type is 'Size'
UI.edge_groups.ui=edge_groups; % smn
UI.exchange.ui=nbs_exchange;



