function [any_significant,con_mat,pval,edge_stats__target]=NBSedge_level_parametric_corr(varargin)
%This script performs Bonferroni or FDR correction over all edges in the
%network, with uncorrected p-values of each edge determined parametrically
%
%   [...]=NBSstats(STATS,H,GLM) calls uses GLM as input to NBSglm to compute
%   statistics. Writes out progress to a listbox with handle H. If writing
%   progress to H fails, progress is written to screen instead.
%
%   A STATS structure contains the following fields:
%       STATS.thresh:     (Not used here.) Primary test statistic threshold
%       STATS.alpha:      Corrected significance (user specified), network 
%                         components not satisfying alpha signficance are 
%                         not reported
%       STATS.N:          Number of nodes in network
%       STATS.test_stat:  1 x J array of the observed test statistics. J is
%                         the number of edges. Test statistics are
%                         computed with NBSglm. Columns are mapped to
%                         edges such that column i=1:J corresponds to the
%                         edge with index ind_upper(i), where ind_upper are
%                         the indexes of the upper trianguler elements. 
%                         ind_upper = find(triu(ones(N,N),1)); 
%       STATS.size        (Not used here.) 'extent' | 'intensity' 
%                         Measure used to assess size of a network 
%                         component  
%                          
%   Outputs:
%       N_CNT:            (Not applicable here.) Number of network components 
%       satisfying alpha 
%                         significance
%       CON_MAT:          single cell containing 1 x J array of significant 
%       edges after correction 
%       PVAL:             1 x J array of corrected p-values 
%   
%   Remarks:
%       None
%
%   Stephanie Noble (@sneuroble on github)
%   adapted from Andrew Zalesky's NBS scripts

STATS=varargin{1}; 
if nargin==2
    %Handle to listbox
    H=varargin{2}; 
elseif nargin==3
    %Handle to GLM
    H=varargin{2};
    GLM=varargin{3}; 
end

%Is BGL available?
if exist('components','file')==2
    %Use components.m provided by MatlabBGL, otherwise use get_components.m
    bgl=1;
else bgl=0;
end

% Specify no permutations for GLM
GLM.perms=1;

%Index of upper triangular elements of connectivity matrix
N=STATS.N; % n nodes
J=N*(N-1)/2; % n edges
ind_upper=find(triu(ones(N,N),1)); 

% Uncorrected p-values
GLM=NBSglm_setup_smn(GLM);
edge_stats__target=NBSglm_smn(GLM);

if strcmp(GLM.test,'onesample') 
    error('Under development.');
elseif strcmp(GLM.test,'ttest') 
    %if strcmp(ttest_type,'paired') % TODO: create
    warning('Assuming paired sample and right-tailed t-test.');
    df = GLM.n_observations/2-1; % observations are 2 per subject
    %elseif strcmp(ttest_type,'unpaired')
    %    df=GLM.n_observations-2; % assuming equal variances
    %end
    p_uncorr = tcdf(-edge_stats__target,df);
elseif strcmp(GLM.test,'ftest') 
    error('Under development.');
    df1=GLM.n_predictors-1; % TODO: check this
    df2=GLM.n_observations-GLM.n_predictors; % TODO: check this
    df1_vec=repmat(df1,size(edge_stats__target));
    df2_vec=repmat(df2,size(edge_stats__target));
    p_uncorr = fcdf(-edge_stats__target,df1_vec,df2_vec);
else
    error('Invalid or no GLM test type specified. Aborting.')
    % NOTE: one-sample statistic currently implemented in NBS is only the mean, not standardized (ie for a t-statistic). 
end

% Corrected p-values
switch STATS.statistic_type
    case 'Parametric_Bonferroni'
        pval = p_uncorr*length(edge_stats__target);
    case 'Parametric_FDR'
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC170937/
        % q value is a function of FDR that guarantees that q increases with p
        [~, pval] = mafdr(p_uncorr);
    otherwise
        error('Invalid or no correction specified. Aborting.');
end

con_mat{1}=pval(:)<STATS.alpha;
any_significant=any(con_mat{1});



