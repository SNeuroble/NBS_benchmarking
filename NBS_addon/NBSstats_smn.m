function [any_significant,con_mat,pval,edge_stats__target,cluster_stats__target]=NBSstats_smn(varargin)
%NBSstats Computes network components among edges that survive a primary 
%test statistic threshold. Assigns a corrected p-value to each component 
%using permuted data that has been supplied as part of the structure STATS. 
%
%   [N_CNT,CON_MAT,PVAL]=NBSstats(STATS) operates on the network and
%   associated statistical data in the structure STATS. 
%
%   [...]=NBSstats(STATS,H) writes out progress to a listbox with handle H. 
%   If writing progress to H fails, progress is written to screen instead. 
%
%   [...]=NBSstats(STATS,H,GLM) GLM is mandatory if the permuted data has 
%   not been precomputed, which is the case when STATS.test_stat is empty. 
%   In this situation, NBSglm is repeatedly called to compute a test
%   statistic for each permutation. Slower than precomputation, but saves 
%   memory. 
%
%   A STATS structure contains the following fields:
%       STATS.thresh:     Primary test statistic threshold    
%       STATS.alpha:      Corrected significance (user specified), network 
%                         components not satisfying alpha signficance are 
%                         not reported
%       STATS.N:          Number of nodes in network
%       STATS.test_stat:  K+1 x J array of test statistics. The first 
%                         row is the oberved test statistic. The remaining 
%                         rows are samples of the test statistic under the 
%                         null hypothesis. Each column corresponds to a 
%                         seperate edge. K is number of permutations. J is
%                         the number of edges. The test statistics can be
%                         computed with NBSglm. Columns are mapped to
%                         edges such that column i=1:J corresponds to the
%                         edge with index ind_upper(i), where ind_upper are
%                         the indexes of the upper trianguler elements. 
%                         ind_upper = find(triu(ones(N,N),1)); 
%       STATS.size        'extent' | 'intensity' 
%                         Measure used to assess size of a network 
%                         component  
%       STATS.[other]     A number of fields have been added to support 
%                         network- and whole brain-level inference
%                          
%   Outputs:
%       N_CNT:            Number of network components satisfying alpha 
%                         significance
%       CON_MAT:          1 x N_CNT cell array of adjacency matrices. Each
%                         cell holds a N x N upper-triangular adjacency 
%                         matrix specifying a network component satisfying 
%                         alpha signifcance
%       PVAL:             1 x N_CNT array of corrected p-values 
%   
%   Remarks:
%       If no network components satisfy alpha significance, CON_MAT and
%       PVAL are returned empty and N_CNT = 0.  
%
%       STATS.test_stat is empty if the number of permutations is too large
%       to precompute. See Limit parameter in NBSrun for details. In this 
%       situation, NBSglm is repeatedly called (J times) to compute test 
%       statistics for each permutation. 
%
%   azalesky@unimelb.edu.au

%% Setup

%Number of most recent permutations to display in listbox
nDisp=5;
new_str={};

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

% Set statistic type to numeric - size (extent/intensity) (0), constrained (1), or TFCE (2)
Intensity=0;
switch STATS.statistic_type
    case 'Size'
        STATS.statistic_type_numeric=1;
        % size of a component measured using extent or intensity?
        Intensity=strcmp(STATS.size,'Intensity');
    case 'TFCE'
        STATS.statistic_type_numeric=2;
    case 'Constrained' % cNBS with FDR correction (Simes)
        STATS.statistic_type_numeric=3; 
    case 'Constrained_FWER' % cNBS with FWER correction (Bonferroni)
        STATS.statistic_type_numeric=4;
    case 'SEA' % SEA
        STATS.statistic_type_numeric=5;
    case 'Omnibus' % Omnibus
        
        STATS.statistic_type_numeric=6;
        switch STATS.omnibus_type
            case 'Threshold_Positive' % Threshold positive
                STATS.omnibus_type_numeric=1;
            case 'Threshold_Both_Dir' % Threshold positive and negative
                STATS.omnibus_type_numeric=2;
            case 'Average_Positive' % Average positive
                STATS.omnibus_type_numeric=3;
            case 'Average_Both_Dir' % Average absolute value of positive and negative
                STATS.omnibus_type_numeric=4;
            case 'Multidimensional_cNBS' % Multidimensional null - uses network-level test stats to calculate Euclidean distance
                STATS.omnibus_type_numeric=5;
            case 'Between_minus_within_cNBS' % Between-network (typically positive in ground truth) minus within-network (typically negative)
                STATS.omnibus_type_numeric=6;
                error([STATS.omnibus_type, ' under development']);
            case 'Multidimensional_all_edges' % Multidimensional null - uses edge-level test stats to calculate Euclidean distance
                STATS.omnibus_type_numeric=7;
            otherwise
                error(sprintf('Omnibus type %s was provided - not a valid omnibus type.',STATS.statistic_type))
        end
    otherwise
        error(sprintf('Statistic type %s was provided - not a valid statistic type. Only ''Size'', ''TFCE'', ''Constrained'', or ''Omnibus'' allowed.',STATS.statistic_type))
end

if ~isempty(STATS.test_stat)
    precomputed=1;
    %Number of permutations
    K=size(STATS.test_stat,1)-1;
else
    precomputed=0;
    %Get number of permutations and set GLM.perms to 1, since NBSglm will be called separately for each permutation
    K=GLM.perms;
    GLM.perms=1;
end

do_parametric_for_2nd_level=NaN;
if STATS.use_edge_groups
    
    % special multidimensional null for cNBS, SEA, Omnibus - cNBS
    n_nulls=length(STATS.edge_groups.unique);
    do_parametric_for_2nd_level=1; % parametric or nonparametric; TODO: move - should be arg from user
    
    % check mask
    [STATS.edge_groups,was_mask_flipped]=check_cNBS_mask(STATS.edge_groups);
else
    n_nulls=1;
end
null_dist=zeros(K,n_nulls);

%Index of upper triangular elements of connectivity matrix
N=STATS.N; % n nodes
J=N*(N-1)/2; % n edges
ind_upper=find(triu(ones(N,N),1)); 


%% Estimate cluster-level statistics - target and null

% 1. Unpermuted cluster statistics
y__target=GLM.y; % save for posterity
GLM=NBSglm_setup_smn(GLM);
edge_stats__target=get_univariate_stats(STATS,GLM,precomputed,1);
[cluster_stats__target,max_stat__target]=get_cluster_stats(edge_stats__target,STATS,ind_upper, N, Intensity,bgl);

% 2, Permuted cluster statistics
%First row of test_stat is the observed test statistics, so start at the second row
for i=2:K+1
    
    GLM.y=permute_signal(GLM);
    edge_stats=get_univariate_stats(STATS,GLM,precomputed,i);
    [~,max_stat]=get_cluster_stats(edge_stats, STATS, ind_upper, N, Intensity,bgl);
    null_dist(i-1,:)=max_stat;
    
    % only print every hundred
    if mod((i-1),100)==0
        str=sprintf('| %5d/%5d perms complete | latest element in null: %6.1f |',i-1,K,null_dist(i-1,end));
        % display no more than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow;
        catch;  fprintf([str,'\n']); end
    end
    
end


%% Perform correction
% This gets cluster-based statistics satisfying FWER threshold
[pval]=perform_correction(null_dist,cluster_stats__target,max_stat__target,do_parametric_for_2nd_level,STATS,K);
any_significant=any(pval(:)<STATS.alpha);
con_mat=0; % this is for original NBS visualization in GUI and needs to be in a specific format; instead, we'll use the stats directly

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [edge_groups,was_mask_flipped]=check_cNBS_mask(edge_groups)
    % TODO: might be more logical to have this in NBSrun_smn
    % Report zeros on the diagonal
    if any(diag(edge_groups.groups))
        warning('Constrained mask has nonzero entries on the diagonal. Including diagonal entries is not advised but will not be stopped.');
    end

    %  Force mask to be upper triangular only
    t=tril(edge_groups.groups,-1);
    was_mask_flipped=0;
    if any(t(:))
        edge_groups.groups=triu(edge_groups.groups');
        was_mask_flipped=1;
    end
end

function y_perm=permute_signal(GLM)
% Permute variable
%
% SMN: all of these said "use the same permutation for every GLM" but I
% don't think this is true... or I don't understand what they meant by this

if isempty(GLM.ind_nuisance) 
    %Permute signal 
    if exist('GLM.blk_ind','var')
        for j=1:GLM.n_blks
            y_perm(GLM.blk_ind(j,:),:)=...
            GLM.y(GLM.blk_ind(j,randperm(GLM.sz_blk)),:);
        end
    else                
        y_perm=GLM.y(randperm(GLM.n_observations)',:);
    end
else
    %Permute residuals
    % TODO: this returns altered GLM, not y_perm
    if exist('blk_ind','var')
        for j=1:GLM.n_blks
            GLM.resid_y(GLM.blk_ind(j,:),:)=...
            GLM.resid_y(GLM.blk_ind(j,randperm(GLM.sz_blk)),:);
        end
    else
        GLM.resid_y=GLM.resid_y(randperm(GLM.n_observations)',:);
    end

    %Add permuted residual back to nuisance signal, giving a realisation 
    %of the data under the null hypothesis (Freedma & Lane)
    
    if ~isempty(GLM.ind_nuisance)
        y_perm=GLM.resid_y+[GLM.X(:,GLM.ind_nuisance)]*GLM.b_nuisance;
    end

    %Flip signs
    %Don't permute first run % TODO: SMN - why is this here and not
    %elsewhere? also, shuffling already happened - whyis this
    %happening again here??
    if strcmp(GLM.test,'onesample')
        y_perm=y_perm.*repmat(sign(rand(GLM.n_observations,1)-0.5),1,GLM.n_GLMs);
    end

end
end


function [test_stat]=get_univariate_stats(STATS,GLM,precomputed,idx)
if precomputed
    %Precomputed test statistics
    test_stat=STATS.test_stat(idx,:); % TODO: here and next section - should not threshold if doing TFCE
else
    %Compute test statistics on the fly
%     test_stat=zeros(2,N*(N-1)/2); %I don't think this preallocation makes a difference since NBSglm is creating the variable then copying?
    test_stat=NBSglm_smn(GLM); % replaced glm so only does one perm, not the orig and then another
end
end

function [cluster_stats,null_stat]=get_cluster_stats(test_stat, STATS, ind_upper, N, Intensity,bgl)
% Returns:  cluster_stats_map: edge-wise map of cluster statistics
%           null_stat: permutation-specific statistic for building the null

% TODO: note that max is the null test stat for first two stat types but not cNBS - 
% should these functions instead just return the cluster vals but then do
% the max outside? Tough bc integrated in the calculation


switch STATS.statistic_type_numeric
    
    case 1 % do NBS +/- Intensity
        
        % Index of edges exceeding the primary threshold
        ind=ind_upper(test_stat>STATS.thresh);
        adj=spalloc(N,N,length(ind)*2);
        adj(ind)=1; adj=adj+adj';
        
        % note: second argument saved for null stat is max cluster size
        [cluster_stats,null_stat]=get_edge_components(adj,Intensity,test_stat,STATS.thresh,N,ind_upper,bgl);
       
    case 2 % do TFCE
        
        test_stat_mat=zeros(N,N);
        test_stat_mat(ind_upper)=test_stat;
        test_stat_mat=(test_stat_mat+test_stat_mat');
        
        cluster_stats=matlab_tfce_transform(test_stat_mat,'matrix');
        
        null_stat=max(cluster_stats(:));
        
    case {3,4} % do cNBS - returns a group-length vector
        
        % TODO: this should be done (here and below - see case 5) directly from the test_stat without making mat - just need network IDs 
        if STATS.use_preaveraged_constrained
            cluster_stats=test_stat;
        else
            test_stat_mat=zeros(N,N);
            test_stat_mat(ind_upper)=test_stat;
            test_stat_mat=(test_stat_mat+test_stat_mat');
            
            cluster_stats=get_constrained_stats(test_stat_mat,STATS.edge_groups);
        end

        null_stat=cluster_stats; % TODO - this is not a max value but instead 1 val for null per group - think of how to name
        
    case 5 % do Set Enrichment Analysis (SEA; GSEA minus the G bc not genetics)
        
        edge_groups_vec=STATS.edge_groups.groups(ind_upper); % TODO: do this during setup
        cluster_stats=get_enrichment_score(edge_groups_vec,STATS.edge_groups.unique,test_stat,1,0);
        
        null_stat=cluster_stats; % TODO - this is not a max value but instead 1 val for null per group - think of how to name
        
    case 6 % do Omnibus
        switch STATS.omnibus_type_numeric
            case 1 % Threshold positive
                cluster_stats=sum(+(test_stat>STATS.thresh));
                null_stat=cluster_stats;
            case 2 % Threshold positive and negative
                cluster_stats=sum(+(test_stat>STATS.thresh | test_stat<(-STATS.thresh)));
                null_stat=cluster_stats;
            case 3 % Average positive
                cluster_stats=sum(test_stat(test_stat>0));
                null_stat=cluster_stats;
            case 4 % Average absolute value of positive and negative
                cluster_stats=sum(abs(test_stat));
                null_stat=cluster_stats;
            case 5 % Multidimensional null
                % TODO: this should be done (here and above - see case 3) directly from the test_stat without making mat - just need network IDs 
                test_stat_mat=zeros(N,N);
                test_stat_mat(ind_upper)=test_stat;
                test_stat_mat=(test_stat_mat+test_stat_mat');

                if STATS.use_preaveraged_constrained
                    cluster_stats=test_stat_mat;
                else
                    cluster_stats=get_constrained_stats(test_stat_mat,STATS.edge_groups);
                end

                null_stat=cluster_stats; % TODO - this is not a max value but instead 1 val for null per group - think of how to name
            case 6 % Multidimensional null
                % DO NOT USE - UNDER DEVELOPMENT
                % TODO: this should be done (here and above - see case 3) directly from the test_stat without making mat - just need network IDs 
                test_stat_mat=zeros(N,N);
                test_stat_mat(ind_upper)=test_stat;
                test_stat_mat=(test_stat_mat+test_stat_mat');

                cluster_stats=get_constrained_stats(test_stat_mat,STATS.edge_groups);

                null_stat=cluster_stats; % TODO - this is not a max value but instead 1 val for null per group - think of how to name
        
            case 7 % Multidimensional null
                % Calculate Euclidean distance for each permutation from null centroid (assumed zeros)
                cluster_stats=sqrt(sum((test_stat).^2));
		null_stat=cluster_stats;
                
        end
        
 end

end


function [pval]=perform_correction(null_dist,target_stat,max_target_stat,do_parametric_for_2nd_level,STATS,K)
% Returns map of FWER-corrected pvals
%
% For multidimensional nulls (cNBS):
%   1st level: get within-group (uncorrected) pvals - nonparametric
%   2nd level: get cross-group FWER-corrected pvals - parametric or nonparametric

switch STATS.statistic_type_numeric
    
    case 1 % NBS
        
        [unique_stats,~,idx_unique_to_target_stat]=unique(target_stat); % added so didn't need to pass
        unique_pvals=arrayfun(@(this_stat) sum(+(this_stat<null_dist))/K, full(unique_stats)); % TODO: return here, prob expensive to convert sparse to full
        pval=unique_pvals(idx_unique_to_target_stat); % TODO: this puts stuff back into a vector but we want it in the original matrix dim to match cluster outputs
        
    case 2 % TFCE
        
        pval=arrayfun(@(this_stat) sum(+(this_stat<null_dist))/K, target_stat);
        
    case {3, 4, 5} % cNBS and SEA
        
        % estimate pvalue for each network
        pval_uncorr=zeros(size(STATS.edge_groups.unique));
        for i=1:length(STATS.edge_groups.unique)
            pval_uncorr(i)=sum(+(max_target_stat(i)<null_dist(:,i)))/K;
        end
        
        % maybe faster:
%         n_groups=length(STATS.edge_groups.unique);
%         pval_uncorr=arrayfun(@(this_group) sum(+(null_dist_multiple(,:)>=test_stat(this_group)))/K, n_groups); 
        
        if do_parametric_for_2nd_level
            if STATS.statistic_type_numeric==4 % Bonferroni
            %if STATS.do_Constrained_FWER_second_level % Bonferroni
                pval=pval_uncorr*length(STATS.edge_groups.unique); % FWER-corrected p-vals
                %sprintf('Trying new Bonferroni 2nd level\n');
            else % FDR
                
                % TODO: think more about procedure/output into pval variable
%                 [~,pval,~,~]=mafdr(pval_uncorr); % FDR-corrected p-vals
                % requires bioinformatics toolbox, so we're using the Simes procedure of FDR as
                % implemented in NBSfdr:
                
                %Simes procedure
                J=length(pval_uncorr);
                ind_srt=zeros(1,J); 
                [pval_uncorr_sorted,ind_srt]=sort(pval_uncorr);
                tmp=(1:J)/J*STATS.alpha;
                ind_sig=pval_uncorr_sorted<=tmp;
                
                pval=ones(1,J);
                pval(ind_srt(ind_sig))=0; % here, binary: 0 means significant (<alpha), 1 is not significant (>alpha) 
                
            end
            pval(pval>1)=1;
            
        else % bootstrapping
            error('Stopping - this section is under development and not ready for use');
            if STATS.statistic_type_numeric==4 % Bonferroni
            %if STATS.do_Constrained_FWER_second_level
                maxima_sorted=sort(max(null_dist_multiple),'descend');
                n_above_alpha=STATS.alpha*K;
                test_stat_threshold=maxima_sorted(n_above_alpha);
                % TODO: still need to calculate p-vals while accounting for non-interchangeability
            else
                ranked_null_dist_multiple=zeros(size(null_dist_multiple));
                for i=1:length(STATS.edge_groups.unique)
                    [~,p] = sort(null_dist_multiple(i,:),'descend');
                    r = 1:length(Data);
                    ranked_null_dist_multiple(i,p) = r;
                end
                % TODO: under development
            end
        end
        
    case 6 % Omnibus
        
        if STATS.omnibus_type_numeric==5 % Special for constrained
            % TODO: under development
            
            % Calculate Euclidean distance for each permutation from null centroid
            null_centroid=mean(null_dist); %  mean of each network across permutations (or assume 0?)
            null_dist__eucl_dist=sqrt(sum((null_dist-null_centroid).^2));
            
            % Calculate Euclidean distance for target stat from null centroid
            target_stat__eucl_dist=sqrt(sum((target_stat-null_centroid).^2));
            
            % Compare target distance to null distances
            pval=sum(+(target_stat__eucl_dist<null_dist__eucl_dist))/K;
        
        else
            pval=sum(+(target_stat<null_dist))/K;
        end
        
    end

end

