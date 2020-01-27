function [n_cnt,con_mat,pval]=NBSstats_smn(varargin)
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

%Number of most recent permutations to display in listbox
nDisp=5;

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
bgl=0;
if exist('components','file')==2
    %Use components.m provided by MatlabBGL, otherwise use get_components.m
    bgl=1;
end

% statistic type can be size (extent/intensity) (0), constrained (1), or TFCE (2)
switch STATS.statistic_type
    case 'Size'
        stat_type=0;
    case 'TFCE'
        stat_type=1;
    case 'Constrained'
        stat_type=2;
    otherwise
        error(sprintf('Statistic type %s was provided, but only ''Size'', ''TFCE'', or ''Constrained'' allowed.',STATS.statistic_type))
end


%Number of nodes
N=STATS.N; 

%Number of edges
J=N*(N-1)/2;

%Index of upper triangular elements of connectivity matrix
ind_upper=find(triu(ones(N,N),1)); 

if stat_type==0
    %Determine whether test statistics have been precomputed and determine
    %index of edges exceeding the primary threshold
    if ~isempty(STATS.test_stat)
        %Precomputed test statistics
        ind=ind_upper(STATS.test_stat(1,:)>STATS.thresh); % TODO: here and next section - should not threshold if doing TFCE
        %Number of permutations
        K=size(STATS.test_stat,1)-1;
    else
        %Compute test statistics on the fly
        test_stat=zeros(2,N*(N-1)/2);
        %Get desired number of permutations
        K=GLM.perms;
        %Set to 1, since NBSglm will be called separately for each permutation
        GLM.perms=1;
        test_stat=NBSglm(GLM);
        ind=ind_upper(test_stat(1,:)>STATS.thresh);
    end
    
    %Size of a component measured using extent or intensity?
    Intensity=0;
    if strcmp(STATS.size,'Intensity')
        %If size measure using intensity, create an N x N matrix cotaining the
        %test statistic for each edge minus the test statistic threshold
        %(primary threshold)
        Intensity=1;
        %Compute a test statistic matrix
        test_stat_mat=zeros(N,N);
        if ~isempty(STATS.test_stat)
            %Precomputed
            test_stat_mat(ind_upper)=STATS.test_stat(1,:)-STATS.thresh;
            test_stat_mat=(test_stat_mat+test_stat_mat');
        else
            %Not precomputed
            test_stat_mat(ind_upper)=test_stat(1,:)-STATS.thresh;
            test_stat_mat=(test_stat_mat+test_stat_mat');
        end
    end
    
    adj=spalloc(N,N,length(ind)*2);
    adj(ind)=1; adj=adj+adj';
    
else % get test stats without thresholding for cNBS and TFCE - smn % TODO - Come back here!
    ind=ind_upper;
    test_stat_mat=zeros(N,N);
    if ~isempty(STATS.test_stat)
        %Number of permutations
        K=size(STATS.test_stat,1)-1;
        test_stat_mat(ind_upper)=STATS.test_stat(1,:);
        test_stat_mat=(test_stat_mat+test_stat_mat');
    else
        %Compute test statistics on the fly
        test_stat=zeros(2,N*(N-1)/2);
        %Get desired number of permutations
        K=GLM.perms;
        %Set to 1, since NBSglm will be called separately for each permutation
        GLM.perms=1;
        test_stat=NBSglm(GLM);
        test_stat_mat(ind_upper)=test_stat(1,:);
        test_stat_mat=(test_stat_mat+test_stat_mat');
    end
    adj=test_stat_mat;
end

switch stat_type
    case 0
        %Only consider components comprising more than one node, equivalent to at
        %least one edge
        if bgl==1
            [a,sz]=components(adj);
        else
            [a,sz]=get_components(adj);
        end
        ind_sz=find(sz>1);
        sz_links=zeros(1,length(ind_sz));
        max_sz=0;
        for i=1:length(ind_sz)
            nodes=find(ind_sz(i)==a);
            if Intensity
                %Measure size as intensity
                sz_links(i)=sum(sum(adj(nodes,nodes).*test_stat_mat(nodes,nodes)))/2;
            else
                %Measure size as extent
                sz_links(i)=sum(sum(adj(nodes,nodes)))/2;
            end
            adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
            if max_sz<sz_links(i)
                max_sz=sz_links(i);
            end
        end
        
        %Subtract one to remove edges not part of a component
        %Although one is also subtracted from edges comprising a component, this is
        %compensated by the (i+1) above
        adj(~~adj)=adj(~~adj)-1;
        
    case 1 % do TFCE
        target_stats = matlab_tfce_transform(adj,'matrix');
        max_sz=max(target_stats(:));
        
    case 2 % do cNBS
        target_stats=get_constrained_stats(adj,STATS.edge_groups);
        
end

%Repeat above for each permutation
%Empirical null distribution of maximum component size
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 
if stat_type==0||stat_type==1
    null_dist=zeros(K,1); 
else
    null_dist_multiple=zeros(length(STATS.edge_groups.unique),K);
%     null_dist_multiple(:,1)=target_stats; % TODO: add this back?? Doesn't
%     seem normal NBS keeps unpermuted stat in null

%     p_approx=zeros(length(edge_groups.unique),1);
%     p_approx_null=zeros(length(edge_groups.unique),K); % TODO: incorporate this for Joe's idea 
end
p_approx=0;
%Store what is already displayed in the listbox
try pre_str=get(H,'string'); catch; end
new_str={};
%First row of test_stat is the observed test statistics, so start at the
%second row
for i=2:K+1
    
    if stat_type==0
        if ~isempty(STATS.test_stat)
            %Precomputed test statistics
            ind=ind_upper(STATS.test_stat(i,:)>STATS.thresh);
        else
            %Compute on the fly
            test_stat=NBSglm(GLM);
            ind=ind_upper(test_stat(2,:)>STATS.thresh);
        end
        if Intensity
            %Compute a test statistic matrix for intensity metric
            test_stat_mat=zeros(N,N);
            if ~isempty(STATS.test_stat)
                test_stat_mat(ind_upper)=STATS.test_stat(i,:)-STATS.thresh;
                test_stat_mat=(test_stat_mat+test_stat_mat');
            else
                test_stat_mat(ind_upper)=test_stat(2,:)-STATS.thresh;
                test_stat_mat=(test_stat_mat+test_stat_mat');
            end
        end
        adj_perm=spalloc(N,N,length(ind)*2);
        adj_perm(ind)=1; adj_perm=adj_perm+adj_perm';
        
    else % get test stats without thresholding for cNBS and TFCE - smn % TODO - Come back here!
        ind=ind_upper;
        test_stat_mat=zeros(N,N);
        if ~isempty(STATS.test_stat)
            test_stat_mat(ind_upper)=STATS.test_stat(i,:);
            test_stat_mat=(test_stat_mat+test_stat_mat');
        else
            %Compute test statistics on the fly
            test_stat=zeros(2,N*(N-1)/2);
            test_stat=NBSglm(GLM);
            test_stat_mat(ind_upper)=test_stat(i,:);
            test_stat_mat=(test_stat_mat+test_stat_mat');
        end
        adj_perm=test_stat_mat;
    end
    
    switch stat_type
        case 0
            if bgl==1
                [a,sz]=components(adj_perm);
            else
                [a,sz]=get_components(adj_perm);
            end
            ind_sz=find(sz>1);
            max_sz_perm=0;
            for j=1:length(ind_sz)
                nodes=find(ind_sz(j)==a);
                if Intensity
                    tmp=sum(sum(adj_perm(nodes,nodes).*test_stat_mat(nodes,nodes)))/2;
                else
                    tmp=sum(sum(adj_perm(nodes,nodes)))/2;
                end
                if tmp>max_sz_perm
                    max_sz_perm=full(tmp);
                end
            end
            null_dist(i-1)=max_sz_perm;
            if max_sz_perm>=max_sz
                p_approx=p_approx+1;
            end
            
        case 1 % TFCE
            target_stats_perm = matlab_tfce_transform(adj_perm,'matrix');
            max_sz_perm=max(target_stats_perm(:));
            
           
            null_dist(i-1)=max_sz_perm;
            if max_sz_perm>=max_sz
                p_approx=p_approx+1;
            end
            
        case 2 % cNBS
            null_dist_multiple(:,i-1)=get_constrained_stats(adj_perm,STATS.edge_groups); % TODO: revisit indexing - shouldn't the null include the nonpermuted stat? but the above null does not
            max_sz=0; % TODO
            max_sz_perm=0; % TODO
            p_approx=0; % TODO
            
            % visualizing results of this perm
            %         summat=zeros(size(STATS.edge_groups.groups));
            %         for n=1:length(null_dist_multiple(:,i-1))
            %             summat(STATS.edge_groups.groups==n)=null_dist_multiple(n,i-1);
            %         end
            %         figure; image(summat*1000)
        
            
    end
    
    if mod((i-1),100)==0 % SMN - only print every hundred
        % str=sprintf('|   %5d/%5d |     %4d |     %4d |   %0.3f |',...
        %v1.1.2 Changed to %6.0f to %6.1f to allow fractional component sizes
        %that arise when component size is measured with intensity.
        str=sprintf('| %5d/%5d |   %6.1f |   %6.1f |  %0.4f |',...
            i-1,K,max_sz_perm,max_sz,p_approx/(i-1));
        %Display no mare than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow;
        catch;  fprintf([str,'\n']); end
    end
    
end
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 

% con_mat will store the results in a cell array. For component-level
% results, each result is in n_cnt separate cells. Otherwise, all edge or network 
% results are in a single cell (n_cnt=1)
n_cnt=0; 
if stat_type==0
    %Determine components satisfying alpha significance threshold
    for i=1:length(sz_links)
        tmp=sum(null_dist>=sz_links(i))/K;
        if tmp<=STATS.alpha
            n_cnt=n_cnt+1;
            ind=find(adj==i);
            con_mat{n_cnt}=spalloc(N,N,length(ind)*2);
            con_mat{n_cnt}(ind)=1; 
            con_mat{n_cnt}=triu(con_mat{n_cnt},1);
            pval(n_cnt)=tmp;
        end
    end
    
elseif stat_type==1
    %Determine edges satisfying alpha significance threshold
%     maxima_sorted=sort(null_dist,'descend');
%     n_above_alpha=nbs.STATS.alpha*nbs.GLM.perms;
%     test_stat_threshold=maxima_sorted(n_above_alpha);
%     
%     tmp=target_stat>test_stat_threshold;
    
    pval=arrayfun(@(x) sum(null_dist>=x)/K, target_stats);
    con_mat{1}=pval<=STATS.alpha;
    if any(con_mat{1}(:))
        n_cnt=n_cnt+1;
    end
    
else % do cNBS
    
    % estimate pvalue for each network
    for i=1:length(STATS.edge_groups.unique)
        pval_uncorr(i)=sum(+(null_dist_multiple(i,:)>=target_stats(i)))/K;
    end
    
    % multiple comparison correction
    do_parametric_correction=1; % TODO: move - should be arg from user
    do_FWER=1; % TODO: move - should be arg from user
    
    if do_parametric_correction
        if do_FWER % Bonferroni
            pval=pval_uncorr*length(STATS.edge_groups.unique); % "adjusted p-vals"
        else % FDR
            [~,pval,~,~]=mafdr(pval_uncorr); % "adjusted p-vals"
        end
        pval(pval>1)=1;
        
    else % bootstrapping
        if do_FWER % TODO: update - this assumes interchangeability!
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
            error('Stopping - this section is under development');
        end
    end
    
    %Determine components satisfying alpha significance threshold
    surviving_group_ids=find(pval<=STATS.alpha); % TODO
    con_mat{1}=spalloc(N,N,1);
    for i=surviving_group_ids
        con_mat{1}(STATS.edge_groups.groups==i)=1;
    end
    if any(surviving_group_ids)
        n_cnt=n_cnt+1;
    end
    
    % Show components in original matrix space
%     summat=zeros(size(STATS.edge_groups.groups));
%     for i=1:length(target_stats)
%       summat(STATS.edge_groups.groups==i)=target_stats(i);
%     end
%     figure; image(summat*1000)

end

if n_cnt==0
    pval=[]; con_mat=[]; 
end