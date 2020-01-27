
%% MAIN SCRIPT

% Set up flags

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

% Set statistic type to numeric - size (extent/intensity) (0), constrained (1), or TFCE (2)
switch STATS.statistic_type
    case 'Size'
        STATS.statistic_type_numeric=0;
        % size of a component measured using extent or intensity?
        Intensity=strcmp(STATS.size,'Intensity');
    case 'TFCE'
        STATS.statistic_type_numeric=1;
    case 'Constrained' % cNBS
        STATS.statistic_type_numeric=2;
    otherwise
        error(sprintf('Statistic type %s was provided, but only ''Size'', ''TFCE'', or ''Constrained'' allowed.',STATS.statistic_type))
end

if STATS.statistic_type==2 % special multidimensional null for cNBS
    n_nulls=length(STATS.edge_groups.unique);
else
    n_nulls=1;
end
null_dist=zeros(K,n_nulls);



% Main

% get unpermuted cluster statistics
y__target=GLM.y; % save for posterity
edge_stat__target=get_univariate_stats(STATS,GLM,precomputed,1);
[cluster_stat__target,max_stat__target]=get_cluster_stats(edge_stat,STATS);
component_map=adj;

% get all other cluster statistics
for i=2:K+1
    
    GLM.y=permute_signal(GLM);
    edge_stat=get_univariate_stats(STATS,GLM,precomputed,i);
    [~,max_stat]=get_cluster_stats(edge_stat,STATS);
    null_dist(i-1,:)=max_sz;
    
    if mod((i-1),100)==0 % SMN - only print every hundred
        str=sprintf('| %5d/%5d perms complete | latest element in null: %6.1f |',i-1,K,null_dist(i,end));
        %Display no mare than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow;
        catch;  fprintf([str,'\n']); end
    end
    
end


