function [cluster_stat]=get_constrained_stats(test_stat,edge_groups)
%GET_CONSTRAINED_STATS
% Summarize stats by network defined by atlas
% TODO: add option to do scaled by weight
% TODO: consider also calculating variance across edges within group
% TODO: consider whether we should make inputs/mask be symmetric or force
% them to be upper/lower triangular (if inputs and mask are symmetric, run
% risk of double counting. Maybe best is symmetric input / lower tril mask)

% Some cleaning
if size(test_stat,1)~=size(test_stat,2)
    error('Input test statistic matrix is not square');
end

if any(test_stat-test_stat')
    error('Input test statistic matrix is not symmetric');
end

% cluster_stat_full=zeros(size(test_stat));
cluster_stat=zeros(1,length(edge_groups.unique));
% get those constrained stats - average within nets
for i=1:length(edge_groups.unique)
    group_ids=edge_groups.groups==edge_groups.unique(i);
    cluster_stat(i)=mean(test_stat(group_ids));
%     cluster_stat_full(group_ids)=cluster_stat_by_group(i);
end
