function [cluster_stats_map,max_sz]=get_edge_components(adj,Intensity,test_stat,alpha_thresh,N,ind_upper,bgl)
%Returns connected components defined as groups of contiguous edges

%If size of component measured with intensity, create N x N matrix containing the test statistic for each edge minus the test statistic threshold (primary threshold)
if Intensity
    test_stat_mat=zeros(N,N);
    test_stat_mat(ind_upper)=test_stat-alpha_thresh;
    test_stat_mat=(test_stat_mat+test_stat_mat');
end

% Get component size in number nodes (does not return links)
if bgl==1 % check for graphics library
    [node_assignments,sz_components_in_nodes]=components(adj);
else
    [node_assignments,sz_components_in_nodes]=get_components2(adj);
end

% Get component size in number of edges and find max copmonent size
% Only components comprising more than one node, equivalent to at least one edge
component_idx=find(sz_components_in_nodes>1);
max_sz=0;
cluster_stats_map=adj;
for i=1:length(component_idx)
    nodes=find(node_assignments==component_idx(i));
    if Intensity
        %Measure size as intensity
        sz_links=sum(sum(cluster_stats_map(nodes,nodes).*test_stat_mat(nodes,nodes)))/2;
    else
        %Measure size as extent
        sz_links=sum(sum(cluster_stats_map(nodes,nodes)))/2;
    end
    % write size of component to each edge (previously unique id for each component)
    %               cluster_stats_map(nodes,nodes)=cluster_stats_map(nodes,nodes)*(i+1);
    cluster_stats_map(nodes,nodes)=cluster_stats_map(nodes,nodes)*(sz_links+1);
    if max_sz<sz_links
        max_sz=sz_links;
    end
end

%Subtract one to remove edges not part of a component
%Although one is also subtracted from edges comprising a component, this is
%compensated by the (i+1) above
% TODO: I'm concerned by this logic. Here, the smallest "sz_links" can be is "1" (similar to before, where the smallest "i" could be was "1"). So i+1 inevitably boosted components in "cluster_stats_map" up to "2", which is then decremented back to "1" in the following. Thus, this does not seem to remove components of size 1 as the following asserts. Consider whether to return this or just orig map without "+1"
cluster_stats_map(~~cluster_stats_map)=cluster_stats_map(~~cluster_stats_map)-1;

