function CC = get_component_IDs(adj,comps,comp_sizes)
% (SMN) Made to return same measures of component size as bwconncomp for use with
% NBS TFCE (bwconncomp developed for images, get_components developed for graphs)
% cc.PixeldxList = IDs of all edges for each cluster
% cc.NumObjects = number of components

if iscell(comps)
    comps=comps{1};
    comp_sizes=comp_sizes{1};
end

% Find all nonzero edges for subsequent shuffling into components - only
% upper triangle - % TODO: should return upper tri but should reflect in
% case lower
adj_triu=triu(adj,1);
[idx]=find(adj_triu);
[xi,yi]=ind2sub(size(adj_triu),idx);

% Measure size as extent
CC = struct(...
    'NumObjects', [], ...
    'PixelIdxList', []);

ind_sz=find(comp_sizes>1);
CC.NumObjects=length(ind_sz);
sz_links=zeros(1,length(ind_sz));
max_sz=0;
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==comps);
    CC.PixelIdxList{i}=idx(ismember(xi,nodes) | ismember(yi,nodes)); % this could prob be sped up
%     size_comp{i}=sum(sum(adj(nodes,nodes)))/2;
%     adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
%     adj=tril(adj,-1);
end
