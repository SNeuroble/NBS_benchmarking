function map=load_atlas_mapping(n_nodes,category)

% syntax: map=load_atlas_mapping(n_nodes,category)
% n_nodes are number of nodes: 268 or 278
% category: 'subnetwork' or 'lobe'
% updated for octave: 1 is new, 2 is old, 3 is category

%load(sprintf('map%0.0f_%s.mat',n_nodes,category));
load(sprintf('map%0.0f_%s__octave.mat',n_nodes,category)); % octave
map_orig=map;
clear map
map.newroi=map_orig(:,1);
map.oldroi=map_orig(:,2);
map.category=map_orig(:,3);
