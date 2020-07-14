function map=load_atlas_mapping(n_nodes,category)

% syntax: map=load_atlas_mapping(n_nodes,category)
% n_nodes are number of nodes: 268 or 278
% category: 'subnetwork' or 'lobe'
% updated for octave: 1 is new, 2 is old, 3 is category

% add atlases to path and turn off warning thatwe found the atlas in the path
[current_path,~,~]=fileparts(mfilename('fullpath')); 
addpath(genpath(current_path));

%load(sprintf('map%0.0f_%s.mat',n_nodes,category));
load(sprintf('map%0.0f_%s.mat',n_nodes,category)); % octave

