function edge_groups=load_atlas_edge_groups(n_nodes,mapping_category)

% returns a n_nodes x n_nodes matrix where the value of each edge reflects
% the network pair it belongs to

% requires atlas loaded from load_atlas_mapping
atlas=load_atlas_mapping(n_nodes,mapping_category);

tmp=atlas.category(2:end)-atlas.category(1:end-1);
lobe_mapping=find(tmp>0);
lobe_mapping=[0;lobe_mapping;length(tmp)+1];
lobe_mapping=lobe_mapping+1;

edge_groups=zeros(lobe_mapping(end)-1,lobe_mapping(end)-1);
this_group=1;
for i=1:length(lobe_mapping)-1
    for j=1:length(lobe_mapping)-1
        if i>=j
            i1=lobe_mapping(i);
            i2=lobe_mapping(i+1)-1;
            j1=lobe_mapping(j);
            j2=lobe_mapping(j+1)-1;
            edge_groups(i1:i2,j1:j2)=this_group;
            this_group=this_group+1;
        end
    end
end

edge_groups=tril(edge_groups)+tril(edge_groups,-1)';
