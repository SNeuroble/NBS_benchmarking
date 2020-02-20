function net_size=count_edges_between_nets(n_nodes)

load(sprintf('map%0.0f_subnetwork.mat',n_nodes));

n_nets=length(unique(map.category));

for i=1:n_nets
    for j=1:n_nets
        len_i=length(find(map.category==i));
        len_j=length(find(map.category==j));
        if i==j
            net_size(i,j)=len_i*(len_j-1)/2;
        else
            net_size(i,j)=len_i*len_j;
        end
    end
end