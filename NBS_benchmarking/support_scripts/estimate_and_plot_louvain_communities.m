function [communities,ids_new]=estimate_and_plot_louvain_communities(mat,gamma,B,nperms)
% SMN: little wrapper around community_louvain for estimating and plotting communities
% Estimates nperms
% See community_louvain.m for original (https://sites.google.com/site/bctnet/measures/list)
% See draw_atlas_boundaries for more plotting guidance
% Example: estimate_and_plot_louvain_communities(dmat_full,1.1,'negative_sym')

M0=[];
for i=1:nperms
    [communities(:,i),~]=community_louvain(mat,gamma,M0,B);
end

% get consensus
communities=mode(communities,2);

% setup for plotting
unique_communities=unique(communities);
[communities_sorted,ids_new]=sort(communities);
mat2=mat(ids_new,ids_new);

% plot
figure;
imagesc(mat2);
hold on
% plot lines
for i=2:length(unique_communities)-1
    this_id=find(communities_sorted==unique_communities(i),1);
    % plot horizontals
    p1 = [this_id-0.5,0];
    p2 = [this_id-0.5,268];
    m=plot([p1(2),p2(2)],[p1(1),p2(1)],'w-','LineWidth',2);
    m.Color(4) = 0.7; % semi-transparent

    % plot verticals
    p1 = [0,this_id-0.5];
    p2 = [268,this_id-0.5];
    m=plot([p1(2),p2(2)],[p1(1),p2(1)],'w-','LineWidth',2);
    m.Color(4) = 0.7; % semi-transparent
end
colormap(bipolar([],0.1));
caxis([-0.8,0.8]);
hold off;