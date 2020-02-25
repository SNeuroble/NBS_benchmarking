function [es,leading_edge]=get_enrichment_score(edge_groups,unique_edge_groups,test_stat,w,get_ledge)
% Calculate enrichment score for Set Enrichment Analysis - 2 complementary variable sets
%
% ---INPUT---
%  ids          : variable list
%  test_stat    : score (doesn't need to be ranked)
%  set1_ids     : positive variable set
%  set2_ids     : negative variable set
%  w            : weight
% get_ledge     : specify whether to return leading edge
% ---OUTPUT---
%  es           : enrichment score
%  leading_edge : leading edge (top variables?)
%
% References:
% Based on the following GSEA implementation:
% https://www.mathwoidss.com/matlabcentral/fileexchange/33599-gsea2
% by Wei Keat Lim
%
% Lim WK, Lyashenko E, & Califano A. (2009) Master regulatotest_stat used as  
%    breast cancer metastasis classifier. Pac Symp Biocomput, 14, 504-515.
% Carro MS*, Lim WK*, Alvarez MJ*, et al. (2010) The transcriptional 
%    netwoids for mesenchymal transformation of brain tumoutest_stat. Nature, 
%    463(7279), 318-325.

% all IDs
ids=[1:length(edge_groups)];

% check input format
if size(ids,1)==1; ids = ids'; end
if size(test_stat,1)==1; test_stat = test_stat'; end
size_tests=size(test_stat);

% augment input variables and rank by score
ids_c       = [ids;ids];
test_stat_c = [test_stat;-test_stat];
dmat        = [ones(size_tests);-ones(size_tests)];
[~,ix]      = sort(test_stat_c,'descend');
ids_c       = ids_c(ix);
test_stat_c = test_stat_c(ix);
dmat        = dmat(ix);
size_tests_combined = size(test_stat_c);

set_ids_positions=1:length(edge_groups);
set_ids_c=[edge_groups;edge_groups];
set_ids_c=set_ids_c(ix);

es=zeros(length(unique_edge_groups),1);
for i=unique_edge_groups;

    % set up set IDs
       % record overlap of variable sets and ranked list 
    %in_set1=ismember(ids_c,set1_ids);
    %in_set2=ismember(ids_c,set2_ids);
    in_set1=set_ids_c==i;
    in_set2=~in_set1;
  
    ids_for_sea = +((dmat>0 & in_set1) | (dmat<0 & in_set2));

    score_hit           = cumsum((abs(test_stat_c.*ids_for_sea)).^w);
    score_hit           = score_hit/score_hit(end);
    score_miss          = cumsum(1-ids_for_sea);
    score_miss          = score_miss/score_miss(end);
    es_all              = score_hit - score_miss;
    [es_max, id_max]    = max(es_all);
    [es_min, id_min]    = min(es_all);
    es(i)               = es_max + es_min;

    % get leading edge if specified
    if get_ledge
        isen = zeros(size_tests_combined);
        if es<0
            sup_id = id_min;
            isen(sup_id:end) = 1; % what to call this?
            leading_edge_precursor = ids_c((isen==1)&(ids_for_sea==1));
            leading_edge(i,:) = leading_edge_precursor(end:-1:1);
        else
            sup_id = id_max;
            isen(1:sup_id) = 1;
            leading_edge(i,:) = ids_c((isen==1)&(ids_for_sea==1));
            
        end
    else
        leading_edge(i)=NaN;
    end
end

end
