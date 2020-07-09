function fullmat=summary_to_full_matrix(summarymat,atlas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will create a full matrix where the value of summary matrix between
% network i and network j will be filled in for all edges between those
% networks in a full matrix.
%
% Example: For ICC, can end up with four 10x10 summary matrices (categorized 
% poor=1, fair=2, good=3, excellent=4) but want to represent as full 268x268 node matrix.
%
% can load atlas with load_atlas_mapping
% atlas=load('/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/general_mri/atlas/atlas mappings/map278_lobe.mat;);
% for i=1:4 % separate binary matrix for each category
%     summat_separatecat{i}=zeros(size(summarymat));
%     summat_separatecat{i}(summarymat==(i+1))=1; % poor, then fair, ... excellent
% end
% fullmat_separatecat=summary_to_full_matrix(summat_separatecat,atlas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if ~iscell(summarymat)
    summarymat={summarymat};
end

% do mapping
tmp=atlas.category(2:end)-atlas.category(1:end-1);
lobe_mapping=find(tmp>0);
lobe_mapping=[0;lobe_mapping;length(tmp)+1];
lobe_mapping=lobe_mapping+1;

for k=1:length(summarymat)
    fullmat{k}=zeros(lobe_mapping(end)-1,lobe_mapping(end)-1);
    for i=1:length(lobe_mapping)-1
        for j=1:length(lobe_mapping)-1
            i1=lobe_mapping(i);
            i2=lobe_mapping(i+1)-1;
            j1=lobe_mapping(j);
            j2=lobe_mapping(j+1)-1;
            fullmat{k}(i1:i2,j1:j2)=summarymat{k}(i,j);
        end
    end
end




% OLD: hard-coded 278 node mapping:
% lobe_mapping = [0,30,41,44,64,84,100,113,129,136,139,167,179,181,202,217,234,249,266,274,278];
% lobe_mapping=lobe_mapping+1;

% OLD: This creates a full matrix of zeros except for a single edge in the
% middle of each network that represents the value in the summary matrix.
% Not sure why this would be useful, but okay.
%
% for i=1:(length(lobe_mapping)-1)
%     nodes(i)=lobe_mapping(i)+ceil((lobe_mapping(i+1)-lobe_mapping(i))/2);
% end
% 
% for k=1:length(summarymat)
%     fullmat{k}=zeros(lobe_mapping(end)-1,lobe_mapping(end)-1);
%     for i=1:length(nodes)
%         for j=1:length(nodes)
%             fullmat{k}(nodes(i),nodes(j))=summarymat{k}(i,j);
%             if nodes(i)==nodes(j) % diag
%                diagnode2=lobe_mapping(i)+ceil((1/3)*(lobe_mapping(i+1)-lobe_mapping(i)));
%                diagnode1=lobe_mapping(i)+ceil((2/3)*(lobe_mapping(i+1)-lobe_mapping(i)));
%                fullmat{k}(diagnode1,diagnode2)=summarymat{k}(i,j);
%             end
%         end
%     end
%     fullmat{k}=tril(fullmat{k})+tril(fullmat{k},-1)'; % reflect - summary usually lower triang
% end
