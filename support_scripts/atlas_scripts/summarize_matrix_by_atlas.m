function summary_matrix=summarize_matrix_by_atlas(mat,varargin)
% load mask, take avg within edge pairs
% IF WANT TO REORDER, MUST NOT ALREADY BE REORDERED!
% inputs: mat
% optionals: atlascategory,datacategory,saveimg,suppressimg,doedgecounts
% e.g., summarize_matrix_by_atlas(mat,'atlascategory','lobe','datacategory','ICC','reorderimg',1);

%% Parse input
p = inputParser;

defaultsaveimg=0;
defaultsuppressimg=0;
defaultreorderimg=0;
defaultdoedgecounts=0;
defaultatlascategory='subnetwork'; % options include subnetwork, lobe
defaultdatacategory='none'; % options include ICC, scandur, none

addParameter(p,'saveimg',defaultsaveimg,@isnumeric);
addParameter(p,'suppressimg',defaultsuppressimg,@isnumeric);
addParameter(p,'reorderimg',defaultreorderimg,@isnumeric);
addParameter(p,'doedgecounts',defaultdoedgecounts,@isnumeric);
addParameter(p,'atlascategory',defaultatlascategory,@ischar);
addParameter(p,'datacategory',defaultdatacategory,@ischar);

parse(p,varargin{:});

saveimg = p.Results.saveimg;
suppressimg = p.Results.suppressimg;
if saveimg; suppressimg = 0; end
reorderimg = p.Results.reorderimg;
doedgecounts = p.Results.doedgecounts;
atlascategory = p.Results.atlascategory;
datacategory = p.Results.datacategory;

clearvars p varargin


%% Setup data

% restructure and reorder
if size(mat,1) ~= size(mat,2)
   mat=structure_data(mat); 
end
if reorderimg
    mat=reorder_matrix_by_atlas(mat,atlascategory);
end

% make mask
matdim=size(mat,1);
mask=ones(size(mat));
% if size(mask,2)==1
%     if sum(mask) == size(mat,1)*(size(mat,2)-1)/2
%         [~,mask]=structure_data(mask,mat);
%     else
%         error('mask dim do not match matrix dim')
%     end
% end

map=load_atlas_mapping(matdim,atlascategory);
lobe_mapping=map.category(2:end)-map.category(1:end-1);
lobe_mapping=[0; find(lobe_mapping==1); length(lobe_mapping)+1];
lobe_mapping=lobe_mapping+1;


%% Count
% if counting edges, binarize to get category of interest

trimask=logical(tril(mask,-1));

%% Summarize lobe edges
for lobe_a=1:(length(lobe_mapping)-1)
    for lobe_b=1:(length(lobe_mapping)-1)
        matrix_ab=mat(lobe_mapping(lobe_a):(lobe_mapping(lobe_a+1)-1),lobe_mapping(lobe_b):(lobe_mapping(lobe_b+1)-1));
        trimask_ab=trimask(lobe_mapping(lobe_a):(lobe_mapping(lobe_a+1)-1),lobe_mapping(lobe_b):(lobe_mapping(lobe_b+1)-1)); %tril
        matrix_ab=matrix_ab(trimask_ab); %tril
 
        if doedgecounts % count edges in a-b; normalize by total edges
            % useful if working with binary map
            edgecount_ab=nansum(matrix_ab);
            edgecount_tot_ab=length(matrix_ab);
            summary_ab=edgecount_ab/edgecount_tot_ab;
        else % do mean edge val in a-b
            % TODO: useful if working with non-binary map...but should give same answer as above if binary
            summary_ab=nanmean(matrix_ab); %tril
        end

        summary_matrix(lobe_a,lobe_b)=summary_ab;
    end
end



%% Draw image
if ~suppressimg
    combined_summary=(tril(summary_matrix));
    draw_atlas_boundaries(combined_summary,'usingsummary',1,'atlascategory',atlascategory,'datacategory',datacategory,'saveimg',saveimg);
end


