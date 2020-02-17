function data_new=get_ctx(data,ctxtype,atlascategory,varargin)
% Example: mat_new=get_ctx(icc{2,6}{51},'cortex','subnetwork')
% ctxtype = 'cortex','subcortex'
% can do from data of "atlas" dim (e.g., 268 x 1), not just matrix.
%   if not matrix,specify fourth argument as 'nomatrix'.
% mat=unstructured, un-ordered matrix


if ~isempty(varargin)
    if strcmp(varargin{1},'nomatrix')
        ismatrix=0;
    else ismatrix=1;
    end
else ismatrix=1;
end

if ismatrix
    % restructure and reorder
    if size(data,1) ~= size(data,2)
       data=structure_data(data); 
    end
    data=reorder_matrix_atlas(data,atlascategory);
end

% get atlas mappings to ctx or subctx
atlasdim=size(data,1);
mask=ones(size(data));
%trimask=logical(tril(mask,-1));


map=load_atlas_mapping(atlasdim,atlascategory);

% find the dividing index btw ctx and subctx
if strcmp(atlascategory,'subnetwork') 
    subctxids=find(strcmp(map.label,'BG'));
else % lobe
    subctxids=find(strcmp(map.label,'Cbl'));
end
subctx_startid=subctxids(1);

if strcmp(ctxtype,'cortex') || strcmp(ctxtype,'ctx') % ctx
    ids=[1:(subctx_startid-1)];
    %maskids=[subctx_startid:matdim];
elseif strcmp(ctxtype,'subcortex') || strcmp(ctxtype,'subctx') % subctx
    ids=[subctx_startid:atlasdim];
    %maskids=[1:(subctx_startid-1)];
else
    error('Specify cortex/ctx or subcortex/subctx.')
end

if ismatrix
    data_new=data(ids,ids);
else
    data_new=data(ids);
end

% mask=trimask(maskids,maskids)=0; % make mask
% % do mean/std
% meandata(1)=mean(mat(trimask));
% meandata(2)=std(mat(trimask));
