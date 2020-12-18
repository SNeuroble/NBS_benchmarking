function [tfced,comps,comp_sizes] = matlab_tfce_transform(img,varargin)
% Generalized for connectivity using code from https://github.com/markallenthornton/MatlabTFCE
% Matrix TFCE parameter defaults from https://github.com/MRtrix3/mrtrix3/blob/master/cmd/connectomestats.cpp
% 
% Specify: 'matrix' for matrix, 'image' for 2 or 3D image
% TODO: other weightings as optional args
%
%OLD MATLAB_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = matlab_tfce_transform(img,H,E,C,dh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity/"neighborhood" (if image)
%   -- dh size of steps for cluster formation
%   -- tfced output tfce statistic (upper triangular)

ndim=length(size(img));
dh = 0.1;
if nargin > 1
    if strcmp(varargin{1},'image')
        is_graph=0;
        if ndim==3 % defaults only tested for fMRI data?
            % defaults parameters for 3D from matlab_fce (Smith et al.) - TODO: move to user-defined section
            % TODO: these param are for 3D - add check for 2D (C=8) (default to full neighborhood)
            H = 2;
            E = 0.5;
            C = 26; % Smith et al. also tested these parameters for C=6
        elseif ndim==2 && ~any(size(img)==1) % defaults only tested for surface/TBSS?
            H = 2;
            E = 1;
            C = 8;
        else
            error('only default params for 2d and 3d')
        end
    elseif strcmp(varargin{1},'matrix')
        is_graph=1;
        if ndim==2
            % default parameters from Mrtrix2 connectomestats.cpp function (see TFCE_*_DEFAULT)- cited Vinokur et al., 2015 in the code and Baggio et al., 2018 in the docs - values differ from the references - close to the values from Baggio et al.: "We therefore recommend E parameter values of 0.5 (combined with H parameter values between 2.25 and 3) or 0.75 (combined with H parameters between 3 and 3.5)"
            H = 3.0;
            E = 0.4;
        else
            error('only default params for 2D matrices currently supported')
        end
    end
else
    error('Must provide input type (image or matrix)')
end
    

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive elements (voxels/edges greater than first threshold)
n_elements = length(img(:));

% get connected components
%   cc.PixeldxList = IDs of all edges for each cluster in a cell (for each dh)
%   cc.NumObjects = number of components (for each dh)
if is_graph
    % network components - ONLY works for upper triangular input
    [comps,comp_sizes] = arrayfun(@(x) get_components2(bsxfun(@ge,img,x),1), threshs);
    cc = arrayfun(@(x) get_component_IDs(bsxfun(@ge,img,threshs(x)),comps{x},comp_sizes{x}), [1:ndh]);
else
    % image components
    cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,img,x),C), threshs);
end

% calculate TFCE statistic
vals = zeros(n_elements,1);
for h = 1:ndh
    clustsize = zeros(n_elements,1);
    ccc = cc(h);
    elempercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:cc(h).NumObjects
        clustsize(ccc.PixelIdxList{c}) = elempercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(img));
tfced(:) = vals.*dh;

comps=comps{1};
comp_sizes=comp_sizes{1};

end


% checked TFCE calculation manually:
% given 33 contiguous edges with value 1
% dh=0.1 so 10 levels from 0 to 1 (not incl. 0)
% height of each edge at each dh = .1*i, weight H=3
% extent at each dh = 33, weighting E = 0.4
%for i=1:10; c(i)=33^.4*(.1*i)^3; end
%sum(c)*.1

% old code for getting numobjects after the fact
% for h = 1:ndh
%     cc(h).PixelIdxList={cc2(h)};
%     cc(h).PixelIdxList={cc2(h)};
% end
% NumObjects=arrayfun(@(x) length(cc(x).PixelIdxList), 1:length(cc));
% cc.NumObjects=NumObjects;

