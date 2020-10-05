function [struct_data,mask]=structure_data(unstruct_data,varargin)
% restructure data that has been pulled out of a mask
% varargin is for structure type. Use defaults if voxelwise or matrix. Otherwise, if 1-D, specify 'none'.
%
% Ex: from glm data
% load glm_data; load gmmask_composite
% ps = get_multiseed_vals(glm_data,'p'); % another example: unstruct_data=pcc_icc_summary{3,1}{1};
% [fdr_ps,qs]=fdr_correct_multiseed_pvals(ps);
% thresh_ps=threshold_multiseed_pvals(qs);
% for i=2:size(thresh_ps,2)
%    struct_data{i-1}=structure_data(thresh_ps{3,i},'mask',gmmask_composite);
% end
% THE FOLLOWING ONLY SAVES FOR SITE 3:
% save_mat_as_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz',struct_data{3},'temp.nii.gz') % visualize
%
%
% tips:
% - provide 3rd argument 'none' if using a array mask and don't want it
% converted to matrix
% - provide 2rd argument 'upper' if giving a array mask and don't want it
% converted to matrix
%
% see run_glm line 30 for obtaining gmmask_composite
%   use make_glm_imgs after to convert to structured niis
% for quick view: nii=load_nii('temp.nii.gz'); view_nii(nii)


%% Parse input
p = inputParser;

defaultmask=0;
defaultismatrix=1;
defaulttriangleside='lower';

addParameter(p,'ismatrix',defaultismatrix,@isnumeric);
addParameter(p,'triangleside',defaulttriangleside,@ischar);
addParameter(p,'mask',defaultmask,@isnumeric);

parse(p,varargin{:});

ismatrix = p.Results.ismatrix;
triangleside = p.Results.triangleside;
mask = p.Results.mask;

using_defaults=p.UsingDefaults;

clearvars p varargin

if ~(ismatrix==0 || ismatrix==1)
    error('ismatrix must be 0 or 1');
end

if ~any(strcmp({'upper','lower'},triangleside))
    error('triangleside must be upper or lower');
end

%% Do masking

% Check whether user input mask
if any(strcmp('mask',using_defaults))
    
    % No mask provided - create mask
    
    if ismatrix==0
        error('No mask provided but vector specified, so mask is assumed to be a vector of ones. Why un-mask if you''re returning the same vector?')
    end
    
    r=int16(roots([1 -1 -2*length(unstruct_data)]));
    
    if any(strcmp('triangleside',using_defaults))
        warning('STRUCTURE_DATA:ASSUME_LOWER_TRI','No mask or triangle side given. ASSUMING MASK IS LOWER TRIANGULAR MATRIX. Not appropriate for NBS corrections.')
    end
    
    if strcmp(triangleside,'lower')
        mask=ones(r(1),r(1));
        mask=tril(mask,-1);
    elseif strcmp(triangleside,'upper')
        mask=ones(r(1),r(1));
        mask=triu(mask,1);
    end
    
    
else
    
    if sum(+mask(:)) ~= length(unstruct_data)
        error('Mask - data dimension mismatch')
    end
    
    if any(size(mask)==1) % mask is array; have to convert to matrix
        
        if ismatrix
            
            r=int16(roots([1 -1 -2*length(mask)]));
            
            if any(strcmp('triangleside',using_defaults))
                warning('Mask is array and specified conversion to matrix, but no triangle side given. ASSUMING MASK IS LOWER TRIANGULAR MATRIX. Not appropriate for NBS corrections.')
            end
            
            if strcmp(triangleside,'lower')
                mask2=ones(r(1),r(1));
                mask2=tril(mask2,-1);
            elseif strcmp(triangleside,'upper')
                mask2=ones(r(1),r(1));
                mask2=triu(mask2,1);
            end
            
            t=find(mask2);
            struct_mask=+mask2;
            struct_mask(t)=mask;
            mask=struct_mask;
            
        end
    end
end

t=find(mask);
struct_data=+mask;
struct_data(t)=unstruct_data;

if any(strcmp('mask',using_defaults)) && ismatrix
    % reflect entries across diag bc matrix - we know the mask was
    % triangular bc we created it (ie, zeros on upper or lower)
    struct_data=struct_data+struct_data';
end

% to subsequently mask by threshold: d=data*(+data>0.74);

end