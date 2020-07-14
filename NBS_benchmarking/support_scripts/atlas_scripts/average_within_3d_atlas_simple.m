function data_atlas = average_within_3d_atlas_simple(data,atlas)
% average within atlas
% can load atlas with load_atlas_mapping
% atlas is typically '~/Google\ Drive\ File\ Stream/My\ Drive/Academic/Lab/Steph\ -\ Lab/Misc/Software/data/bioimagesuite/images/shenetal_neuroimage2013/shen_1mm_268_parcellation.nii.gz'

if ~isequal( size(data) , size(atlas))
    error(['Sizes don''t match. Size data is ',sprintf('%d ',size(data)),', size atlas is ',sprintf('%d ',size(atlas))]);
end

if min(atlas(atlas~=0)) < 0
    warning('Negative entries found in atlas. Proceeding anyways.')
end

data_atlas=atlas;
for i = min(atlas(atlas~=0)):max(atlas(:))
        ind=find(atlas==i);
        data_atlas(ind)=mean(data(ind));
end
