function summarize_within_surface_atlas(data,this_hemisphere)
% takes volumetric data, extracts takes mean within regions of volumetric Glasser atlas, and saves as gifti in surface Glasser atlas
% data should be volume registered to 1mm MNI Colin space
% OR 32k or 64k matrix
% thishemisphere='L', 'R', 'both'
% e.g.,
%
% %(mount test_retest)
% load('~/Documents/data/mnt/myelin_MRTool/results_reliability/icc_psd__single_Dcoeff.mat','icc_single');
% load('~/Documents/data/mnt/myelin_MRTool/results_reliability/icc_21H43M3__MRTool_psd__mask_and_var.mat','gmmask_composite')
% data=structure_data(icc_single,'mask',gmmask_composite,'ismatrix',0);
% 
% OR
%data=gifti('meanmyelin_L.func.gii'); data=data.cdata;

% check whether volume
if length(size(data))==3
    isvol=1;
else
    isvol=0;
end

if strcmp(this_hemisphere,'L');         hem=1;
elseif strcmp(this_hemisphere,'R');     hem=2;
elseif strcmp(this_hemisphere,'both');
    hem=[1,2];
    this_hemisphere=['L','R'];
    if ~isvol; error('Must R or L (not both) for surface data inputs.'); end
end

volume_atlas=load_nii('../fmri/Glasser_atlas/icbm152_t1_tal_nlin_asym_09a__masked__MMP_in_MNI_corr__improved_asym__Colin__reorient.nii.gz');
% volume is in Colin MNI space
surf_atlas=ft_read_cifti('Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
mapping_ids=surf_atlas.corticalareas_dil_final_final_areas_group_colors;
nnodes_in_hem=length(mapping_ids)/2;
% roi_ids=unique(mapping_ids(~isnan(mapping_ids)));

roi_ids_L=unique(mapping_ids(1:nnodes_in_hem));
roi_ids(:,1)=roi_ids_L(~isnan(roi_ids_L));
roi_ids_R=unique(mapping_ids(nnodes_in_hem+1:nnodes_in_hem*2));
roi_ids(:,2)=roi_ids_R(~isnan(roi_ids_R));
nrois_in_hem=length(roi_ids(:,1));

clearvars surf_atlas

% if volume,extract 64k data
hemids=[1 nnodes_in_hem nnodes_in_hem+1 nnodes_in_hem*2];
if isvol
    newdata_surface=nan(size(mapping_ids));
    for i=roi_ids(:)'
        newdata_surface(mapping_ids==i)=mean(data(volume_atlas.img==i));
    end
    
else
    offset=0;
    if length(data)==nnodes_in_hem
        hemids=[1 nnodes_in_hem]; % this isn't a true hemid; it's about where to assign data
        if hem==2; offset=nnodes_in_hem; end
    end
    
    newdata_surface=nan(size(data));
    for i=roi_ids(:,hem)'
        newdata_surface(find(mapping_ids==i)-offset)=mean(data(find(mapping_ids==i)-offset));
    end
end

clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
    
for i=1:length(hem)
    template_gifti=gifti(sprintf('meanmyelin_%s.func.gii',this_hemisphere(i)));
    template_gifti.cdata=newdata_surface(hemids((i-1)*2+1):hemids(i*2));
    save(template_gifti,sprintf('tmp_surf_inGlasserparc_%s_%s.func.gii',this_hemisphere(i),timeID))
end
