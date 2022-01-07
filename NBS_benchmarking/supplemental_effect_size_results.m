function supplemental_effect_size_results(supp_summary_type)
% Run analyses for supplemental results examining the spatial extent of ground truth 
% effects and applicability of the Shen268-10 community partition

% do_paired_task_contrasts: compare effect sizes (well, t-statistics at this point) across tasks
% do_subnetwork_anova
% do_community_detection
% calc_ground_truth_pvals
% calc_ground_truth_pvals_mv


%% PARSE USER INPUT
% p = inputParser;
% addOptional(p,'do_paired_task_contrasts',0);
% addOptional(p,'do_subnetwork_anova',0);
% addOptional(p,'do_community_detection',0);
% addOptional(p,'calc_ground_truth_pvals',0);
% addOptional(p,'calc_ground_truth_pvals_mv',0);
% addOptional(p,'get_clusters',0);
% parse(p,varargin{:});
% 
% do_paired_task_contrasts=p.Results.do_paired_task_contrasts;
% do_subnetwork_anova=p.Results.do_subnetwork_anova;
% do_community_detection=p.Results.do_community_detection;
% calc_ground_truth_pvals=p.Results.calc_ground_truth_pvals;
% calc_ground_truth_pvals_mv=p.Results.calc_ground_truth_pvals_mv;
% get_clusters=p.Results.get_clusters;

%% PARAMS

% check whether to use rest trimmed to shortest task (EMOTION)
use_trimmed_rest=1; % 0 by default
trimmed_rest_run='REST_176frames';

% plot params
save_figs=1;

margin=.01;
caxis_single=[-25,25];
caxis_contrast=[-7,7];
yaxislim=[-2.5,2.5];
%     plot_dim=[0,0,pp.fig_width_1plt_norm_big, pp.fig_height_1plt_norm_big]

%% SETUP

% [current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarking
% addpath(genpath(['../',current_path])); % need to add the complete script dir to get the config files
addpath(genpath('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking'));
addpath('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/general/') % fdr corr
setpaths;
do_fpr=0; % needed for setparams
setparams_summary;

grsize=grsize_default; % needed for set_datetimestr
summary_type='visualize_gt'; % needed for set_datetimestr
stat_type=stat_type_gt; % needed for set_datetimestr
set_datetimestr_and_files;
% use local copies instead of remote for quick access
% data_dir='/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/tmp/';
file_prefix='/Volumes/GoogleDrive/My Drive/Lab/NBS_benchmarking/NBS manuscript/figs/intermediate files/supplemental/';

% set task names (overwriting setup above) and switch things around if using trimmed
all_tasks={'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM','REST'};
rest_id=contains(all_tasks,'REST');
if use_trimmed_rest
    all_tasks{rest_id}=trimmed_rest_run;
end
rest_str=all_tasks{rest_id};

n_tasks=length(all_tasks);
% stat_type={'Size_Extent'};

% task duration in num frames
n_frames.EMOTION=176; % shortest task
n_frames.GAMBLING=253;
n_frames.LANGUAGE=316;
n_frames.MOTOR=284;
n_frames.RELATIONAL=232;
n_frames.SOCIAL=274;
n_frames.WM=405; % longest task
n_frames.REST=1200;
n_frames.REST2=1200;
n_frames.REST_176frames=176;

task_duration=[];
% task_duration_min=[2+16/60, 3+12/60, 3+57/60, 3+34/60, 2+56/60, 3+27/60, 5+01/60, 14+33/60];
% task duration in min:sec: 2x 2:16, 3:12, 3:57, 3:34, 2:56, 3:27, 5:01
% REST1: 2x 14:33


%% MAIN

% pre-load data and make triu mask for certain procedures: mean dcoeff, partition, triangle mask
if any(strcmp(supp_summary_type,{'do_subnetwork_anova','do_community_detection','get_clusters','calc_ground_truth_pvals'}))
    
%     load('/Volumes/GoogleDrive/My Drive/Steph-Lab/Misc/Software/scripts/Matlab/myscripts/tmp/ground_truth__EMOTION_Size_Extent_03182021_2251.mat');
    load('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/tmp/mean_dcoeff.mat')
    edge_stats=mean_dcoeff.Parametric_FDR;

    % make tril matrix for plotting and full matrix for partitioning
    triu_msk=triu(true(268),1);
    dmat=+triu_msk;
    dmat(triu_msk)=edge_stats;
    dmat_full=dmat+dmat';
    
    % load original mapping table (for community detection comparison) and edge groups matrix (for spatial perm)
    shen268_mapping=load_atlas_mapping(268,'subnetwork');
    shen268_mat=load_atlas_edge_groups(268,'subnetwork');
    shen268_triu=shen268_mat(triu_msk);
    shen268_nodecommunities=shen268_mapping.category;
    shen268_n_nodecommunities=length(unique(shen268_nodecommunities));
    shen268_n_subnets=length(unique(shen268_triu));
    
    % special masks for within and between nets
    unique_within_subnets=unique(diag(shen268_mat));
    shen268_triu_within_msk=ismember(shen268_triu,unique_within_subnets);
    shen268_triu_btw_msk=~shen268_triu_within_msk;
    
    ids_reordered_within_first=[find(shen268_triu_within_msk);find(shen268_triu_btw_msk)];
    
end
    
switch supp_summary_type

%% Test for differences between subnetworks via nonparametric ANOVA
% with and without within-community structure
case 'do_subnetwork_anova'
    
    anova_file_prefix=[file_prefix,'network variance/'];
    mkdir(anova_file_prefix)
    
    % setup for perms (TODO: may have to *not* preallocate for parfor?)
    n_perms=10000;
    F_perm=zeros(n_perms);
    F_perm_within=zeros(n_perms);
    F_perm_btw=zeros(n_perms);
    F_perm_within__keeping_within_structure=zeros(n_perms);
    F_perm_within__keeping_btw_structure=zeros(n_perms);
        
    % Target F-stat for all subnetworks
    [p_parametric,tbl,stats_parametric]=anova1(edge_stats(ids_reordered_within_first),num2str(shen268_triu(ids_reordered_within_first)));
    F_target=tbl{2,5};
    
    
    % plot verticals

%     hAx=gca; 
%     xtk=hAx.XTick; 
%     hold on;
    pause(1) % not too fast so we access the correct plot 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.8, 0.5])
    hold on;
    m=plot([10.5,10.5],[-3,3],'k-','LineWidth',0.5);
    m.Color(4) = 0.4; % semi-transparent
    m=plot([0,100],[0,0],'k-','LineWidth',0.5);
    m.Color(4) = 0.4; % semi-transparent
    ylim(yaxislim);
    
    if save_figs; print(gcf,[anova_file_prefix,'anova_boxplots'],'-dpng','-r300'); end
    
% hAx=gca;                                   % retrieve the axes handle
% xtk=hAx.XTick;                             % and the xtick values to plot() at...
% hold on
% hL=plot(xtk,rand(size(xtk)),'k-*')
% hL=plot([10.5,10.5],[-3,3]);
    

    % Target F-stat for subnetworks, separately within-community subnetworks and and between-community subnetworks
    [p_within_parametric,tbl_within,stats_within_parametric]=anova1(edge_stats(shen268_triu_within_msk),shen268_triu(shen268_triu_within_msk),'off');
    F_target_within=tbl_within{2,5};
%     ylim(yaxislim); print(gcf,[file_prefix,'anova_boxplots_within_communities'],'-dpng','-r300'); 
    
    [p_btw_parametric,tbl_btw,stats_btw_parametric]=anova1(edge_stats(shen268_triu_btw_msk),shen268_triu(shen268_triu_btw_msk),'off');
    F_target_btw=tbl_btw{2,5};
%     ylim(yaxislim); print(gcf,[file_prefix,'anova_boxplots_btw_communities'],'-dpng','-r300'); 
    
    
    
    
    parfor k=1:n_perms

        % Random parcellation
        ids_perm=randperm(268);
        dmat_full_perm=dmat_full(ids_perm,ids_perm);
        edge_stats_perm=dmat_full_perm(triu_msk);

        [~,tbl_perm,~]=anova1(edge_stats_perm,shen268_triu,'off');
        F_perm(k)=tbl_perm{2,5};
        [~,tbl_perm_within,~]=anova1(edge_stats_perm(shen268_triu_within_msk),shen268_triu(shen268_triu_within_msk),'off');
        F_perm_within(k)=tbl_perm_within{2,5};
        [~,tbl_perm_btw,~]=anova1(edge_stats_perm(shen268_triu_btw_msk),shen268_triu(shen268_triu_btw_msk),'off');
        F_perm_btw(k)=tbl_perm_btw{2,5};

        % Random reassignment of edges, keeping within-community edges in a within-community subnetwork
        ids_perm_within=shen268_triu_within_msk(randperm(length(shen268_triu_within_msk)));
        ids_perm_btw=shen268_triu_btw_msk(randperm(length(shen268_triu_btw_msk)));

        [~,tbl_perm_within,~]=anova1(edge_stats(ids_perm_within),shen268_triu(shen268_triu_within_msk),'off');
        F_perm_within__keeping_within_structure(k)=tbl_perm_within{2,5};
        [~,tbl_perm_btw,~]=anova1(edge_stats(ids_perm_btw),shen268_triu(shen268_triu_btw_msk),'off');
        F_perm_btw__keeping_btw_structure(k)=tbl_perm_btw{2,5};


    end

    % Test whether subnetworks show substantial modularity relative to random parcellation
    F_pval=sum(+(F_perm > F_target))/n_perms;

    % Test whether subnetworks are modular relative to a random reassignment of edges,
    % keeping within-community edges in a within-community subnetwork
    % (gets at whether all within-community edges could be pooled together)
    F_pval_within=sum(+(F_perm_within > F_target_within))/n_perms;
    F_pval_btw=sum(+(F_perm_btw > F_target_btw))/n_perms;
    F_pval_within__keeping_within_structure=sum(+(F_perm_within__keeping_within_structure > F_target_within))/n_perms;
    F_pval_btw__keeping_btw_structure=sum(+(F_perm_btw__keeping_btw_structure > F_target_btw))/n_perms;
        
        %{
        
            %{
                for g=1:shen268_n_subnets
                this_msk=find(shen268_triu==g);
                edges_in_net_perm=edge_stats_perm(this_msk);
%                 mean_perm(g,k)=mean(edges_in_net_perm);
%                 var_perm(g,k)=var(edges_in_net_perm);
            end
            %}
        
%         p_right=sum(+(mean_target(g)<mean_perm(g,:)))/n_perms;
%         p_left=sum(+(mean_target(g)>mean_perm(g,:)))/n_perms;
%         if p_right<p_left
%             target_is_different_pval(g)=p_right;
%         else target_is_different_pval(g)=p_left;
%         end


        for g=1:shen268_n_subnets
            target_is_less_variable_pval(g)=sum(+(var_target(g)>var_perm(g,:)))/n_perms;
        end
        

        [~, ~, ~, target_is_different_p_corr]=fdr_bh(target_is_different_pval,0.05);
        n_different=sum(+target_is_different_p_corr<0.05);
        [~, ~, ~, target_is_less_variable_p_corr]=fdr_bh(target_is_less_variable_pval,0.05);
        n_more_homogeneous=sum(+target_is_less_variable_p_corr<0.05);
    %     [~,p_fdr]=mafdr(target_is_different_pval);

        test_group=10;
        figure;
        histogram(var_perm(test_group,:)); hold on; histogram(var_target(test_group));
        %}

    
    % 1.2. Pairwise tests to determine which subnetworks differ?
    

%% Re-estimate communities (louvain) from edge-level ground truth data
% Partition of NEGATIVE dcoeff (bc negative within-community subnets)
case 'do_community_detection'
    
    community_file_prefix=[file_prefix,'communities/'];
    mkdir(community_file_prefix)
    
    
    % plot dcoeffs in original communities
    figure; draw_atlas_boundaries(-dmat_full);
    colormap(bipolar([],0.1));
    caxis([-0.8,0.8]);
    if save_figs; print(gcf,[community_file_prefix,'shen_partition'],'-dpng','-r300'); end  
    
    % Estimate and plot partition via louvain (Bassett)
    louvain_thresholds=[3,4]; % [3,4,5];
    nperms_for_community_est=10000;
    cmap_overlap_map=[linspace(0, 0.98, 256).^3', linspace(0, 0.98, 256).^3', linspace(0, 0.98, 256).^1.5'];
    
    for this_louvain_thresh=louvain_thresholds
        
        % estimate communities from consensus (~7 sec for 1,000 perms)
        [louvain_communities,~]=estimate_and_plot_louvain_communities(-dmat_full,this_louvain_thresh,'negative_asym',nperms_for_community_est);
        n_louvain_communities=max(louvain_communities);
        caxis([-0.8,0.8]);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.5])
        if save_figs; print(gcf,[community_file_prefix,'louvain_thresh',num2str(this_louvain_thresh,2)],'-dpng','-r300'); end
        
        % plot dcoeffs node-level communities
        figure;  
        colormap(cmap_overlap_map)
        subp(1)=subplot(1,shen268_n_nodecommunities+2,1);
        imagesc(shen268_nodecommunities)
        set(gca,'XTickLabel',[]);
        subp(2)=subplot(1,shen268_n_nodecommunities+2,2);
        imagesc(louvain_communities)
        text(1, length(louvain_communities)+20, sprintf('%1.0f',n_louvain_communities),'HorizontalAlignment','Center')
        axis off;
        colormap(subp(1),'jet'); colormap(subp(2),'jet')

        for i=1:shen268_n_nodecommunities
            new_assignments=louvain_communities(shen268_nodecommunities==i);
            unique_assignments=unique(new_assignments);

            overlap_map=zeros(size(louvain_communities)); % will highlight the most overlapping communities
            for j=1:2
                overlapping_community(i,j)=mode(new_assignments);
                overlap_map=overlap_map+(louvain_communities==overlapping_community(i,j))*(4-j);
                overlap_norm(i,j)=sum(+(new_assignments==overlapping_community(i,j)))/sum(+(louvain_communities==overlapping_community(i,j))); % proportion overlap
    
                new_assignments(new_assignments==overlapping_community(i,j))=NaN; % remove this mode before searching for the next mode
            end

            subplot(1,shen268_n_nodecommunities+2,i+2);
            imagesc(overlap_map)
            if i==shen268_n_nodecommunities/2
                labelstr='Overlap amongst all new nodes';
            else
                labelstr='';
            end
            text(1, length(louvain_communities)+20, sprintf(['%1.0f\n',labelstr,'\n%0.3f'],i,overlap_norm(i,1)),'HorizontalAlignment','Center')
            axis off; box on; %set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        end

        orig_net_labels=unique(shen268_mapping.label,'stable');
        group_for_each_orig_net=[orig_net_labels';num2cell(overlapping_community(:,1))'];
        
%         fprintf('Amt overlap: %0.4f\n',amt_overlap(:,1));

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5])
        if save_figs; print(gcf,[community_file_prefix,'louvain-shen_overlap_thresh',num2str(this_louvain_thresh,2)],'-dpng','-r300'); end
    end
    %{
    % Estimate partition using https://www.mathworks.com/matlabcentral/fileexchange/52699-normalized-cut-segmentation-using-color-and-texture-data
    sNcut = 100; %50; % The smallest Ncut value (threshold) to keep partitioning 
    sArea = 15; % The smallest size of area (threshold) to be accepted as a segment
    [segments segment_ids segment_ncut_vals] = NcutPartition([1:268], dmat_full, sNcut, sArea, 'test');
    length(segments)
    % plot arranged by new segments
    ids_new=cell2mat(segments);
    dmat_full2=dmat_full(ids_new,ids_new);
    figure; draw_atlas_boundaries(dmat_full2);
    colormap(bipolar([],0.1));
    caxis([-0.8,0.8]);
    
    %}  



%% Threshold gt map and get clusters
case 'get_clusters'
    % Note: Only get multiple clust (and higher for neg = within-net) high thresh, max n clusters @ thresh=1.1
    
    cluster_file_prefix=[file_prefix,'clusters/'];
    mkdir(cluster_file_prefix)
    
    thresholds=[0.2,0.5,0.8,1.0,1.1]; % Cohen's d CDT
    new_way=1;
    
    if new_way
        for this_cdt=thresholds
    
        sign_multiplier=[1,-1];
        figure;
        for i=1:2 % each sign
            [cluster_stats_map{i},max_sz(i)]=get_edge_components(+(sign_multiplier(i)*dmat>this_cdt),0,edge_stats,this_cdt,268,find(triu_msk),0);
            unique_clust{i}=unique(cluster_stats_map{i}(cluster_stats_map{i}>0.5)); % 0.5=island (single edge)
            n_clust(i)=length(unique_clust{i});
            for j=1:n_clust(i)
                cluster_stats_map{i}(cluster_stats_map{i}==unique_clust{i}(j))=j;
            end
        end
        
        cluster_stats_map_bipolar=cluster_stats_map{1}-cluster_stats_map{2};
        draw_atlas_boundaries(cluster_stats_map_bipolar')
        
        colormap(bipolar([],0.1));
        max_caxis=max(abs(cluster_stats_map_bipolar(:)));
        caxis([-max_caxis,max_caxis]);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 1])
        if save_figs; print(gcf,[cluster_file_prefix,'cluster_thresh_d',num2str(this_cdt,2),'_tmp'],'-dpng','-r300'); end
    end
    else
        for this_cdt=thresholds

            sign_multiplier=[1,-1];
            figure;
            for i=1:2 % each sign
                [cluster_stats_map{i},max_sz(i)]=get_edge_components(+(sign_multiplier(i)*dmat>this_cdt),0,edge_stats,this_cdt,268,find(triu_msk),0);
                unique_clust{i}=unique(cluster_stats_map{i}(cluster_stats_map{i}>0.5));
                n_clust(i)=length(unique_clust{i}); % 0.5=island
                for j=1:n_clust(i)
                    cluster_stats_map{i}(cluster_stats_map{i}==unique_clust{i}(j))=j;
                end
                subplot(2,1,i); draw_atlas_boundaries(cluster_stats_map{i}')
            end

            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 1])
            if save_figs; print(gcf,[cluster_file_prefix,'cluster_thresh_d',num2str(this_cdt,2),'_tmp'],'-dpng','-r300'); end
        end
    end

%% Get edge-level pvals for ground truth maps: d->t->pvals
case 'calc_ground_truth_pvals'

    sig_gt_file_prefix=[file_prefix,'significant ground truth/'];
    mkdir(sig_gt_file_prefix)
    
    load('/Volumes/GoogleDrive/My Drive/Lab/Misc/Software/scripts/Matlab/myscripts/NBS_benchmarking/tmp_imgs/dcoeff.mat')
    p_thresh=0.025; %p=0.05 2-sided
    
    do_nets=0;
    do_appx=0;
    
    if do_nets
        scale_idx=2;
        scale_str='edges';
    else
        scale_idx=1;
        scale_str='networks';
    end
    
    if ~do_appx
        n_subs_per_task=[1022,1057,1021,1058,1016,1027,1058];
        for i=1:7
            n_subs=n_subs_per_task(i);
            tstat=dcoeff.Parametric_FDR{scale_idx}(:,i)*sqrt(n_subs);
            tstat_task(:,i)=tstat;
    
            ids_pos=tstat>0;
            pval(ids_pos)=1-tcdf(tstat(ids_pos),n_subs-1);
            pval(~ids_pos)=1-tcdf(-tstat(~ids_pos),n_subs-1);
            [~, ~, ~, p_corr_task(:,i)]=fdr_bh(pval,p_thresh); % TODO: fdr_bh rather than mafdr in case bioinformatics toolbox not available
            
            sig_mat_task(:,i)=p_corr_task(:,i)<p_thresh; % two-sided
            prop_sig_task(i)=mean(+sig_mat_task(:,i));
            
            % NEW: repeated for  bonf
            p_corr_task_bonf(:,i)=pval*length(pval(:));
            sig_mat_task_bonf(:,i)=p_corr_task_bonf(:,i)<p_thresh; % two-sided
            prop_sig_task_bonf(i)=mean(+sig_mat_task_bonf(:,i));
            
        end
        
        avg_sig_mat_task=mean(sig_mat_task,2);
        avg_prop_sig_task=mean(prop_sig_task);
        
        
        % calc from avg tstat
        tstat_avg=mean(tstat_task,2);
        pval(ids_pos)=1-tcdf(tstat_avg(ids_pos),n_subs-1);
        pval(~ids_pos)=1-tcdf(-tstat_avg(~ids_pos),n_subs-1);
        [~, ~, ~, p_corr]=fdr_bh(pval,p_thresh); % TODO: fdr_bh rather than mafdr in case bioinformatics toolbox not available

        sig_mat_from_mean_tstat=p_corr<p_thresh; % two-sided, from avg tstat
        prop_sig=mean(+sig_mat_from_mean_tstat);

        % define sig mat for plotting
        sig_mat=avg_sig_mat_task;
            
%         sig_mat=prod(sig_mat_it,2);
%         prop_sig=mean(+sig_mat);

        % NEW: repeated for bonf
        
        avg_sig_mat_task_bonf=mean(sig_mat_task_bonf,2);
        avg_prop_sig_task_bonf=mean(prop_sig_task_bonf);
        
        % calc from avg tstat
        tstat_avg=mean(tstat_task,2);
        pval(ids_pos)=1-tcdf(tstat_avg(ids_pos),n_subs-1);
        pval(~ids_pos)=1-tcdf(-tstat_avg(~ids_pos),n_subs-1);
        p_corr_bonf=pval*length(pval(:));

        sig_mat_from_mean_tstat_bonf=p_corr_bonf<p_thresh; % two-sided, from avg tstat
        prop_sig_bonf=mean(+sig_mat_from_mean_tstat_bonf);

        % define sig mat for plotting
        sig_mat_bonf=avg_sig_mat_task_bonf;



    else
        % approximation based on avg dcoeff and avg num subs
        n_subs=1037;
        tstat=edge_stats*sqrt(n_subs);
    
        ids_pos=tstat>0;
        pval(ids_pos)=1-tcdf(tstat(ids_pos),n_subs-1);
        pval(~ids_pos)=1-tcdf(-tstat(~ids_pos),n_subs-1);
        [~, ~, ~, p_corr]=fdr_bh(pval,p_thresh); % TODO: why using fdr_bh rather than mafdr??

        sig_mat=p_corr<p_thresh; % two-sided
        prop_sig=mean(+sig_mat);
    end

    if do_nets
         triu_msk_net=triu(true(10)); % confirmed that this is a triu situation
        sig_mat_struct=structure_data(sig_mat,'mask',triu_msk_net);
%         tstat_avg_struct=structure_data(tstat_avg,'mask',triu_msk_net);
        imagesc(sig_mat_struct');
%         imagesc(tstat_avg_struct');
        axis('square')
%         ids_triu=find(triu_msk);
%         summary_to_full_matrix(sig_mat,shen268_mapping);
    else
        dmat=+triu_msk;
        dmat(triu_msk)=sig_mat; % two-sided
    %     dmat(triu_msk)=sig_mat'.*edge_stats; % mask dcoeffs by sig
        figure; draw_atlas_boundaries(dmat');
    %     colormap(bipolar([],0.1));
    end
    
    colormap(gray)
    caxis([0,1]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 1])
    if save_figs; print(gcf,[sig_gt_file_prefix,'fdr_corr_p_0.025_',scale_str,'_tmp'],'-dpng','-r300'); end
   
    

%% IN PROGRESS: Get multivariate pval for ground truth map
% Mahal d -> HotellingT2 -> F -> N for which this T2 is critical (p=0.05)
%rule of thumb: in univ case, X subs for small, X subs for med, X subs for large
case 'calc_ground_truth_pvals_mv'

    % hard coded results from gambling task
    d_omnibus=21.6472; % from working memory task, see %nbs_net_mv.mahald;
    n=1061;
    p=55;
    
    % d to T2
    T2=(n/(n-1))*d_omnibus^2; % https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution#Hotelling_t-squared_statistic
    F=((n-p)/(p*(n-1)))*T2; % F appx: https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.3  see also: https://www.mathworks.com/matlabcentral/fileexchange/2844-hotellingt2
    pval=1-fcdf(F,p,n-p);
    


%% Do paired task contrasts, with and without trimmed rest
case 'do_paired_task_contrasts'

    task_contrast_file_prefix=[file_prefix,'task contrasts/'];
    mkdir(task_contrast_file_prefix)
    
    if use_trimmed_rest; trimmed_str='_trimmed'; else trimmed_str=''; end
        
    % Load data for each task
    for i=1:n_tasks
       this_task=all_tasks{i};
       task.(this_task)= load([output_dir,'ground_truth__',this_task,'_',stat_type_gt,'_',date_time_str_ground_truth.(this_task),'.mat']);
       task_duration=[task_duration,n_frames.(this_task)];
    end

    % Check rest is last
    if ~contains(all_tasks{n_tasks},'REST')
        error('The rest of the script expects REST to be last...')
    end

    % Calculate mean across tasks
    msk=triu(true(268),1);
    it=1;
    for this_task=1:n_tasks
        d_all(:,this_task)=task.(all_tasks{this_task}).edge_stats;
        if this_task~=n_tasks
            d_all_no_rest(:,it)=d_all(:,this_task);
            it=it+1;
        end
    end

    d_all_mean=mean(d_all,2);
    d_all_no_rest_mean=mean(d_all_no_rest,2);

    % Plot data

    
    % 1. All tasks with contrasts
    figure;

    % 1.1. all versus 0

    for t=1:n_tasks
        d=+msk;
        d(msk)=task.(all_tasks{t}).edge_stats;
        subplot_tight(4,n_tasks,t,margin);
    %     draw_atlas_boundaries(d')
        summarize_matrix_by_atlas(d');
        task_magnitude_wb(t)=mean(abs(task.(all_tasks{t}).edge_stats));
        caxis(caxis_single)
    end

    [r_magn_v_scandur(1),r_magn_v_scandur(2)]=corr(task_duration',task_magnitude_wb');
    % clear correlation btw scan duration and magnitude of effect size (as expected)

    % 1.2. all versus rest

    d_contrast=+msk;
    d_contrast(msk)=task.(rest_str).edge_stats;
    % figure;
    for t=1:n_tasks
        d=+msk;
        d(msk)=task.(all_tasks{t}).edge_stats;
        if t==n_tasks
%             d_contrast(msk)=-d_all_mean;
            d=0; % plot 0-rest in the final row for comparison
        end
        subplot_tight(4,n_tasks,n_tasks+t,margin);
        d2=d-d_contrast;
    %     draw_atlas_boundaries(d2')
        summarize_matrix_by_atlas(d2');
        [r_v_rest(t,1),r_v_rest(t,2)]=corr(d2(msk),-d_contrast(msk));
        caxis(caxis_contrast)
    end

    [r_contrast_v_scandur(1),r_contrast_v_scandur(2)]=corr(task_duration(1:end-1)',r_v_rest(1:end-1,1));
    % return

    % 1.3 all versus mean task - incl rest

    d_contrast=+msk;
    d_contrast(msk)=d_all_mean;
    for t=1:length(all_tasks)
        d=+msk;
        d(msk)=task.(all_tasks{t}).edge_stats;
        subplot_tight(4,n_tasks,n_tasks*2+t,margin);
    %     draw_atlas_boundaries(d'-d_contrast')
        summarize_matrix_by_atlas(d'-d_contrast');
        caxis(caxis_contrast)
    end

    % 1.4 all versus mean task - NOT incl rest

    d_contrast=+msk;
    d_contrast(msk)=d_all_no_rest_mean;
    for t=1:length(all_tasks)
        d=+msk;
        d(msk)=task.(all_tasks{t}).edge_stats;
        subplot_tight(4,n_tasks,n_tasks*3+t,margin);
    %     draw_atlas_boundaries(d'-d_contrast')
        summarize_matrix_by_atlas(d'-d_contrast');
        caxis(caxis_contrast)
    end

    colormap(bipolar([],0.1));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1])
%     if save_figs; print(gcf,[task_contrast_file_prefix,'task_v_0_rest_meantask_meantaskminusrest',trimmed_str],'-dpng','-r300'); end
%     
%     fprintf('Corr between scan duration and magnitude of effect size: r=%0.4f, p=%0.4f\n Corr between scan duration and bias towards rest: r=%0.4f, p=%0.4f\n', ...
%     r_magn_v_scandur(1),r_magn_v_scandur(2),r_contrast_v_scandur(1),r_contrast_v_scandur(2));
    

    % mean task (not incl rest) minus rest - using the previously defined
    % d_contrast (mean not incl rest) and d (rest)
    figure;
    summarize_matrix_by_atlas(d_contrast'-d');
    caxis(caxis_contrast)
    colormap(bipolar([],0.1));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5])
    if save_figs; print(gcf,[task_contrast_file_prefix,'meantaskminusrest_v_rest',trimmed_str],'-dpng','-r300'); end


    % 2. Paired contrasts

    figure;
    for this_task=1:n_tasks
        d=+msk;
        d(msk)=task.(all_tasks{this_task}).edge_stats;
        for this_other_task=1:n_tasks
            d_contrast=+msk;
            d_contrast(msk)=task.(all_tasks{this_other_task}).edge_stats;
            subplot_tight(n_tasks,n_tasks,n_tasks*(this_task-1)+this_other_task,margin);
    %         draw_atlas_boundaries(d'-d_contrast')
            summarize_matrix_by_atlas(d'-d_contrast');
            caxis(caxis_contrast)
        end
    end

    colormap(bipolar([],0.1));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1])
    if save_figs; print(gcf,[task_contrast_file_prefix,'paired_contrasts',trimmed_str],'-dpng','-r300'); end
    

% Get run time (not polished but this is the idea)
% requires flattenStruct2Cell; data must be mounted
case 'get_run_time'
    
    error('This procedure is not ready for \"run time\" ;) ')
    
    all_stat_types={'Parametric_Bonferroni','Parametric_FDR','Size_Extent','TFCE','Constrained_FWER','Constrained','Omnibus_Multidimensional_cNBS'};
    all_grsizes={'40','80','120'};
    
    for g=1:length(all_grsizes)
%         fprintf(['\nGroup size ',all_grsizes{g},'\n']);
        load(['~/Documents/data/mnt/combined_figs_alias/combined_grsize',all_grsizes{g},'_05102021_summary.mat'],'log_data')
%         addpath('../../general/flattenStruct2Cell/')
        [log_data2,log_data2_names]=flattenStruct2Cell(log_data);
        for st=1:length(all_stat_types)
            names1=contains(log_data2_names,all_stat_types{st}) & contains(log_data2_names,'run_time');
            run_times(st,g)=nanmean(cell2mat(log_data2(names1)));
%             fprintf('%s: %1.4f\n',st{1},nanmean(cell2mat(log_data2(names1))));
        end
    end
    
    fprintf('%.1f, %.1f, %.1f\n',run_times')

otherwise
    errro(sprintf('Specified command ''%s'' does not exist.',supp_summary_type))
end