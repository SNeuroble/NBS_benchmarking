function compare_methods(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script summarizes and visualizes methods, using data already
% combined across tasks
%
% Usage: compare_methods('tasks',{'SOCIAL_v_REST'},'stat_types',{'Size_Extent'},'grsize',40,'make_figs',0);
%   Task choices: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION; MOTOR; GAMBLING
%
% Required: Ground truth summary (see calculate_ground_truth), benchmarking results
% 
% STEPS
% Summarization: fits spline to effect size vs. mean TPR
% Plot: d v. TPR spline, d v. TPR residual map
%
% Recommended to first run with local access to data to obtain intermediate 
% summary file--this intermediate file is slow to create but can then be
% used to recreate any summaries/visualizations. Note that this step
% doesn't rely on the ground truth
% (When using remote data, mount data dir: sshfs smn33@172.23.202.124:d3_smn33/ mnt/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARSE PARAMETERS

all_stat_types={'Size_Extent','TFCE','Constrained'}; % TODO: add FDR and Omnibus_Threshold_Both_Dir
all_stat_scalings=[1,1,2]; % 1=edges, 2=networks, 3=omni
grsize_default=40; % TODO: change to 80
save_figs_default=1; % TODO: change back to 1
save_log_default=1; % TODO: change back to 1

p = inputParser;
addOptional(p,'stat_types',all_stat_types);
addOptional(p,'stat_scalings',all_stat_scalings);
addOptional(p,'grsize',grsize_default);
addOptional(p,'save_figs',save_figs_default);
addOptional(p,'save_log',save_log_default);
parse(p,varargin{:});

stat_types=p.Results.stat_types;
stat_scalings=p.Results.stat_scalings;
grsize=p.Results.grsize;
save_figs=p.Results.save_figs;
save_log=p.Results.save_log;

%% Setup

% setup
task='all_tasks';
stat_type='all_stats';
[current_path,~,~]=fileparts(mfilename('fullpath')); % assuming current folder is NBS_benchmarking
addpath(genpath(current_path));
setpaths;
setparams_metasummary;
warning('off', 'SPLINES:CHCKXYWP:NaNs');

date_time_str=datestr(now,'mmddyyyy');
summary_output_dir=[output_dir,task,'_',stat_type,'_summary/'];
method_comparison_basename_prefix=[summary_output_dir,'results__',task,'_full_comparison','_','grsize',num2str(grsize),'_',date_time_str];
% method_comparison_basename_prefix=['results__',task,'_full_comparison','_','grsize',num2str(grsize),'_',date_time_str];
% method_comparison_filename=[output_dir,comparison_basename_prefix,'_summary.mat'];

% make summary output dir
if ~exist(summary_output_dir,'dir'); mkdir(summary_output_dir); end

% TODO: check whether to save

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET METHOD-SPECIFIC SUMMARY (COMBINED TASKS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('* Summarizing true positive benchmarking results.\n');

for s=1:length(stat_types)

    fprintf(['Loading method ',stat_types{s},'\n'])

    % TODO: define date_time_str_combined.(stat_type) in "setparams_summary"
    combined_summary_basename_prefix=['results__',task,'_',stat_types{s},'_','grsize',num2str(grsize),'_',date_time_str_combined.(stat_types{s})];
    combined_summary_filename=[output_dir,combined_summary_basename_prefix,'_summary.mat'];

    % TODO: bring this back
    load(combined_summary_filename,'dcoeff_scaled_all','tpr_scaled_all','log_data_combined');
%     this_summary{s}=load(combined_summary_filename);
    dcoeff.(stat_types{s})=dcoeff_scaled_all{1};
    tpr.(stat_types{s})=tpr_scaled_all{1};
    log.(stat_types{s})=log_data_combined;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO VISUALIZATION COMPARING METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbins=ceil((pp.ax_ymax_tp-pp.ax_ymin)/pp.tpr_bin_width);
bin_edges=linspace(pp.ax_ymin,pp.ax_ymax_tp,nbins+1);


% 0. Plot TPR histogram (cumulative histogram? how many effects have x power or greater?)

figure;
hold on;

for s=1:length(stat_types)
    tpr_histcounts=histcounts(tpr.(stat_types{s}){stat_scalings(s)}(:),bin_edges,'Normalization','probability');
    bin_centers=bin_edges(1:end-1)+pp.tpr_bin_width/2;
    ndist = tpr_histcounts / sum(tpr_histcounts);
    cdist = cumsum(ndist,'reverse'); % percent of effects showing X power or greater
    plot(bin_centers,cdist,'-','LineWidth',2,'MarkerEdgeColor',stat_color_order(s,:));

    avg_tpr(s)=mean(tpr.(stat_types{s}){stat_scalings(s)}(:));
end

legend(stat_types,'Interpreter', 'none');
set(gca,'fontsize',pp.fontsz)
hold off;

if save_figs
    saveas(gcf,[method_comparison_basename_prefix,'_TPR_revcumhist'],'png')
end

% bar plot of average power

figure;
b=bar(avg_tpr,'FaceColor','flat');
% bar(diag(avg_tpr),'stacked');
set(gca,'xticklabel',stat_types,'TickLabelInterpreter','none','fontsize',pp.fontsz);
for s = 1:length(avg_tpr) % colors
    b.CData(s,:) = stat_color_order(s,:);
end
ylim([pp.ax_ymin pp.ax_ymax_tp]);

if save_figs
    saveas(gcf,[method_comparison_basename_prefix,'_avg_TPR'],'png')
end

% 1. Plot effect size vs. TPR
% TODO: get only absolute value

for do_abs=0:1        

    figure;
    hold on;

    for s=1:length(stat_types)

        if do_abs
            this_dcoeff=abs(dcoeff.(stat_types{s}){stat_scalings(s)}(:));
            this_tpr=abs(tpr.(stat_types{s}){stat_scalings(s)}(:));
        else
            this_dcoeff=dcoeff.(stat_types{s}){stat_scalings(s)}(:);
            this_tpr=tpr.(stat_types{s}){stat_scalings(s)}(:);
        end

        [tpr_fit,res_scaled,dcoeff_windowed,tpr_windowed,tpr_std,~]=...
                fit_spline(this_dcoeff,this_tpr,pp.spline_smoothing{stat_scalings(s)},pp.window_sz{stat_scalings(s)});
        [~,ind]=sort(this_dcoeff); 
        plot(this_dcoeff(ind),tpr_fit(ind),'-','LineWidth',2,'Color',stat_color_order(s,:))

        if exist('shadedErrorBar','file')
            [~,min_val]=min(this_dcoeff);
%             shadedErrorBar(dcoeff_windowed,tpr_windowed,tpr_std,'noLine',1,'lineProps',{'-','color',stat_color_order(s,:)})
            shadedErrorBar([this_dcoeff(min_val) dcoeff_windowed],[tpr_fit(min_val) tpr_windowed],[tpr_std(1) tpr_std],'noLine',1,'lineProps',{'-','color',stat_color_order(s,:)})
        else
            warning('Can''t find shadedError function, so won''t draw shaded error bars.')
        end
    end

    if (do_abs); ax_xmin=0; else ax_xmin=pp.ax_xmin; end  
    axis([ax_xmin,pp.ax_xmax,pp.ax_ymin,pp.ax_ymax_tp]);
    
    if save_figs
        if do_abs
            abs_str='_abs';
        else
            abs_str='';
        end
        saveas(gcf,[method_comparison_basename_prefix,'_TPR_v_esz',abs_str],'png')
    end
    hold off;

end

% Log

% total % effects detected


if save_log
    logfile=[method_comparison_basename_prefix,'_log','.txt'];
    fprintf('Saving log in %s.\n',logfile);
    fid=fopen(logfile,'w');

    for s=1:length(stat_types)

        fprintf(fid,'Mean TPR =%1.3f\n',mean(tpr.(stat_types{s}){stat_scalings(s)}));
        fprintf(fid,'Run time: %1.2f hours',log_data_combined.run_time_h); % toc is in sec
%         fprintf(fid,'Mean TPR between d=+/-%1.1f: %1.3f\n',[thresholds; tpr_lt_thr]);
%         fprintf(fid,'Mean TPR between d=%1.1f and %1.1f: %1.3f\n',[thresholds(2:end); thresholds(1:end-1); tpr_btw_thr_and_thr_below]);
%         fprintf(fid,'Mean TPR at d=+/-%1.1f: %f (+), %f (-), %f (mean)\n',[thresholds; reshape(tpr_at_thr,3,length(thresholds))]);
        % TODO: should just check that the following are the same across methods
        fprintf(fid,'\n%d total repetitions',log_data_combined.n_repetitions);
        fprintf(fid,'\n%s total permutations',log_data_combined.n_perms);
        fprintf(fid,'\n%d subjects sampled out of %d total subjects',log_data_combined.n_subs_subset,round(log_data_combined.n_subs_total));
        fprintf(fid,'\n\n');
      
    end

    fclose(fid);

end




