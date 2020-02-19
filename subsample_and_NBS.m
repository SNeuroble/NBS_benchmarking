function [edge_stats,edge_stats_neg,cluster_stats,cluster_stats_neg,...
pvals,pvals_neg,sig_result,sig_result_neg] = subsample_and_NBS(...
this_repetition,m_test,UI,n_subs,n_subs_subset,do_TPR,nbs_contrast_neg)

% subsamples data and runs NBS

    fprintf('\n* Repetition %d - positive contrast\n',this_repetition)

    % shuffle data
    %global m;
%	ids=randperm(n_subs,n_subs_subset);
%    if do_TPR
%        ids=[ids;ids+n_subs];
%    end
%
%    m_test=m(:,:,ids);

    % simulate effects
%    if rep_params.do_simulated_effect
%        effect=+ismember(edge_groups,rep_params.networks_with_effects);
%        effect=effect+effect';
%        m_with_effect=m_test;
%        m_with_effect(:,:,1:1:n_subs/2)=m_with_effect(:,:,1:n_subs/2)+effect;
%        m_test=m_with_effect;
%    end

    UI_new=UI;
    UI_new.matrices.ui=m_test;
    nbs=NBSrun_smn(UI_new);
    
    % re-run with negative
    fprintf('\n* Repetition %d - negative contrast\n',this_repetition)
    UI_new_neg=UI_new;
    UI_new_neg.contrast.ui=nbs_contrast_neg;
    nbs_neg=NBSrun_smn(UI_new_neg);
    

    % check for any positives (if there was no ground truth effect, this goes into the FWER calculation)
    if nbs.NBS.n>0; sig_result=1;
	else; sig_result=0;
	end
    if nbs_neg.NBS.n>0; sig_result_neg=1;
	else; sig_result_neg=0;
    end

    % record everything
    edge_stats=nbs.NBS.edge_stats;
    cluster_stats=full(nbs.NBS.cluster_stats);
    pvals=nbs.NBS.pval(:); % TODO: had to vectorize for TFCE... should give all outputs in same format tho 

	edge_stats_neg=nbs_neg.NBS.edge_stats;
    cluster_stats_neg=full(nbs_neg.NBS.cluster_stats);
    pvals_neg=nbs_neg.NBS.pval(:); % TODO: same as above

end
