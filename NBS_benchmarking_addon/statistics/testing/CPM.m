% CPM-based inference

% Setup
resampling_param.n_reps=1000;
resampling_param.aggregation_thresh=0.75; % proportion of times edges included in model to be included in aggregation

%% Estimate model and performance

% 1. Set up for resampling
y_orig=GLM.y;
X_orig=GLM.X;
% some vector for predictions and models

% 2. Do CPM
% Note: unlike our typical CPM, this is set up so behavior = categorical predictor X, not continuous outcome Y
% Note: unlike the NBS corrections, here CPM is about resampling subset of data rather than permuting all data

for i=1:K % SMN - new range
    
    GLM=resample_data(GLM,y_orig,X_orig,resampling_param); % SMN - new function based on permute_signal, new input, new output (new y, new X, etc)
    edge_stats=get_univariate_stats(STATS,GLM,precomputed,i);
    
    % apply threshold and calculate summary stat - new function, new output
    [summary_stat(i),model{i}]=get_summary_stat(edge_stats, STATS, ind_upper, N, Intensity,bgl);
    
    % fit model
    X_predicted=fit_lm(summary_stat, GLM.X);
    
    % only print every hundred
    if mod((i-1),100)==0
        str=sprintf('| %5d/%5d perms complete | latest element in null: %6.1f |',i-1,K,null_dist(i-1,end));
        % display no more than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow;
        catch;  fprintf([str,'\n']); end
    end
    
end


%% Perform cross-validation and output final model (correction analogue)

% This gets cluster-based statistics satisfying FWER threshold
performance=estimate_performance(X_predicted,X_orig); % note that X_orig is design matrix and X_predicted is vector - convert X_orig to vec

% Aggregate model
% Can we test if the model is validated, and only return the aggregate model if yes?
% -> What is the threshold for saying the model has been validated?? Even a cross-validated model can have a lucky day.
validated=performance>0;
any_significant=validated;
if validated
   positives=aggregate_model(models,resampling_param.aggregation_thresh);
end
con_mat=0;

% Multiple comparison correction allows you to make inferences over a set 
% of tests while limiting a false positive rate to a desired level. For a
% set of features, this means selecting those features which, in isolation,
% are most likely to be associated w an outcome.
% Many machine learning approaches are not so much corrections as ways of
% allowing you to select/combine features to maximize generalizability.
% The goals of minimizing type I error at the feature level and maximizing 
% generalizability are related, but not directly. At least one example is a
% case where only a single FP is allowed in every iteration but 1000 large 
% effect size TPs are allowed, st FWER is 100% but generalizability is also
% almost 100%.

% Since the goal of machine learning is generalizable behavior and not
% specified control of a feature-level error rate (i.e., 5% of positive 
% edges allowed to be FP or 5% of iterations allowed to have at least one 
% FP edge), it's difficult to define the expected behavior.
% We might be able to kind of estimate it for CPM because we aggregate 
% results across folds to a set of features - but it seems important that 
% the aggregation step itself isn't validated or theoretically motivated.
% Other ML techniques that do not result in a selected set of features 
% (e.g., neural nets) are probably more comparable to other omnibus tests. 
% I think having a 
% well-defined aggregation step is probably important here.




