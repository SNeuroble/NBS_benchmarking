function test_stat=NBSglm_smn(GLM)
%NBSglm  Fits a general linear model (GLM).
%
%   Test_Stat=NBSglm(GLM) operates on each GLM defined by the structure
%   GLM. 
%
%   % OLD - now no uibar update - Test_Stat=NBSglm(GLM,H) attempts to write out progress to uiwaitbar 
%   with handle H.
%
%   A GLM structure contains the following fields:
%       % OLD - new and previous version of NBS does not seem to use the loop over perms in this function - GLM.perms:        Number of permutations to generate
%       GLM.X:            n x p design matrix, including a column of ones
%                         if necessary. p is the number of independent
%                         variables, n is the number of observations. 
%       GLM.y:            n x M array, where each column stores the 
%                         dependent variables for a seperate GLM
%       GLM.contrast:     1 x p contrast vector, only one contrast can be
%                         specified
%       GLM.test:         Type of test: {'onesample','ttest','ftest'}
%       GLM.exchange:     n x 1 vector specifying valid exchange blocks for 
%                         a repeated measures design. An exchange block of
%                         [1 2 3 1 2 3] confines permutation to
%                         observations 1 & 4, 2 & 5 and 3 & 6 [optional]
%       TODO: add stuff from NBSglm_setup
%
%   See also GLMexamples
%
%   Test_Stat is a 1 x M array of test-statistics representing the resulting test statistic.
%
%   Remarks: 
%       A column of ones is always included in the F-statistic reduced 
%       model, unless doing so results in a rank deficient design matrix.
%       However, a column of ones is not included in the full model unless
%       specified by the user in GLM.X. 
%
%       Any independent variable with a zero contrast is treated as a 
%       nuisance regressor. The method of Freedman & Lane is used to deal 
%       with nuisance regressors. This method is described in Anderson & 
%       Robinson (2001) Permutation tests for linear models. 43(1):75-88
%
%       The same permutation sequence is applied to every GLM. This is
%       important if spatial dependencies between each GLM are used in 
%       subsequent cluster-based inference. 
%       
%   azalesky@unimelb.edu.au


test_stat=zeros(1,GLM.n_GLMs);

beta=zeros(GLM.n_predictors,GLM.n_GLMs);
beta=GLM.X\GLM.y;

%Compute statistic of interest
% TODO: consider moving to switch/case in glm_setup
if strcmp(GLM.test,'onesample')
        
        % TODO:  why isn't a standardized effect being calculated? this would influence the NBS extent/intensity calculation bc the initial threshold is a std effect size 
        test_stat(:)=mean(GLM.y); 

elseif strcmp(GLM.test,'ttest')
    resid=zeros(GLM.n_observations,GLM.n_GLMs);
    mse=zeros(GLM.n_observations,GLM.n_GLMs);
    resid=GLM.y-GLM.X*beta;
    mse=sum(resid.^2)/(GLM.n_observations-GLM.n_predictors);
    se=sqrt(mse*(GLM.contrast*inv(GLM.X'*GLM.X)*GLM.contrast'));
    test_stat(:)=(GLM.contrast*beta)./se;
elseif strcmp(GLM.test,'ftest')
    sse=zeros(1,GLM.n_GLMs);
    ssr=zeros(1,GLM.n_GLMs);
    %Sum of squares due to error
    sse=sum((GLM.y-GLM.X*beta).^2);
    %Sum of square due to regression
    ssr=sum((GLM.X*beta-repmat(mean(GLM.y),GLM.n_observations,1)).^2);
    if isempty(GLM.ind_nuisance)
        test_stat(:)=(ssr/(GLM.n_predictors-1))./(sse/(GLM.n_observations-GLM.n_predictors));
    else
        %Get reduced model
        %Column of ones will be added to the reduced model unless the
        %resulting matrix is rank deficient
        X_new=[ones(GLM.n_observations,1),GLM.X(:,GLM.ind_nuisance)];
        %+1 because a column of 1's will be added to the reduced model
        b_red=zeros(length(GLM.ind_nuisance)+1,GLM.n_GLMs);
        %Number of remaining variables
        v=length(find(GLM.contrast))-1;
        [GLM.n_observations,ncolx]=size(X_new);
        [Q,R,perm]=qr(X_new,0);
        rankx = sum(abs(diag(R)) > abs(R(1))*max(GLM.n_observations,ncolx)*eps(class(R)));
        if rankx < ncolx
            %Rank deficient, remove column of ones
            X_new=GLM.X(:,GLM.ind_nuisance);
            b_red=zeros(length(GLM.ind_nuisance),GLM.n_GLMs);
            v=length(find(GLM.contrast));
        end

        sse_red=zeros(1,GLM.n_GLMs);
        ssr_red=zeros(1,GLM.n_GLMs);
        b_red=X_new\GLM.y;
        sse_red=sum((GLM.y-X_new*b_red).^2);
        ssr_red=sum((X_new*b_red-repmat(mean(GLM.y),GLM.n_observations,1)).^2);
        test_stat(:)=((ssr-ssr_red)/v)./(sse/(GLM.n_observations-GLM.n_predictors));
    end
end

%Added to v1.1.2
%Covers the case where the dependent variable is identically zero for all
%observations. The test statistic in this case in NaN. Therefore, force any
%NaN elements to zero. This case is typical of connectivity matrices
%populated using streamline counts, in which case some regional pairs are
%not interconnected by any streamlines for all subjects. 
test_stat(isnan(test_stat))=0; 
    
