function GLM=NBSglm_setup_smn(GLM)
%NBSglm_setup_smn  Sets up for running the general linear model (GLM) in
%NBSglm, include performing nuisance regression
%
%   GLM=NBSglm_setup_smn(GLM) updates the structure GLM.
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
%
%   This function adds the following count variables and nuisance specific
%   variables:
%         GLM.n_predictors
%         GLM.n_GLMs
%         GLM.n_observations
%         GLM.ind_nuisance
%         GLM.blks
%         GLM.sz_blk
%         GLM.blk_ind
%         GLM.b_nuisance
%         GLM.resid_y
%
%

%Number of predictors (including intercept)
GLM.n_predictors=length(GLM.contrast);

%Number of independent GLM's to fit
GLM.n_GLMs=size(GLM.y,2);

%Number of observations
GLM.n_observations=size(GLM.y,1);

%Determine nuisance predictors not in contrast
GLM.ind_nuisance=find(~GLM.contrast);

if isfield(GLM,'exchange')
    %Set up exchange blocks
    GLM.blks=unique(GLM.exchange); 
    %Number of blocks
    GLM.n_blks=length(GLM.blks);
    %Number of observations per block
    GLM.sz_blk=GLM.n_observations/GLM.n_blks; 
    GLM.blk_ind=zeros(GLM.n_blks,GLM.sz_blk);
    for i=1:length(GLM.blks)
        GLM.blk_ind(i,:)=find(GLM.blks(i)==GLM.exchange);
    end
end

if isempty(GLM.ind_nuisance)
    %No nuisance predictors
else
    %Regress out nuisance predictors and compute residual
    GLM.b_nuisance=zeros(length(GLM.ind_nuisance),GLM.n_GLMs);
    GLM.resid_y=zeros(GLM.n_observations,GLM.n_GLMs); 
    GLM.b_nuisance=GLM.X(:,GLM.ind_nuisance)\GLM.y;
    GLM.resid_y=GLM.y-GLM.X(:,GLM.ind_nuisance)*GLM.b_nuisance; 
end



