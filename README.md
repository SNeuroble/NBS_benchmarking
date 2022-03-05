# Inferential Procedures for Networks & Empirical Benchmarking

**NOTE:** The main repository has been updated to reflect recent benchmarking and summarization procedures as of 01/01/2022. For previous version, see branch "old-master".

## Purpose:
1. Perform inference in networks at various scales and from the Matlab command line
2. Empirically benchmark and compare performance of inferential procedures

### Inferential procedures currently include:
- edge-level: (FWER, parametric)
- edge-level: (FDR, parametric)
- edge-level: (FDR)
- component/cluster-level: Network-Based Statistic (NBS; FWER; Options: Size or Intensity; Zalesky, Fornito, & Bullmore, 2010)
- component/cluster-level: Threshold-Free Cluster Enhancement (TFCE; FWER; Smith & Nichols, 2009)
- network-level: Constrained NBS (cNBS; FWER)
- network-level: Constrained NBS (FDR)
- whole brain-level/omnibus: Options: Threshold_Positive, Threshold_Both_Dir, Average_Positive, Average_Both_Dir, Multidimensional_cNBS, Multidimensional_all_edges

(P-values obtained nonparametrically unless noted otherwise.)

Note: All procedures besides *NBS* and *edge-level (nonparametric FDR)* are implemented here (so any mistakes are mine!), relying in part on underlying functionality in the NBS toolbox (see NBS_addon for extending scripts). cNBS and multidimensional cNBS are introduced here (Noble & Scheinost, 2020). Connectome-based empirical benchmarking is introduced here (Noble & Scheinost, 2020).



## Prerequisites

[Matlab](https://www.mathworks.com/products/matlab.html)

[NBS toolbox](https://sites.google.com/site/bctnet/comparison/nbs)

## Usage

### Network-Based Inference

0. Open Matlab
1. Set paths and parameters in `setparams.m`
    - Example material for testing can be found in the NBS toolbox and NBS_benchmarking toolbox (this toolbox):
        - NBS toolbox "SchizophreniaExample" directory: example data and design matrix for schizophrenia study
        - NBS_benchmarking toolbox "NBS_addon" directory: simple and Shen edge groups
    - See below for tips for the construction of design matrices, contrasts, exchangeability, and two-sided tests
2. Add main NBS_benchmarking folder and subfolders (e.g., `addpath(genpath('~/NBS_benchmarking'))`)
3. Run `run_NBS_cl.m`
4. View results are all in the nbs variable (e.g., p-values are in nbs.NBS.pval). A sample visualization of the results is provided for cNBS.


### Empirical Benchmarking of Accuracy Metrics

1. Set paths and parameters
    - Set script and data paths in `setparams_bench.m`
        - Optional: If want system-dependent paths, set paths for each system in `setpaths.m`. Must set system_dependent_paths=1 in `setparams_bench.m` to use. This will overwrite paths in `setparams_bench.m`, so no need to set paths in `setparams_bench.m`.
    - Set parameters and script/data paths in `setparams_bench.m` (e.g., do_TPR, use_both_tasks, etc.)
2. Run resampling procedure
    - Run `run_benchmarking.m`
3. Calculate ground truth
    - Set `task_gt` in `setparams_bench.m`
    - Run `calculate_ground_truth.m`
3. Summarize accuracy & other results
    - Set parameters for resampling results to be summarized in `setparams_summary.m`
    - If doing summary from another workstation, mount these directories and re-define paths for resampling results and ground truth data paths. This is where `system_dependent_paths` will come in handy (see Step 1.)
    - Set date/time info for resampling results to be summarized in `set_datetimestr_and_files.m`
    - Run `summarize_tprs.m` or `summarize_fprs.m`


### Tips for constructing design matrices, contrasts, exchangeability, and two-sided tests

- Example 1. Two-sample test for 6 subjects split into 2 groups (S1G1 S2G1 S3G1 S4G2 S5G2S6G2)
    - design matrix: [1 0; 1 0; 1 0; 0 1; 0 1; 0 1];
    - contrasts: [1,-1]
    - exchangeability: None
    - Note: for now, parametric edge-level inference only performs t->p estimation for paired, not two-sample, t-test (you will receive a warning to this effect)
- Example 2. Paired-sample test for 2 subjects with 2 measurements each (S1M1 S2M1 S1M2 S2M2) 
    - design matrix: [1 1 0; 1 0 1; -1 1 0; -1 0 1];
    - contrasts: [1, -1]
    - exchangeability: [1, 2, 1, 2]; - required for proper permutation for paired-sample
- Example 3. Correlation for 4 subjects with a continuous measure (S1 S2 S3 S4)
    - design matrix: [1 5; 1 0.3; 1 4; 1 3.4]
    - contrasts: [0, 1]
    - exchangeability: None
- All functions are designed to perform one-sided tests. To perform a two-sided test, set the alpha parameters to your desired alpha value divided by two, run tests for the original contrast and the opposite, and combine results.
    - For example, for desired alpha=0.05, set alpha=0.025 and run two contrasts: [1, -1] and [-1, 1]).
- For reference, some excellent guides for constructing models can be found here:
    - [FSL GLM page](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM)
    - [Freesurfer GLM tutorial](http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/centos6/fsl_507/doc/wiki/attachments/GLM/JMglm.pdf)

(H/t Raimundo Rodriguez for pointing out areas to clarify motivating this section.)

## References

- Zalesky, A., Fornito, A. and Bullmore, E.T., 2010. Network-based statistic: identifying differences in brain networks. Neuroimage, 53(4), pp.1197-1207.

- Smith, S.M. and Nichols, T.E., 2009. Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), pp.83-98.

- Noble, S. and Scheinost, D., 2020. The Constrained Network-Based Statistic: A New Level of Inference for Neuroimaging. In International Conference on Medical Image Computing and Computer-Assisted Intervention (pp. 458-468). Springer, Cham.

