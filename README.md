# Expanded NBS tools for new statistics and benchmarking

Purpose: 1. Perform command-line inference with new NBS statistics (cNBS, TFCE); 2. Evaluate performance across statistics

## Getting Started

### Prerequisites

Matlab
NBS toolbox: https://sites.google.com/site/bctnet/comparison/nbs

### Usage

#### 1. Minimal command line usage (available on devel branch)

1. Set files and parameters in setparams.m
2. Run run_NBS_cl.m
3. Results are all in the nbs variable (e.g., p-values are in nbs.NBS.pval). A sample visualization of the results is provided for cNBS.

#### 2. Benchmarking

1. Set parameters in setparams_bench.m (e.g., do_TPR, use_both_tasks, etc.)
2. Optional: If want system-dependent paths, set script and data paths for each system in setpaths.m. Must also define system_dependent_paths=1 in setparams_bench.m to use.
3. Run run_benchmarking.m
4. Set params for results to be summarized in setparams_summary.m
4. Run summarize_tprs.m or summarize_fprs.m

### References

- Zalesky, A., Fornito, A. and Bullmore, E.T., 2010. Network-based statistic: identifying differences in brain networks. Neuroimage, 53(4), pp.1197-1207.

- Smith, S.M. and Nichols, T.E., 2009. Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), pp.83-98.

- Noble, S., and Scheinost, D., 2020. The constrained network-based statistic: a new level of inference for neuroimaging. MICCAI. (Full reference forthcoming)


