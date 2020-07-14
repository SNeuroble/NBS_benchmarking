# New NBS statistics (principally cNBS) and benchmarking scripts

Purpose: 1. Perform command-line inference with new NBS statistics (TFCE, cNBS); 2. Evaluate performance across statistics

## Getting Started

### Prerequisites

Matlab
NBS toolbox: https://sites.google.com/site/bctnet/comparison/nbs

### Usage

#### 1. Minimal command line usage (available on devel branch)

1. Set files and parameters in setparams.m
2. Run run_NBS_cl.m

#### 2. Benchmarking

1. Set parameters in setparams_bench.m
2. Optional: If want system-dependent paths, set script and data paths for each system in setpaths.m. Must also define system_dependent_paths=1 in setparams_bench.m to use.
3. Run run_benchmarking.m

### References

- Zalesky, A., Fornito, A. and Bullmore, E.T., 2010. Network-based statistic: identifying differences in brain networks. Neuroimage, 53(4), pp.1197-1207.

- Smith, S.M. and Nichols, T.E., 2009. Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), pp.83-98.

- Noble, S., and Scheinost, D., (Accepted 2020). The constrained network-based statistic: a new level of inference for neuroimaging. MICCAI.


