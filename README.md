# New NBS statistics (principally cNBS) and benchmarking scripts

Purpose: 1. Perform command-line inference with new NBS statistics (TFCE, cNBS); 2. Evaluate performance across statistics

## Getting Started

### Prerequisites

Matlab
NBS toolbox: https://sites.google.com/site/bctnet/comparison/nbs

### Usage

#### 1. Minimal command line usage (available on devel branch)

1. Set files and parameters in setparams_cl.m
2. Run run_NBS_cl.m

#### 2. Benchmarking

1. Set script and data paths in setpaths.m (this file is intended to enable different path presets for multiple systems)
2. Set parameters in setparams.m
3. Run run_benchmarking.m


For more about the Constrained Network Based Statistic (cNBS), please see

"The constrained network-based statistic: a new level of inference for neuroimaging" (MICCAI 2020, reference forthcoming).
