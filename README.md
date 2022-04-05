# blipr: An R implementation of BLiP

## Introduction

In many applications, we can tell that a signal of interest exists but cannot perfectly "localize" it. For example, when regressing an outcome *Y* on highly correlated covariates (X<sub>1</sub>, X<sub>2</sub>), the data may suggest that *at least* one of (X<sub>1</sub>, X<sub>2</sub>) influences *Y*, but it may be challenging to tell which of (X<sub>1</sub>, X<sub>2</sub>) is important. Likewise, in genetic fine-mapping, biologists may have high confidence that a gene influences a disease without knowing precisely which genetic variants cause the disease. Similar problems arise in many settings with spatial or temporal structure, including change-point detection and astronomical point-source detection.

*Bayesian Linear Programming* (BLiP) is a method which jointly detects as many signals as possible while localizing them as precisely as possible. BLiP can wrap on top of nearly any Bayesian model or algorithm, and it will return a set of regions which each contain at least one signal with high confidence. For example, in regression problems, BLiP might return the region (X<sub>1</sub>, X<sub>2</sub>), which suggests that at least one of (X<sub>1</sub>, X<sub>2</sub>) is an important variable. BLiP controls the false discovery rate while also making these regions as narrow as possible, meaning that (roughly speaking) it will perfectly localize signals whenever this is possible! 

``blipr`` is an efficient python implementation of BLiP, which is designed so that BLiP can wrap on top of the Bayesian model in only one or two lines of code.

## Installation

You can install ``blipr`` from GitHub:

```R
# install.packages("remotes")
remotes::install_github("amspector100/blipr")
```

## Minimal example

Here, we apply BLiP to perform variable selection in a sparse linear regression. The first step is to generate synthetic data and fit a Bayesian model using, e.g., the [NPrior package](https://github.com/rabbitinasubmarine/NPrior).

```R
# Generate sparse linear regression data
set.seed(123); n <- 100; p <- 200
data <- blipr::generate_regression_data(n=n, p=p)

# Fit a Bayesian spike-and-slab model
nprior <- NPrior::NPrior_run(
   X=data$X, y=data$y, N=5000, prior='SpSL-G'
)
post_samples <- t(nprior$ThetaSamples)

# Run BLiP on the posterior samples
detections <- blipr::BLiP(
  samples=post_samples, q=0.1, error='fdr'
)
for (j in 1:length(detections)) {
  group <- paste(detections[[j]]$group, collapse=', ')
  cat("BLiP has detected a signal in ", group, ".\n", sep='')
}
```

## Documentation

Documentation and tutorials are available at [amspector100.github.io/blipr](https://amspector100.github.io/blipr).

## Development notes

Currently, there are utility two functions (``lattice_peps`` and ``hierarchical_groups``) which are slower than the version in the python pacakge, [pyblip](https://github.com/amspector100/pyblip). Of course, these are just helper functions and blipr's core functionality is quite fast. Nonetheless, we hope to substantially speed up the implementation of these functions in the near future, which will speed up the ``BLiP_cts`` wrapper.

## Reference

If you use ``blipr`` or BLiP in an academic publication, please consider citing Spector and Janson (2022). The bibtex entry is below:

```
@article{AS-LJ:2022,
  title={Controlled Discovery and Localization of Signals via Bayesian Linear Programming},
  author={Spector, Asher and Janson, Lucas},
  journal={arXiv preprint arXiv:2203.17208},
  url={https://arxiv.org/abs/2203.17208},
  year={2022}
}
```
