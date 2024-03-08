# Multivariate inverse Gaussian

This **R** package consists of utilities for multivariate inverse Gaussian (MIG) models with mean $\boldsymbol{\xi}$ and scale matrix $\boldsymbol{\Omega}$ defined over the halfspace $\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^\top\boldsymbol{x} > 0\}$, including density evaluation and random number generation and kernel smoothing.


### Distributions

- `mig` for the MIG distribution(`rmig` for random number generation and `dmig` for density)
- `tellipt` (`rtellipt` for random vector generation and `dtellipt` the density) for truncated Student-$t$ or Gaussian distribution over the half space $\{\boldsymbol{x}: \boldsymbol{\beta}^\top\boldsymbol{x}>\delta\}$ for $\delta \geq 0$.
- `fit_mig` to estimate the parameters of the MIG distribution via maximum likelihood (`mle`) or the method of moments (`mom`).

### Kernel density estimation

- `mig_kdens_bandwidth` to estimate the bandwidth matrix minimizing the asymptotic mean integrated squared error (AMISE) or the leave-one-out likelihood cross validation, minimizing the Kullback--Leibler divergence. The `amise` estimators are estimated by drawing from a `mig` or truncated Gaussian vector via Monte Carlo
- `normalrule_bandwidth` for the normal rule of Scott for the Gaussian kernel
- `mig_kdens` for the kernel density estimator
- `tellipt_kdens` for the truncated Gaussian kernel density estimator
