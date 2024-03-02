# Multivariate inverse Gaussian

This **R** package consists of utilities for multivariate inverse Gaussian (MIG) models with mean $\boldsymbol{\xi}$ and scale matrix $\boldsymbol{\Omega}$ defined over the halfspace $\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^\top\boldsymbol{x} > 0\}$, including density evaluation and random number generation and kernel smoothing.


**To-do list**

- [ ] Add tests for functions
- [ ] Implement likelihood leave-one-out cross validation in Cpp and reduce calculations by exploiting the symmetry.
- [ ] Add robust likelihood LOO-CV
- [ ] Add support for heterogeneous compound symmetry model
- [ ] Code simulation study to fit the MIG and models to different mixtures of elliptical distributions (4 models) in different dimension ($d \in \{2,3,5\}$) and with varying sample sizes $n \in \{250, 500, 1000\}$. Measure discrepancy to the true model using a Monte Carlo estimate of the Kullback--Leibler divergence.
- [ ] Data illustration: draw posterior samples from a generalized Pareto distribution for exceedances with few exceedances. Fit kernel model to these and compare to (truncated) Gaussian approximation to the posterior with density plots.
