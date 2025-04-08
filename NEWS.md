# mig 2.0  (Release date 2025-04-10)

## New:

- `kdens_bandwidth`: new kernel density estimators: truncated normal (`tnorm`) and Gaussian on rotated half-space (`hsgauss)
- `kdens_bandwidth`: new approximation method (`kernel`) using a plug-in from the sample
- `mle_truncgauss`: numerical optimization routine for truncated Gaussian vector
- `proj_hs`: projection matrix to and from the half-space

## Changes:

- `approx="tnorm"` for `kdens_bandwidth` now uses the maximum likelihood estimator (rather than the sample mean and variance).
- faster calculations for leave-one-out cross validation.
- Breaking change: `mig_kdens_bandwidth` renamed to `kdens_bandwidth`, arguments names also changed


# mig 1.0 (Release date: 2024-07-14)

Initial release
