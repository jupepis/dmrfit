## Update

This update follows version 0.1.0 closely because it fixes a bug that makes constrained fits fail.

This is a bug-fix update of version 0.1.0. Constrained fits, `dmrfit(data, structure = )`, could fail and return NULL
for some network structures, because the trust-region step produced NaN values. The step is now computed on the free
parameters only. Regression tests were added. Unconstrained fits give the same results as in version 0.1.0.

The default initial step size of the adaptive sampler in `dmrfit_bayes()` is also changed, from 0.1 to 1.0, the value
used in the accompanying paper.

The Savage-Dickey Bayes factors of `dmrfit()` and `dmrfit_bayes()` are now more stable across seeds: when few
posterior draws lie near zero, the density at zero is estimated from a normal approximation of the draws, instead of a
kernel estimate from those few draws or a floor.

## Test environments

* local: macOS 26.6 (aarch64), R 4.5.2
* GitHub Actions: Windows (R release), Ubuntu (R devel, release, oldrel-1)
* win-builder: R devel
* mac-builder: R release

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: the update fixes a bug in constrained fits (see above).

## Reverse dependencies

There are no reverse dependencies.
