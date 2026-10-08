## Update

This is a bug-fix update of version 0.1.0. Constrained fits, `dmrfit(data, structure = )`, could fail and return NULL
for some network structures, because the trust-region step produced NaN values. The step is now computed on the free
parameters only. Regression tests were added. Unconstrained fits give the same results as in version 0.1.0.

## Test environments

* local: macOS 26.6 (aarch64), R 4.5.2
* GitHub Actions: Windows (R release), Ubuntu (R devel, release, oldrel-1)
* win-builder: R devel
* mac-builder: R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.
