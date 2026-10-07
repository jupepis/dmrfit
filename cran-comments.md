## Test environments

* local: macOS 26.6 (aarch64), R 4.5.2
* GitHub Actions: Windows (R release), Ubuntu (R devel, release, oldrel-1)
* win-builder: R devel
* R-hub: linux, windows, macos-arm64, clang-asan, gcc14 and gcc16 (Fedora), all R devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

  Possibly misspelled words in DESCRIPTION: GHW, Godambe, Ising, Langevin, Liang, Marsman, Skare, Titsias, et al. These
  are names of authors and of methods (GHW is the Godambe-Huber-White covariance), and are spelled correctly.

## Comments

* The vignette is precompiled: `vignettes/dmrfit.Rmd` is generated from `vignettes/dmrfit.Rmd.orig` (excluded from the
  package build), because its model fits take about 30 seconds.
* The dataset `rads2` is distributed under the CC BY 4.0 license, as stated in the `Copyright` field of DESCRIPTION,
  in `?rads2` and in `inst/COPYRIGHTS`. The rest of the package is under the MIT license.
* Parallel computation through OpenMP is controlled by the argument `ncores`, with default 1. Examples and tests use at
  most 2 cores.
* `igraph` (in Suggests) provides the default network layout of `plot()`. Examples and tests run without it.
