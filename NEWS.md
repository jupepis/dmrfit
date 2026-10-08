# dmrfit 0.1.1

* Constrained fits, `dmrfit(data, structure = )`, no longer fail and return NULL ("Optimization failed") for some
  network structures. The trust-region step is now computed on the free parameters only (the thresholds and the
  interactions present in the structure). Unconstrained fits are unchanged.
* `dmrfit_bayes()`: the default initial step size of the adaptive FisherMALA sampler is now `sigma2 = 1.0` (was 0.1),
  the value used in Arena and Marsman (2026). The step size adapts during the burn-in, so the target distribution is
  unchanged, but the draws for a given seed differ from version 0.1.0.
* The Savage-Dickey Bayes factors of `dmrfit()` and `dmrfit_bayes()` are now more stable across seeds: when few
  posterior draws lie near zero, the density at zero is estimated from a normal approximation of the draws, instead of
  a kernel estimate from those few draws or a floor. The estimator used for each
  interaction is in `savage_dickey$density_at_zero`, and `summary()` reports how many Bayes factors use the normal
  approximation.

# dmrfit 0.1.0

* First CRAN release.
