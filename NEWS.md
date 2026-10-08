# dmrfit 0.1.1

* Constrained fits, `dmrfit(data, structure = )`, no longer fail and return NULL ("Optimization failed") for some
  network structures. The trust-region step is now computed on the free parameters only (the thresholds and the
  interactions present in the structure). Unconstrained fits are unchanged.
* `dmrfit_bayes()`: the default initial step size of the adaptive FisherMALA sampler is now `sigma2 = 1.0` (was 0.1),
  the value used in Arena and Marsman (2026). The step size adapts during the burn-in, so the target distribution is
  unchanged, but the draws for a given seed differ from version 0.1.0.

# dmrfit 0.1.0

* First CRAN release.
