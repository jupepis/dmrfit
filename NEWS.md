# dmrfit 0.1.1

* Constrained fits, `dmrfit(data, structure = )`, no longer fail and return NULL ("Optimization failed") for some
  network structures. The trust-region step is now computed on the free parameters only (the thresholds and the
  interactions present in the structure). Unconstrained fits are unchanged.

# dmrfit 0.1.0

* First CRAN release.
