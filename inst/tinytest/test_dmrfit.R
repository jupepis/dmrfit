# Tests for dmrfit(): input checks, the point estimates, constrained estimation and the Savage-Dickey Bayes factors
# computed by sampling importance resampling.

data(rads2, package = "dmrfit")
X <- rads2[, c("D3", "D6", "D7", "D8")]

# --- input checks
expect_error(dmrfit(X, savage_dickey = TRUE), "requires with_prior = TRUE")
expect_error(dmrfit(-as.matrix(X)), "non-negative integers")
expect_error(dmrfit(as.matrix(X) + 0.5), "non-negative integers")
expect_error(dmrfit(X, parinit = rep(0, 3)), "Length of parinit")
expect_error(dmrfit(X, structure = diag(3)), "square matrix")
S_asym <- matrix(1, 4, 4); S_asym[1, 2] <- 0
expect_error(dmrfit(X, structure = S_asym), "symmetric")
expect_error(dmrfit(X, structure = matrix(2, 4, 4)), "only 0s and 1s")

# --- point estimates: named vector, thresholds first, then the interactions in lower-triangle order
fit <- dmrfit(X, with_prior = TRUE)
expect_inherits(fit, "dmrfit")
expect_equal(length(fit$argument), 12 + 6)
expect_equal(names(fit$argument)[c(1, 13, 18)], c("mu[1,1]", "sigma[2,1]", "sigma[4,3]"))
expect_true(all(is.finite(fit$argument)))
# data coded from 1 (as in rads2) and from 0 give the same fit (the baseline category is the minimum)
expect_equal(as.vector(dmrfit(as.matrix(X) - 1L, with_prior = TRUE)$argument), as.vector(fit$argument))

# --- parallel computation (OpenMP): two cores give the same estimates, Hessian and sandwich covariance as one
fit_1 <- dmrfit(X, with_prior = TRUE, ncores = 1)
fit_2 <- dmrfit(X, with_prior = TRUE, ncores = 2)
expect_identical(fit_2$argument, fit_1$argument)
expect_identical(fit_2$utils, fit_1$utils)

# --- constrained estimation: absent edges stay at zero
S <- matrix(1, 4, 4); diag(S) <- 0; S[1, 2] <- S[2, 1] <- 0
fit_s <- dmrfit(X, structure = S)
expect_equal(unname(fit_s$argument["sigma[2,1]"]), 0)
# Savage-Dickey Bayes factors only for the free interactions
fit_s_sd <- dmrfit(X, structure = S, with_prior = TRUE, savage_dickey = TRUE, M = 500)
expect_equal(names(fit_s_sd$savage_dickey$bf_01), setdiff(names(fit_s$argument)[-(1:12)], "sigma[2,1]"))
expect_stdout(print(summary(fit_s_sd)), "constrained")

# --- Savage-Dickey Bayes factors
set.seed(10)
before <- .Random.seed
fit_sd <- dmrfit(X, with_prior = TRUE, savage_dickey = TRUE, M = 500, oversampling = 10)
# the caller's random state is left untouched
expect_identical(.Random.seed, before)
sd <- fit_sd$savage_dickey
expect_equal(length(sd$bf_01), 6L)
expect_true(all(sd$bf_01 > 0))
expect_equal(sd$pr_null, sd$bf_01 / (1 + sd$bf_01))
expect_true(is.logical(sd$zero_beyond_draws))
# Bayes factors of the flagged interactions sit at the floor 1 / (10 M)
expect_true(all(sd$bf_01[sd$zero_beyond_draws] == sd$bf_floor))
expect_equal(sd$bf_floor, 1 / (10 * 500))
# the importance sample is not degenerate
expect_true(sd$ess > 500 && sd$ess <= sd$n_proposals)
# the same seed gives the same Bayes factors
expect_identical(dmrfit(X, with_prior = TRUE, savage_dickey = TRUE, M = 500, oversampling = 10)$savage_dickey$bf_01, sd$bf_01)

# --- print and summary
expect_stdout(print(fit_sd))
expect_stdout(print(summary(fit_sd)), "Savage-Dickey")
