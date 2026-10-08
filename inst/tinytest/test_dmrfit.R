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
expect_equal(names(fit$argument)[c(1, 13, 18)], c("mu[1,1]", "theta[2,1]", "theta[4,3]"))
expect_true(all(is.finite(fit$argument)))
# data coded from 1 (as in rads2) and from 0 give the same fit (the baseline category is the minimum)
expect_equal(as.vector(dmrfit(as.matrix(X) - 1L, with_prior = TRUE)$argument), as.vector(fit$argument))

# --- parallel computation (OpenMP): two cores give the same estimates, Hessian and GHW covariance as one
fit_1 <- dmrfit(X, with_prior = TRUE, ncores = 1)
fit_2 <- dmrfit(X, with_prior = TRUE, ncores = 2)
expect_identical(fit_2$argument, fit_1$argument)
expect_identical(fit_2$utils, fit_1$utils)

# --- constrained estimation: absent edges stay at zero
S <- matrix(1, 4, 4); diag(S) <- 0; S[1, 2] <- S[2, 1] <- 0
fit_s <- dmrfit(X, structure = S)
expect_equal(unname(fit_s$argument["theta[2,1]"]), 0)
# Savage-Dickey Bayes factors only for the free interactions
fit_s_sd <- dmrfit(X, structure = S, with_prior = TRUE, savage_dickey = TRUE, M = 500)
expect_equal(names(fit_s_sd$savage_dickey$bf_01), setdiff(names(fit_s$argument)[-(1:12)], "theta[2,1]"))
expect_stdout(print(summary(fit_s_sd)), "constrained")

# --- constrained fit: the GHW covariance comes from the free parameters only, zero for the absent edges
free <- c(rep(TRUE, 12), S[lower.tri(S)] == 1)
Hf_inv <- solve(fit_s$utils$hessian[free, free])
expect_equal(fit_s$utils$HW[free, free], Hf_inv %*% fit_s$utils$FisherInfo[free, free] %*% Hf_inv)
expect_true(all(fit_s$utils$HW[!free, ] == 0) && all(fit_s$utils$HW[, !free] == 0))
# with every edge present, the constrained construction gives the GHW covariance of the unconstrained fit (with and without
# the prior)
S_all <- matrix(1, 4, 4); diag(S_all) <- 0
for (pr in c(FALSE, TRUE)) {
    expect_equal(dmrfit(X, structure = S_all, with_prior = pr)$utils$HW, dmrfit(X, with_prior = pr)$utils$HW, tolerance = 1e-8)
}
# with every edge present, the constrained trust-region steps are those of the unconstrained fit
expect_identical(dmrfit(X, structure = S_all)$argument, dmrfit(X)$argument)

# --- constrained fits whose start gave a NaN trust-region step (exact zero eigenvalues of the absent edges in the masked
# Hessian), and returned NULL
for (s in c(1, 22, 30)) {
    set.seed(s); P <- sample(4:8, 1)
    X_s <- matrix(rbinom(300 * P, 1, 0.4), 300, P)
    A <- matrix(rbinom(P * P, 1, 0.4), P, P); A[upper.tri(A, TRUE)] <- 0; A <- A + t(A)
    fit_ns <- dmrfit(X_s, structure = A)
    expect_false(is.null(fit_ns))
    # the gradient vanishes on the free parameters, and the absent edges stay at zero
    free_ns <- c(rep(TRUE, P), A[lower.tri(A)] == 1)
    expect_true(max(abs(fit_ns$utils$gradient[free_ns])) < 1e-4)
    expect_true(all(fit_ns$argument[!free_ns] == 0))
}
set.seed(1)
X_s <- matrix(rbinom(800, 1, 0.4), 100, 8)
A <- matrix(0, 8, 8); A[2, 1] <- A[1, 2] <- 1
expect_false(is.null(dmrfit(X_s, structure = A)))

# --- profile likelihood-ratio intervals: the adjusted profile statistic equals the chi-square quantile at both bounds
fit_lrt <- dmrfit(X, with_prior = TRUE, lrt_intervals = TRUE, level = 0.9)
li <- fit_lrt$lrt_intervals
expect_equal(rownames(li), names(fit_lrt$argument))
expect_equal(li$C, diag(solve(fit_lrt$utils$hessian)) / diag(fit_lrt$utils$HW))
expect_true(all(li$lower < li$estimate & li$estimate < li$upper))
Xs <- cbind(as.matrix(X) - 1, 2 * t(apply(as.matrix(X) - 1, 1, function(x) { M <- x %*% t(x); M[lower.tri(M)] })))
profile_nll <- function(k, value) {
    dmrfit:::cpp_optimize_profile(data = Xs, parinit = fit_lrt$argument, which_parconstr = k - 1, parconstr = value,
                                  n_categories = rep(4, 4), P = 4, f_term = sqrt(.Machine$double.eps),
                                  m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1, rmax = 10,
                                  with_prior = TRUE, epsilon = 1e-6, ncores = 1L, thresholds_alpha = 0.5,
                                  thresholds_beta = 0.5, interactions_location = 0, interactions_scale = 2.5)$utils$value
}
stat <- function(k, value) li$C[k] * 2 * (profile_nll(k, value) - fit_lrt$utils$value)
for (k in c(1, 12, 13, 18)) {
    expect_equal(c(stat(k, li$lower[k]), stat(k, li$upper[k])), rep(qchisq(0.9, 1), 2), tolerance = 1e-3)
}
# the profile fixes theta_k and optimizes the other parameters: their gradient is close to zero
prof <- dmrfit:::cpp_optimize_profile(data = Xs, parinit = fit_lrt$argument, which_parconstr = 12, parconstr = li$upper[13],
                                      n_categories = rep(4, 4), P = 4, f_term = sqrt(.Machine$double.eps),
                                      m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1, rmax = 10,
                                      with_prior = TRUE, epsilon = 1e-6, ncores = 1L, thresholds_alpha = 0.5,
                                      thresholds_beta = 0.5, interactions_location = 0, interactions_scale = 2.5)
expect_equal(prof$argument[13], li$upper[13])
expect_true(max(abs(prof$utils$gradient[-13])) < 1)
# a subset of parameters by name
expect_equal(rownames(dmrfit(X, lrt_intervals = c("theta[2,1]", "mu[1,1]"))$lrt_intervals), c("theta[2,1]", "mu[1,1]"))
expect_error(dmrfit(X, lrt_intervals = "theta[9,1]"), "unknown parameter")
expect_true(is.null(dmrfit(X)$lrt_intervals))
# fewer than 10 observations per free parameter: a warning suggests the likelihood-ratio intervals
expect_warning(dmrfit(X[1:100, ]), "lrt_intervals = TRUE")
expect_silent(dmrfit(X[1:100, ], lrt_intervals = "theta[2,1]"))
expect_silent(dmrfit(X))
# constrained fit: intervals for the free parameters only
expect_equal(rownames(dmrfit(X, structure = S, lrt_intervals = TRUE)$lrt_intervals), names(fit_s$argument)[free])
expect_error(dmrfit(X, lrt_intervals = NA), "TRUE, FALSE or a vector")
expect_error(dmrfit(X, level = 1), "between 0 and 1")

# --- confint(): stored likelihood-ratio intervals, or Wald intervals from the GHW covariance
expect_equal(unname(confint(fit_lrt, level = 0.9)), unname(as.matrix(li[, c("lower", "upper")])))
expect_equal(colnames(confint(fit_lrt, level = 0.9)), c("5 %", "95 %"))
wald <- confint(fit_lrt, parm = "theta[2,1]", method = "wald")  # at the level of the fit, 0.9
expect_equal(unname(wald[1, ]), unname(fit_lrt$argument["theta[2,1]"] + c(-1, 1) * qnorm(0.95) * sqrt(fit_lrt$utils$HW[13, 13])))
expect_equal(rownames(confint(fit_lrt, parm = c(1, 13), level = 0.9)), names(fit_lrt$argument)[c(1, 13)])
expect_equal(confint(fit_lrt), confint(fit_lrt, level = 0.9))  # default: the level of the fit
expect_error(confint(fit_lrt, level = 0.95), "computed at level 0.9")
# summary(): Wald and likelihood-ratio intervals next to each other, at the level of the fit
sm <- summary(fit_lrt)
expect_equal(sm$level, 0.9)
expect_equal(colnames(sm$intervals), c("Wald lower", "Wald upper", "LRT lower", "LRT upper"))
expect_equal(unname(sm$intervals[, 3:4]), unname(as.matrix(li[, c("lower", "upper")])))
expect_equal(unname(sm$intervals[, 1:2]), unname(confint(fit_lrt, method = "wald")))
expect_equal(colnames(summary(dmrfit(X))$intervals), c("Wald lower", "Wald upper"))
expect_stdout(print(sm), "Intervals \\(90%\\)")
expect_error(confint(dmrfit(X)), "lrt_intervals = TRUE")
expect_error(confint(fit_lrt, parm = "theta[9,1]", level = 0.9), "unknown parameter")
expect_error(confint(fit_s, parm = "theta[2,1]", method = "wald"), "absent edges")

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
