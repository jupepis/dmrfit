# Tests for dmrfit_bayes(): argument checks, every sampler and rescaling covariance, the stored output, the
# Savage-Dickey Bayes factors and the effective sample sizes. Short chains on three items keep the tests fast.

data(rads2, package = "dmrfit")
X <- rads2[, c("D3", "D6", "D7")]
quiet_fit <- function(...) suppressWarnings(dmrfit_bayes(X, nsim = 600, burnin = 300, progress = FALSE, ...))

# --- argument checks
expect_error(dmrfit_bayes(X, method = "pseudo"))
expect_error(dmrfit_bayes(X, scale = "other"))
expect_error(dmrfit_bayes(X, control = list(foo = 1)), "Unknown 'control' settings: foo")
expect_error(dmrfit_bayes(rads2[, 1:12], method = "exact"), "control\\$max_states")
expect_warning(dmrfit_bayes(X, method = "exact", scale = "mch", nsim = 50, burnin = 50, progress = FALSE), "only used when method")
expect_warning(dmrfit_bayes(X, sigma2 = 2, nsim = 50, burnin = 50, progress = FALSE), "greater than 1.0")

# --- every method and rescaling covariance runs and returns finite draws (parameters by rows)
fits <- list(core = quiet_fit(), adacore = quiet_fit(method = "adacore"), exact = quiet_fit(method = "exact"),
             dmh = quiet_fit(method = "dmh", control = list(dmh_aux = 300)),
             mch = quiet_fit(scale = "mch", control = list(mc_draws = 3000)),
             rm = quiet_fit(scale = "rm", control = list(mc_draws = 3000, rm_iter = 5)))
for (k in names(fits)) {
    f <- fits[[k]]
    expect_inherits(f, "dmrfit_bayes", info = k)
    expect_equal(dim(f$draws), c(9L + 3L, 600L), info = k)
    expect_true(all(is.finite(f$draws)), info = k)
    expect_equal(rownames(f$draws)[10], "sigma[2,1]", info = k)
}
expect_equal(fits$core$method, "core")
expect_equal(fits$core$scale, "ghw")
expect_equal(fits$rm$scale, "rm")
expect_true(is.na(fits$exact$scale))
expect_true(fits$adacore$adacore_updates >= 0)

# --- the progress bar (written to the error stream) is shown by default and not with progress = FALSE
bar_with <- capture.output(invisible(suppressWarnings(dmrfit_bayes(X, nsim = 600, burnin = 300))), type = "message")
bar_without <- capture.output(invisible(suppressWarnings(dmrfit_bayes(X, nsim = 600, burnin = 300, progress = FALSE))), type = "message")
expect_true(any(grepl("%", bar_with)))
expect_equal(length(bar_without), 0L)

# --- very short chains: a warning, not an error, when the effective sample size cannot be computed or is small
expect_warning(dmrfit_bayes(X, nsim = 20, burnin = 10, progress = FALSE, control = list(adaptive_stage = 10)), "Savage-Dickey")

# --- the posterior targeted by the fit is stated by print() and summary()
expect_stdout(print(fits$exact), "full-likelihood posterior")
expect_stdout(print(summary(fits$core)), "coordinate-rescaled pseudo-posterior")

# --- Savage-Dickey Bayes factors and effective sample sizes
sd <- fits$core$savage_dickey
expect_equal(length(sd$bf_01), 3L)
expect_true(all(sd$bf_01[sd$zero_beyond_draws] == 1 / (10 * 600)))
expect_equal(length(sd$ess), 12L)
expect_true(all(sd$ess > 0))
expect_true(is.finite(fits$core$mess))
# too few batches of floor(sqrt(nsim)) draws for the number of parameters: no multivariate effective sample size
expect_true(is.na(suppressWarnings(dmrfit_bayes(rads2[, 1:6], nsim = 100, burnin = 100, progress = FALSE))$mess))

# --- the same seed gives the same draws; the caller's random state is left untouched
set.seed(20)
before <- .Random.seed
f1 <- quiet_fit(seed = 7)
expect_identical(.Random.seed, before)
expect_identical(quiet_fit(seed = 7)$draws, f1$draws)
expect_false(identical(quiet_fit(seed = 8)$draws, f1$draws))
