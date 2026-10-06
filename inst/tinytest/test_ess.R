# Tests for the effective sample sizes reported by dmrfit_bayes(): .ess_mcmc() (one parameter, autocorrelation-based)
# and .mess_mcmc() (multivariate, lugsail batch means). Checked on chains whose effective sample size is known.

ess_mcmc <- dmrfit:::.ess_mcmc
mess_mcmc <- dmrfit:::.mess_mcmc

# --- .ess_mcmc: independent draws give about n, an AR(1) chain about n (1 - phi) / (1 + phi)
set.seed(1)
n <- 4000
expect_true(abs(ess_mcmc(rnorm(n)) / n - 1) < 0.15)
ar <- as.numeric(stats::arima.sim(list(ar = 0.9), n))
expect_true(abs(ess_mcmc(ar) / (n * 0.1 / 1.9) - 1) < 0.5)
expect_true(ess_mcmc(ar) < ess_mcmc(rnorm(n)))

# constant or too short chains have no effective sample size
expect_true(is.na(ess_mcmc(rep(1, 100))))
expect_true(is.na(ess_mcmc(c(1, 2, 3))))

# --- .mess_mcmc: draws are passed as a p x n matrix (parameters by rows)
p <- 5
iid <- matrix(rnorm(p * n), p, n)
expect_true(abs(mess_mcmc(iid) / n - 1) < 0.3)
ar_draws <- t(sapply(seq_len(p), function(j) as.numeric(stats::arima.sim(list(ar = 0.9), n))))
expect_true(mess_mcmc(ar_draws) < mess_mcmc(iid))
expect_true(abs(mess_mcmc(ar_draws) / (n * 0.1 / 1.9) - 1) < 0.5)

# the asymptotic covariance needs more batches (floor(sqrt(n)) draws each) than parameters
expect_true(is.na(mess_mcmc(matrix(rnorm(80 * 400), 80, 400))))
