# Tests for the cpp_* utilities that have no exported R wrapper. They are called by dmrfit() and dmrfit_bayes(), so
# a silent change inside them (indexing of the parameters, sign conventions) would not raise any error at R level.

# a small ordinal network: three RADS-2 items, four categories each, recoded to 0, ..., 3
data(rads2, package = "dmrfit")
X <- as.matrix(rads2[, c("D3", "D6", "D7")]) - 1L
P <- ncol(X)
n_categories <- rep(4L, P)
n_thresholds <- sum(n_categories - 1)
n_pars <- n_thresholds + P * (P - 1) / 2
# data with the cross-product columns 2 x_i x_j, as dmrfit() passes it to the pseudo-likelihood functions
X_suff <- cbind(X, 2 * t(apply(X, 1, function(x) { S <- x %*% t(x); S[lower.tri(S)] })))

# --- cpp_npseudologlik_draws: same value as cpp_npseudologlik at every draw, with and without prior
set.seed(3)
draws <- matrix(rnorm(n_pars * 25, 0, 0.5), n_pars, 25)
for (with_prior in c(TRUE, FALSE)) {
    batched <- dmrfit:::cpp_npseudologlik_draws(draws, X, n_categories, with_prior, 0.5, 0.5, 0, 2.5)
    single <- apply(draws, 2, function(pars) dmrfit:::cpp_npseudologlik(pars, X_suff, P, n_categories, with_prior, 1L,
                                                                       0.5, 0.5, 0, 2.5))
    expect_equal(as.vector(batched), single, tolerance = 1e-10)
}

# duplicated rows only change the frequencies of the unique response patterns, not the value per observation
doubled <- dmrfit:::cpp_npseudologlik_draws(draws[, 1:3], rbind(X, X), n_categories, FALSE, 0.5, 0.5, 0, 2.5)
once <- dmrfit:::cpp_npseudologlik_draws(draws[, 1:3], X, n_categories, FALSE, 0.5, 0.5, 0, 2.5)
expect_equal(as.vector(doubled), 2 * as.vector(once), tolerance = 1e-10)

# the stabilized log-sum-exp keeps the value finite at large parameter values
expect_true(all(is.finite(dmrfit:::cpp_npseudologlik_draws(matrix(30, n_pars, 1), X, n_categories, TRUE, 0.5, 0.5, 0, 2.5))))

# --- dmrf_deriv: curvature of the negative log pseudo-posterior and GHW covariance at the mode
fit <- dmrfit:::cpp_optimize(data = X_suff, parinit = rep(0, n_pars), n_categories = n_categories, P = P,
                             f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100,
                             rinit = 1, rmax = 10, with_prior = TRUE, epsilon = 1e-06, ncores = 1L, thresholds_alpha = 0.5,
                             thresholds_beta = 0.5, interactions_location = 0, interactions_scale = 2.5)
d <- dmrfit:::dmrf_deriv(as.vector(fit$argument), X_suff, P, n_categories, TRUE, 1L, 0.5, 0.5, 0, 2.5)
expect_equal(dim(d$hessian), c(n_pars, n_pars))
expect_equal(d$hessian, t(d$hessian), tolerance = 1e-10)
expect_equal(d$HW, t(d$HW), tolerance = 1e-10)
# both are positive definite (the CoRe samplers take their Cholesky factors)
expect_true(all(eigen(d$hessian, symmetric = TRUE, only.values = TRUE)$values > 0))
expect_true(all(eigen(d$HW, symmetric = TRUE, only.values = TRUE)$values > 0))
# the optimizer stops close to the mode: one more Newton step moves each parameter by less than 0.1 posterior SDs
# (the gradient itself is a sum over the observations, so its raw size is not a useful criterion)
newton_step <- solve(d$hessian, d$gradient)
expect_true(max(abs(newton_step) / sqrt(diag(solve(d$hessian)))) < 0.1)

# --- cpp_mvnrnd_arma: draws by columns, with the requested covariance
set.seed(4)
Sigma <- matrix(c(1, 0.5, 0.5, 2), 2)
Z <- dmrfit:::cpp_mvnrnd_arma(c(0, 0), Sigma, 20000L)
expect_equal(dim(Z), c(2L, 20000L))
expect_equal(cov(t(Z)), Sigma, tolerance = 0.05)

# --- cpp_build_permutations_stats: one row per response pattern, one column per parameter
patterns <- as.matrix(expand.grid(lapply(n_categories, function(m) 0:(m - 1))))
S <- dmrfit:::cpp_build_permutations_stats(permutations = patterns, n_pars = n_pars, n_thresholds = n_thresholds,
                                           n_categories = n_categories)
expect_equal(dim(S), c(prod(n_categories), n_pars))
# each pattern has one threshold indicator per variable not in the baseline category
expect_equal(rowSums(S[, seq_len(n_thresholds)]), rowSums(patterns > 0))

# --- Gibbs sampler: without interactions the nodes are independent, with category probabilities proportional to
# exp(threshold of the category), the baseline category 0 having threshold 0
mu <- matrix(c(-1, 0.5, 0, 1, -0.5, 0.5), nrow = 2, byrow = TRUE)
set.seed(1)
gibbs <- dmrfit:::cpp_gibbs_sampler_omrf(mu = mu, sigma = matrix(0, 2, 2), n_categories = c(4L, 4L), N = 20000L,
                                         P = 2L, iter = 5L, X_start = NULL, save_iter = FALSE)
expect_equal(dim(gibbs$X), c(20000L, 2L))
for (p in 1:2) {
    expect_equal(as.vector(table(factor(gibbs$X[, p], levels = 0:3))) / 20000,
                 exp(c(0, mu[p, ])) / sum(exp(c(0, mu[p, ]))), tolerance = 0.02)
}

