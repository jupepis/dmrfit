# Tests for the C++ samplers called by dmrfit_bayes(): output format, valid draws, reproducibility with a seed, and the
# scale of the CoRe target. Short chains on a three-item network keep the tests fast.

data(rads2, package = "dmrfit")
X <- as.matrix(rads2[, c("D3", "D6", "D7")]) - 1L
P <- ncol(X)
n_categories <- rep(4L, P)
n_thresholds <- sum(n_categories - 1)
n_pars <- n_thresholds + P * (P - 1) / 2
X_suff <- cbind(X, 2 * t(apply(X, 1, function(x) { S <- x %*% t(x); S[lower.tri(S)] })))
fit <- dmrfit:::cpp_optimize(data = X_suff, parinit = rep(0, n_pars), n_categories = n_categories, P = P,
                             f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100,
                             rinit = 1, rmax = 10, with_prior = TRUE, epsilon = 1e-06, ncores = 1L, thresholds_alpha = 0.5,
                             thresholds_beta = 0.5, interactions_location = 0, interactions_scale = 2.5)
mode <- as.vector(fit$argument)
current_scale <- chol(fit$utils$hessian)
new_scale <- t(chol(fit$utils$HW))
patterns <- as.matrix(expand.grid(lapply(n_categories, function(m) 0:(m - 1))))
X_states <- dmrfit:::cpp_build_permutations_stats(permutations = patterns, n_pars = n_pars, n_thresholds = n_thresholds,
                                                  n_categories = n_categories)
nsim <- 1500
args <- list(data = X, pars = mode, n_categories = n_categories, nsim = nsim, burnin = 500, progress = FALSE)
run <- list(
    core    = function() do.call(dmrfit:::cpp_core_sampler, c(args, list(pmles = mode, current_scale = current_scale, new_scale = new_scale))),
    adacore = function() do.call(dmrfit:::cpp_adacore_sampler, c(replace(args, "data", list(X_suff)), list(pmles = mode))),
    exact   = function() do.call(dmrfit:::cpp_exact_sampler, c(args, list(X = X_states))),
    dmh     = function() do.call(dmrfit:::cpp_dmh_sampler, c(args, list(L = 300L))),
    pseudo  = function() do.call(dmrfit:::cpp_pseudo_sampler, args)
)

# --- every sampler: draws by rows (one per kept iteration), finite, with an acceptance rate in (0, 1)
out <- lapply(run, function(f) { set.seed(5); f() })
for (k in names(out)) {
    expect_equal(dim(out[[k]]$draws), c(nsim, n_pars), info = k)
    expect_true(all(is.finite(out[[k]]$draws)), info = k)
    expect_true(out[[k]]$acceptance > 0 && out[[k]]$acceptance < 1, info = k)
}
expect_true(out$adacore$counter_update >= 0)

# --- the same seed gives the same draws (the samplers use R's random number generator)
set.seed(5)
expect_identical(run$core()$draws, out$core$draws)

# --- CoRe targets the GHW covariance; the pseudo-posterior is narrower (the deficit CoRe corrects)
sd_core <- apply(out$core$draws, 2, sd)
sd_pseudo <- apply(out$pseudo$draws, 2, sd)
expect_true(abs(median(sd_core / sqrt(diag(fit$utils$HW))) - 1) < 0.2)
expect_true(median(sd_pseudo / sd_core) < 0.95)

# --- the exact posterior and CoRe agree in location (both centred close to the pseudo-posterior mode)
expect_true(median(abs(colMeans(out$core$draws) - colMeans(out$exact$draws)) / apply(out$exact$draws, 2, sd)) < 0.3)
