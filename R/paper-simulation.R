# Edge selection with spike-and-slab priors (paper code): data preparation, the CoRe ingredients, and a wrapper for
# the pseudo-likelihood and CoRe edge-selection samplers (src/paper-pseudo-edge.cpp, src/paper-core-edge.cpp).
# Simulation drivers follow in this file.

#' Prepare data for the edge-selection samplers
#'
#' Shifts each variable to the categories 0, ..., C - 1 (baseline: its minimum), removes rows with missing values,
#' and collapses the data to its unique response patterns (in order of first occurrence) with their frequencies.
#'
#' @param data matrix or data frame of integer responses (rows: observations, columns: variables)
#' @return a list with the shifted data, the unique patterns \code{X}, their frequencies \code{freq}, and the number
#'   of categories per variable \code{n_categories}
#' @noRd
.paper_edge_prep <- function(data) {
    data <- as.matrix(data)
    data <- data[!is.na(rowSums(data)), , drop = FALSE]
    for (i in seq_len(ncol(data))) data[, i] <- data[, i] - min(data[, i])
    key <- apply(data, 1, paste, collapse = "_")
    first <- !duplicated(key)
    list(data = data,
         X = data[first, , drop = FALSE],
         freq = as.numeric(table(factor(key, levels = key[first]))),
         n_categories = apply(data, 2, max) + 1)
}

#' Data with the cross-product columns used by the optimizer and dmrf_deriv
#'
#' @param data matrix of responses coded 0, ..., C - 1
#' @return \code{data} with the columns 2 x_i x_j appended, in dmrfit's interaction order (column-major lower triangle)
#' @noRd
.paper_crossproducts <- function(data) {
    idx <- which(lower.tri(diag(ncol(data))), arr.ind = TRUE)
    cbind(data, 2.0 * data[, idx[, 1], drop = FALSE] * data[, idx[, 2], drop = FALSE])
}

#' Sufficient statistics of the full likelihood
#'
#' @param data matrix of responses coded 0, ..., C - 1
#' @param n_categories number of categories per variable
#' @return the sufficient statistics of the data in dmrfit's parameter order: the count of each category h > 0 of
#'   each variable, then the sums of x_i x_j over the observations (column-major lower triangle)
#' @noRd
.paper_data_statistics <- function(data, n_categories) {
    P <- ncol(data)
    S_thr <- unlist(lapply(seq_len(P), function(p) vapply(seq_len(n_categories[p] - 1),
                                                           function(h) sum(data[, p] == h), numeric(1))))
    idx <- which(lower.tri(diag(P)), arr.ind = TRUE)
    S_int <- colSums(data[, idx[, 1], drop = FALSE] * data[, idx[, 2], drop = FALSE])
    as.numeric(c(S_thr, S_int))
}

#' All response patterns and the sufficient statistics of the exact likelihood
#'
#' @param data matrix of responses coded 0, ..., C - 1
#' @param n_categories number of categories per variable
#' @return a list with \code{Y}, all response patterns (integer matrix), and \code{S}, the sufficient statistics of
#'   the data (\code{.paper_data_statistics()})
#' @noRd
.paper_exact_statistics <- function(data, n_categories) {
    Y <- as.matrix(expand.grid(lapply(n_categories, function(m) 0:(m - 1))))
    storage.mode(Y) <- "integer"
    dimnames(Y) <- NULL
    list(Y = Y, S = .paper_data_statistics(data, n_categories))
}

#' Anchor and rescaling matrix of the CoRe edge-selection sampler
#'
#' The anchor \code{eta_star} is the pseudo-posterior mode of the full model (default priors). The rescaling matrix
#' is \code{A = Gamma L'}, with \code{L L' = -H} (\code{H} the Hessian of the log pseudo-likelihood at
#' \code{eta_star}) and \code{Gamma Gamma' = V}, so that the adjusted pseudo-likelihood
#' \code{l_PL(eta_star + A^{-1} (beta - eta_star))} has curvature \code{V^{-1}} at \code{eta_star}. All quantities are
#' likelihood-only: the prior of the edge-selection model is placed on \code{beta}.
#'
#' @param data matrix of responses coded 0, ..., C - 1 (as returned by \code{.paper_edge_prep()})
#' @param n_categories number of categories per variable
#' @param scale target covariance \code{V}: \code{"godambe_mc"} (default) the expected Godambe covariance
#'   \code{Hbar^{-1} J Hbar^{-1}} of the pseudo-likelihood estimator under the fitted model, with \code{Hbar} and the
#'   centered \code{J} from one simulation of \code{mc_draws} observations (\code{godambe_gibbs_iter} Gibbs iterations from
#'   resampled observed rows). \code{"ghw"}: the data-based Godambe-Huber-White covariance. \code{"mch"}: the inverse
#'   Fisher information of the full likelihood, by Monte Carlo (\code{mc_draws} observations, \code{mc_gibbs_iter}
#'   Gibbs iterations), without the prior curvature
#' @param v_at point at which the model-based \code{V} is evaluated (\code{"godambe_mc"}, \code{"mch"}). Default
#'   \code{eta_star}. In a simulation study, the true parameter gives the oracle \code{V}.
#' @param seed seed of the Monte Carlo simulation
#' @return a list with \code{eta_star}, \code{A}, \code{Ainv}, \code{H} (the negative Hessian), \code{V} and the scale
#' @noRd
paper_edge_ingredients <- function(data, n_categories, scale = c("godambe_mc", "ghw", "mch"), parinit = NULL,
                                   mc_draws = 1e5, mc_gibbs_iter = 5, godambe_gibbs_iter = 20, v_at = NULL, seed = 1) {
    scale <- match.arg(scale)
    P <- ncol(data)
    n_thr <- sum(n_categories - 1)
    n_pars <- n_thr + P * (P - 1) / 2
    if (is.null(parinit)) parinit <- rep(0, n_pars)
    data_cp <- .paper_crossproducts(data)

    # --- anchor: pseudo-posterior mode of the full model ---
    fit <- cpp_optimize(data = data_cp, parinit = parinit, n_categories = n_categories, P = P,
                        f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps),
                        n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = TRUE, epsilon = 1e-06,
                        ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5,
                        interactions_location = 0.0, interactions_scale = 2.5)
    eta_star <- as.vector(fit$argument)

    # --- negative Hessian of the log pseudo-likelihood at the anchor (likelihood only) ---
    d <- dmrf_deriv(pars = eta_star, data = data_cp, P = P, n_categories = n_categories, with_prior = FALSE,
                    ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5,
                    interactions_location = 0.0, interactions_scale = 2.5)
    H <- d$hessian
    at <- if (is.null(v_at)) eta_star else v_at

    # --- target covariance V ---
    if (scale == "ghw") {
        if (!is.null(v_at)) stop("v_at applies only to the model-based scales.")
        V <- d$HW
    } else if (scale == "godambe_mc") {
        set.seed(seed)
        n <- nrow(data)
        mu <- matrix(0, P, max(n_categories) - 1)
        pos <- 0
        for (p in seq_len(P)) for (h in seq_len(n_categories[p] - 1)) { pos <- pos + 1; mu[p, h] <- at[pos] }
        sigma <- matrix(0, P, P)
        sigma[lower.tri(sigma)] <- at[-seq_len(n_thr)]
        sigma <- sigma + t(sigma)
        start <- data[sample.int(n, mc_draws, replace = TRUE), , drop = FALSE] # warm start at observed rows
        storage.mode(start) <- "integer"
        Xs <- cpp_gibbs_sampler_omrf(mu = mu, sigma = sigma, n_categories = n_categories, N = mc_draws, P = P,
                                     iter = godambe_gibbs_iter, X_start = start, save_iter = FALSE)$X
        ds <- dmrf_deriv(pars = at, data = .paper_crossproducts(Xs), P = P, n_categories = n_categories,
                         with_prior = FALSE, ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5,
                         interactions_location = 0.0, interactions_scale = 2.5)
        H_bar <- ds$hessian * n / mc_draws
        g_bar <- ds$gradient / mc_draws # mean score (its sign is irrelevant below)
        J <- (ds$FisherInfo - mc_draws * tcrossprod(g_bar)) * n / mc_draws
        H_bar_inv <- solve(H_bar)
        V <- H_bar_inv %*% J %*% H_bar_inv
    } else {
        set.seed(seed)
        H_mc <- cpp_compute_mc_hessian(data = data, pars = at, n_categories = n_categories, L = mc_draws,
                                       sampler_n_iter = mc_gibbs_iter, thresholds_alpha = 0.5,
                                       thresholds_beta = 0.5, interactions_location = 0.0,
                                       interactions_scale = 2.5)
        # cpp_compute_mc_hessian adds the log-prior curvature on the diagonal: remove it (likelihood only)
        thr <- at[seq_len(n_thr)]
        s <- stats::plogis(thr)
        int <- at[-seq_len(n_thr)]
        prior_d2 <- c(-(0.5 + 0.5) * s * (1 - s), 2 * (int^2 - 2.5^2) / (2.5^2 + int^2)^2)
        V <- solve(-(H_mc - diag(prior_d2)))
    }
    V <- (V + t(V)) / 2
    A <- t(chol(V)) %*% chol(H)
    list(eta_star = eta_star, A = A, Ainv = solve(A), H = H, V = V, scale = scale)
}

#' Spike-and-slab edge selection with the pseudo-likelihood, CoRe, exact-likelihood or DMH-MoMS sampler
#'
#' @param data matrix or data frame of integer responses
#' @param method \code{"core"} (CoRe sampler), \code{"pseudo"} (pseudo-likelihood sampler), \code{"exact"}
#'   (exact likelihood, normalizing constant by enumerating all response patterns, at most \code{max_states}
#'   patterns) or \code{"dmh"} (DMH-MoMS: full likelihood, ratio of normalizing constants estimated from auxiliary
#'   draws as in dmrfit's DMH sampler)
#' @param nsim number of kept iterations. \code{burnin}: burn-in iterations (phases 1-3, edge selection from
#'   \code{burnin / 2}). \code{adaptive_stage}: iterations of phase 1. \code{thin}: thinning
#' @param edge_prior \code{"beta-bernoulli"} (\code{beta_alpha}, \code{beta_beta}) or \code{"bernoulli"}
#'   (\code{inclusion_prob})
#' @param slab \code{"cauchy"} or \code{"normal"}, with scale \code{slab_scale}
#' @param threshold_prior \code{"betaprime"} or \code{"normal"} (\code{threshold_mean}, \code{threshold_sd})
#' @param ingredients output of \code{paper_edge_ingredients()} (computed with \code{scale} if \code{NULL})
#' @param pars initial parameters (default: the anchor \code{eta_star})
#' @param max_states largest number of response patterns for \code{method = "exact"} (default 2^20, i.e. 10
#'   variables with 4 categories)
#' @param dmh_c for \code{method = "dmh"}, the number of auxiliary draws per observation: \code{L = dmh_c * n}.
#'   \code{dmh_L}: the number of auxiliary draws, which overrides \code{dmh_c}. \code{dmh_gibbs_iter}: the number of
#'   Gibbs iterations per auxiliary draw
#' @param seed seed of the sampler
#' @return the sampler output with the draws (iterations by rows), the inclusion probabilities (\code{inclusion},
#'   frequency of presence, and \code{inclusion_rb}, Rao-Blackwellized) and the edge indicators \code{gamma}
#' @noRd
paper_edge_sampler <- function(data, method = c("core", "pseudo", "exact", "dmh"), nsim = 2000, burnin = 5000, adaptive_stage = 500,
                               thin = 1, sigma2 = 0.1, edge_prior = c("beta-bernoulli", "bernoulli"),
                               beta_alpha = 1, beta_beta = 1, inclusion_prob = 0.5,
                               slab = c("cauchy", "normal"), slab_scale = 2.5,
                               threshold_prior = c("betaprime", "normal"), threshold_mean = NULL, threshold_sd = 1,
                               scale = c("godambe_mc", "ghw", "mch"), ingredients = NULL, pars = NULL,
                               use_likelihood = TRUE, max_states = 2^20, dmh_c = 50, dmh_L = NULL, dmh_gibbs_iter = 5,
                               seed = NULL, progress = TRUE) {
    method <- match.arg(method)
    edge_prior <- match.arg(edge_prior)
    slab <- match.arg(slab)
    threshold_prior <- match.arg(threshold_prior)
    scale <- match.arg(scale)
    d <- .paper_edge_prep(data)
    P <- ncol(d$X)
    if (P < 2) stop("edge selection needs at least 2 variables.")
    n_thr <- sum(d$n_categories - 1)
    n_int <- P * (P - 1) / 2
    if (method == "exact" && prod(as.numeric(d$n_categories)) > max_states)
        stop("method = \"exact\" enumerates ", prod(as.numeric(d$n_categories)), " response patterns, more than max_states.")
    if (is.null(ingredients)) {
        # the full-likelihood samplers need only the anchor (starting values). "ghw" avoids the Monte Carlo simulation
        ingredients <- paper_edge_ingredients(d$data, d$n_categories,
                                              scale = if (method %in% c("exact", "dmh")) "ghw" else scale)
    }
    if (is.null(pars)) pars <- ingredients$eta_star
    if (!is.null(seed)) set.seed(seed)
    args <- list(X = d$X, freq = d$freq, n_categories = d$n_categories, pars = pars, nsim = nsim, burnin = burnin,
                 adaptive_stage_n_iter = adaptive_stage, thin = thin, sigma2 = sigma2,
                 edge_prior = edge_prior, inclusion_prob = inclusion_prob, beta_alpha = beta_alpha,
                 beta_beta = beta_beta, slab = slab, interactions_location = 0, interactions_scale = slab_scale,
                 threshold_prior = threshold_prior, thresholds_alpha = 0.5, thresholds_beta = 0.5,
                 threshold_mean = threshold_mean, threshold_sd = threshold_sd,
                 use_likelihood = use_likelihood, progress = progress)
    out <- if (method == "core") {
        do.call(cpp_core_edge_sampler, c(args, list(eta_star = ingredients$eta_star, Ainv = ingredients$Ainv)))
    } else if (method == "pseudo") {
        do.call(cpp_pseudo_edge_sampler, args)
    } else if (method == "exact") {
        ex <- .paper_exact_statistics(d$data, d$n_categories)
        args$freq <- NULL
        do.call(cpp_exact_edge_sampler, c(args, list(Y = ex$Y, S = ex$S, n = nrow(d$data))))
    } else {
        L <- if (is.null(dmh_L)) ceiling(dmh_c * nrow(d$data)) else dmh_L
        do.call(cpp_dmh_edge_sampler, c(args, list(S = .paper_data_statistics(d$data, d$n_categories), data = d$data,
                                                   n = nrow(d$data), L = L, inner_sampler_n_iter = dmh_gibbs_iter)))
    }
    colnames(out$edge_counts) <- c("add_try", "add_acc", "del_try", "del_acc", "within_try", "within_acc")
    gamma <- out$draws[, n_thr + seq_len(n_int), drop = FALSE] != 0
    out$inclusion <- colMeans(gamma)
    out$inclusion_rb <- as.vector(out$inclusion_rb)
    out$gamma <- gamma
    out$method <- method
    out$n_thresholds <- n_thr
    out$n_categories <- d$n_categories
    out
}
