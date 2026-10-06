#' @title dmrfit_bayes
#'
#' @description Bayesian inference for a discrete Markov Random Field (Ising or ordinal) by Markov chain Monte Carlo.
#' By default, the posterior is approximated by the coordinate-rescaled (CoRe) pseudo-posterior of Arena and Marsman
#' (2026): the pseudo-posterior rescaled around its mode so that its covariance matches the sandwich
#' (Godambe-Huber-White) covariance, which corrects the underestimated posterior variability of the pseudo-likelihood.
#' All samplers use an adaptive Fisher-preconditioned Metropolis-adjusted Langevin algorithm (FisherMALA).
#'
#' @param data data matrix, with rows as samples and columns as variables. Each variable should be rescaled to the range of 0 to m-1, where m is the number of categories for that variable. The baseline category is always the minimum value in the variable. The internal processing will check if the variables are rescaled and will rescale them if necessary. If there are any NAs in the data, they will be removed before optimization (listwise deletion).
#' @param parinit parinit, initial parameter estimates for the optimization that finds the pseudo-posterior mode, a vector of length equal to the number of parameters in the model. If NULL, it will be initialized to a vector of zeros. The number of parameters in the model is calculated as sum(P * (n_categories - 1)) + P * (P - 1) / 2, where P is the number of variables.
#' @param method the sampler: \code{"core"} (default) samples the CoRe pseudo-posterior with a fixed rescaling
#'   matrix; \code{"adacore"} adapts the rescaling matrix during burn-in at the running mean of the draws;
#'   \code{"exact"} samples the posterior based on the full likelihood, whose normalizing constant is computed by
#'   enumerating all response patterns (only for small networks, see \code{control}); \code{"dmh"} samples the
#'   posterior based on the full likelihood with the double Metropolis-Hastings algorithm (Liang, 2010). DMH simulates auxiliary data at every
#'   iteration and is much slower than the other methods (minutes to hours instead of seconds, depending on the
#'   network size and \code{control$dmh_aux}); \code{"core"} runs at the cost of the pseudo-posterior sampler.
#' @param scale the covariance to which the pseudo-posterior is rescaled, only used when \code{method = "core"}:
#'   \code{"ghw"} (default) the sandwich (Godambe-Huber-White) covariance; \code{"mch"} the inverse of the negative
#'   Hessian of the full log-posterior at the pseudo-posterior mode, estimated by Monte Carlo simulation; \code{"rm"}
#'   the same Hessian at the full-posterior mode, estimated by a Newton-type Robbins-Monro algorithm (the sampler
#'   remains centred at the pseudo-posterior mode). \code{"mch"} and \code{"rm"} simulate data from the model and are
#'   slower.
#' @param nsim number of posterior draws kept after burn-in. Default is 1000.
#' @param burnin number of burn-in iterations, after an initial adaptive stage of \code{control$adaptive_stage} iterations. Default is 1000.
#' @param ncores ncores, number of cores to use for parallel computing of the gradient, hessian and likelihood values of the optimization. If ncores is larger than the number of available cores, it will be set to the number of available cores minus one. Default is 1.
#' @param thresholds_alpha alpha parameter for the Beta-Prime prior on thresholds (default is 0.5).
#' @param thresholds_beta beta parameter for the Beta-Prime prior on thresholds (default is 0.5).
#' @param interactions_location location parameter for the Cauchy prior on pairwise interactions (default is 0.0).
#' @param interactions_scale scale parameter for the Cauchy prior on pairwise interactions (default is 2.5).
#' @param sigma2 initial value for the adaptive variance parameter in the FisherMALA sampler (default is 0.1).
#' @param seed random seed for reproducibility of the MCMC sampler (default is 123). The caller's random number stream is restored on exit.
#' @param progress logical, whether to show a progress bar while sampling (default is TRUE).
#' @param control a list of tuning settings, each with a default: \code{adaptive_stage} (500) iterations of the
#'   initial adaptive stage; \code{max_states} (1e6) largest number of response patterns for \code{method = "exact"};
#'   \code{dmh_aux} (25000) and \code{dmh_gibbs_iter} (5) auxiliary draws and Gibbs sweeps per iteration of
#'   \code{method = "dmh"}; \code{mc_draws} (1e5) and \code{mc_gibbs_iter} (5) simulated observations and Gibbs sweeps
#'   for \code{scale = "mch"} and \code{"rm"}; \code{rm_iter} (50) Robbins-Monro iterations for \code{scale = "rm"}.
#'
#' @return an object of class \code{dmrfit_bayes} (also of class \code{dmrfit}), including the posterior draws (\code{draws}), the Savage-Dickey Bayes factors with the effective sample size of each parameter (\code{savage_dickey}), and the multivariate effective sample size of the draws (\code{mess}; Vats, Flegal and Jones, 2019), which is \code{NA} when there are fewer than \code{P + 1} batches of \code{floor(sqrt(nsim))} draws per parameter.
#'
#' @references Arena, G. and Marsman, M. (2026). Bayesian inference for discrete Markov random fields through
#' coordinate rescaling. Manuscript submitted for publication.
#'
#' Liang, F. (2010). A double Metropolis-Hastings sampler for spatial models with intractable normalizing constants.
#' \emph{Journal of Statistical Computation and Simulation}, 80(9), 1007-1022.
#'
#' @examples
#'
#' # simulate data from a 3-node Ising model with 2 categories per node
#' set.seed(123)
#' n <- 1000
#' P <- 3
#' data <- matrix(rbinom(n * P, size = 1, prob = c(0.2, 0.5, 0.3)[rep(1:P, each = n)]), nrow = n, ncol = P)
#' fit <- dmrfit_bayes(data = data, nsim = 500, burnin = 500)
#' summary(fit)
#'
#' # small network: the exact posterior is also available
#' fit_exact <- dmrfit_bayes(data = data, method = "exact", nsim = 500, burnin = 500)
#'
#' @export
#'
dmrfit_bayes <- function(data, parinit = NULL, method = c("core", "adacore", "exact", "dmh"), scale = c("ghw", "mch", "rm"), nsim = 1e03, burnin = 1e03, ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5, interactions_location = 0.0, interactions_scale = 2.5, sigma2 = 0.1, seed = 123, progress = TRUE, control = list()) {

    # save the matched call for print and summary methods
    cl <- match.call()

    # --- Random seed for the sampler (the caller's random state is restored on exit) ---
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
        if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)

    # --- Method, scale and tuning settings ---
    method <- match.arg(method)
    scale_given <- !missing(scale)
    scale <- match.arg(scale)
    if (method != "core" && scale_given && scale != "ghw")
        warning("'scale' is only used when method = \"core\"; it is ignored for method = \"", method, "\".")
    ctrl <- list(adaptive_stage = 500, max_states = 1e6, dmh_aux = 25000, dmh_gibbs_iter = 5,
                 mc_draws = 1e5, mc_gibbs_iter = 5, rm_iter = 50)
    unknown <- setdiff(names(control), names(ctrl))
    if (length(unknown)) stop("Unknown 'control' settings: ", paste(unknown, collapse = ", "), ".")
    ctrl[names(control)] <- control

    # processing input arguments

    if(sigma2 > 1.0) {
        warning("sigma2 is set to a value greater than 1.0. This may lead to unstable sampling. Consider using a smaller value (e.g., 1.0, 0.1 or 0.01) for better performance.")
    }
    if(sigma2 <= 0.0) {
        sigma2 <- 1.0
        warning("sigma2 must be positive. It is reset to 1.0.")
    }

    # check if data is a matrix, if not convert it to a matrix
    if(!is.matrix(data)) {
        data <- as.matrix(data)
    }

     # remove NAs from data if any exist
    data <- data[!is.na(rowSums(data)), ]

    # check that columns of data are integer and non-negative
    if(any(data < 0) || any(data != floor(data))) {
        stop("All values in data must be non-negative integers.")
    }

    # rescale the data to the range of 0 to m-1 if necessary
    for(i in 1:ncol(data)) {
        data[, i] <- data[, i] - min(data[, i]) # baseline category is always the minimum value in the variable
    }

    # n_categories is calculated as the maximum value in each column of data plus 1, since the categories are assumed to be coded from 0 to m-1
    n_categories <- apply(data, 2, max) + 1
    # number of variables
    P <- ncol(data)
    n_thresholds <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2
    n_pars <- n_thresholds + n_interactions

    # --- Exact sampling enumerates all response patterns: only for small networks ---
    n_states <- prod(as.numeric(n_categories))
    if (method == "exact" && n_states > ctrl$max_states)
        stop("method = \"exact\" enumerates all ", format(n_states, big.mark = ",", scientific = FALSE), " response patterns, ",
             "more than control$max_states = ", format(ctrl$max_states, big.mark = ",", scientific = FALSE), ". ",
             "Use method = \"core\" or \"dmh\", or increase control$max_states (memory and time grow with the number of patterns).")

    # if parinit is not provided, initialize it to a vector of zeros
    if(is.null(parinit)) {
        parinit <- rep(0.0, n_pars)
    }
    else if(length(parinit) != n_pars) {
        stop(paste("Length of parinit must be equal to the number of parameters in the model:", n_pars))
    }

    # number of cores: the value given by the user, capped at the number of available cores minus one
    ncores <- max(1L, min(as.integer(ncores), parallel::detectCores() - 1L), na.rm = TRUE)

    # cross-product terms for the pairwise associations
    data_raw <- data # the P variables only (input of the samplers)
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    # finding the PMLEs via optimization of the pseudo-likelihood with trust region method
    pmles <- suppressWarnings(tryCatch(expr = cpp_optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = TRUE, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    if(is.null(pmles)) {
        warning("Optimization failed. Returning NULL.")
        return(NULL)
    }
    mode_pseudo <- as.vector(pmles$argument)

    # --- Sampler arguments shared by all methods ---
    args <- list(data = data_raw, pars = mode_pseudo, n_categories = n_categories, nsim = nsim, burnin = burnin,
                 adaptive_stage_n_iter = ctrl$adaptive_stage, sigma2 = sigma2,
                 thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                 interactions_location = interactions_location, interactions_scale = interactions_scale,
                 verbose = FALSE, progress = progress)

    # --- Run the sampler ---
    out <- tryCatch({
        if (method == "core") {
            # rescaling matrices: Cholesky factors of the pseudo-posterior curvature and of the target covariance
            current_scale <- chol(pmles$utils$hessian)
            target_cov <- switch(scale,
                ghw = pmles$utils$HW,
                mch = qr.solve(-cpp_compute_mc_hessian(data = data_raw, pars = mode_pseudo, n_categories = n_categories,
                                                       L = ctrl$mc_draws, sampler_n_iter = ctrl$mc_gibbs_iter,
                                                       thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                                                       interactions_location = interactions_location,
                                                       interactions_scale = interactions_scale)),
                rm = qr.solve(-cpp_compute_robbins_monro(data = data_raw, pars_init = mode_pseudo, n_categories = n_categories,
                                                         thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                                                         interactions_location = interactions_location,
                                                         interactions_scale = interactions_scale, L = ctrl$mc_draws,
                                                         sampler_n_iter = ctrl$mc_gibbs_iter, rm_max_iter = ctrl$rm_iter)$hessian))
            new_scale <- t(chol((target_cov + t(target_cov)) / 2))
            do.call(cpp_core_sampler, c(args, list(pmles = mode_pseudo, current_scale = current_scale, new_scale = new_scale)))
        } else if (method == "adacore") {
            # AdaCoRe takes the data with the cross-product columns (it recomputes the sandwich covariance)
            args_ada <- args
            args_ada$data <- data
            do.call(cpp_adacore_sampler, c(args_ada, list(pmles = mode_pseudo)))
        } else if (method == "exact") {
            patterns <- as.matrix(expand.grid(lapply(n_categories, function(m) 0:(m - 1))))
            X <- cpp_build_permutations_stats(permutations = patterns, n_pars = n_pars, n_thresholds = n_thresholds,
                                              n_categories = n_categories)
            do.call(cpp_exact_sampler, c(args, list(X = X)))
        } else {
            do.call(cpp_dmh_sampler, c(args, list(L = ctrl$dmh_aux, inner_sampler_n_iter = ctrl$dmh_gibbs_iter)))
        }
    }, error = function(e) {
        warning("Sampler (method = \"", method, "\") failed: ", conditionMessage(e), " Returning NULL.")
        NULL
    })
    if(is.null(out)) return(NULL)
    out$draws <- t(out$draws) # parameters by rows

    # metadata needed by print/summary
    pmles$call <- cl
    pmles$P <- P
    pmles$n_categories <- n_categories
    pmles$N <- nrow(data)
    pmles$with_prior <- TRUE
    pmles$structured <- FALSE
    pmles$structure <- NULL
    pmles$ncores <- ncores
    pmles$acceptance <- out$acceptance
    pmles$seconds_elapsed <- out$seconds_elapsed
    pmles$method <- method
    pmles$scale <- if (method == "core") scale else NA_character_
    if (!is.null(out$counter_update)) pmles$adacore_updates <- out$counter_update # updates of the rescaling matrix during burn-in

    # label the parameter vector
    n_thresholds <- sum(n_categories - 1)
    thresh_names <- unlist(lapply(seq_len(P), function(p) {
        paste0("mu[", p, ",", seq_len(n_categories[p] - 1), "]")
    }))
    inter_names <- unlist(lapply(1:(P - 1), function(j) {
        lapply((j + 1):P, function(i) paste0("sigma[", i, ",", j, "]"))
    }))
    par_names <- c(thresh_names, inter_names)
    pmles$argument <- setNames(as.vector(pmles$argument), par_names) # plain named vector (the optimizer returns a one-column matrix)

    # label the posterior draws
    rownames(out$draws) <- par_names

    # posterior summary from MCMC draws
    inter_idx <- (n_thresholds + 1):n_pars
    pars <- pmles$argument
    se <- apply(out$draws, 1, sd)
    names(se) <- par_names

    # Savage-Dickey via posterior draws
    log_prior_at_zero <- dcauchy(0, location = interactions_location, scale = interactions_scale, log = TRUE)

    bf_01 <- numeric(length(inter_idx))
    names(bf_01) <- par_names[inter_idx]

    # --- Savage-Dickey density ratio at zero for each interaction ---
    # When zero lies beyond all posterior draws, the density at zero cannot be estimated from the draws (a kernel
    # estimate there only extrapolates the tail of the nearest kernel): BF_01 is then reported at a floor and flagged
    bf_floor <- 1 / (10 * nsim)
    zero_beyond_draws <- logical(length(inter_idx))
    names(zero_beyond_draws) <- names(bf_01)

    for (k in seq_along(inter_idx)) {
        j <- inter_idx[k]
        draws_j <- out$draws[j, ]
        zero_beyond_draws[k] <- 0 < min(draws_j) || 0 > max(draws_j)
        if (zero_beyond_draws[k]) { bf_01[k] <- bf_floor; next }
        bf_01[k] <- tryCatch({
            # Gaussian kernel density estimate of the marginal posterior at zero, on the log scale (log-sum-exp)
            h <- stats::bw.nrd0(draws_j)
            log_k <- stats::dnorm(0, mean = draws_j, sd = h, log = TRUE)
            log_post_at_zero <- max(log_k) + log(mean(exp(log_k - max(log_k))))
            exp(log_post_at_zero - log_prior_at_zero)
        }, error = function(e) {
            warning("Savage-Dickey: density estimation failed for ", par_names[j],
                    ". The 'BF_01' is set to NA.")
            NA_real_
        })
    }

    pr_null <- bf_01 / (1 + bf_01) # this is the Pr(=0|x)

    # --- Effective sample size of the posterior draws, per parameter (autocorrelation-based) ---
    ess <- apply(out$draws, 1, .ess_mcmc)
    names(ess) <- par_names
    ess_inter <- ess[inter_idx]
    if (anyNA(ess_inter)) {
        # NA: the draws of an interaction never changed (the chain did not move), typically for very short chains
        warning("Savage-Dickey: the draws of ", sum(is.na(ess_inter)), " interaction(s) are constant; the Bayes factors ",
                "are unreliable. Consider increasing 'nsim' and 'burnin'.")
    } else if (min(ess_inter) < 100) {
        warning("Savage-Dickey: the smallest effective sample size of the interactions is ", round(min(ess_inter), 1),
                "; the Bayes factors may be unreliable. Consider increasing 'nsim'.")
    }

    pmles$savage_dickey <- list(
        bf_01 = bf_01,
        pr_null = pr_null,
        estimate = pars[inter_idx],
        se = se[inter_idx],
        interactions_location = interactions_location,
        interactions_scale = interactions_scale,
        ess = ess,  # effective sample size of the posterior draws of each parameter
        M = nsim,
        zero_beyond_draws = zero_beyond_draws,
        bf_floor = bf_floor
    )

    # store posterior draws
    pmles$draws <- out$draws

    # --- Multivariate effective sample size of the posterior draws (Vats, Flegal and Jones, 2019) ---
    pmles$mess <- .mess_mcmc(out$draws)

    return(structure(pmles, class = c("dmrfit_bayes", "dmrfit")))
}

# summary.dmrfit_bayes
#' @title Summary of a \code{dmrfit_bayes} object
#' @rdname summary.dmrfit_bayes
#' @description Produces a detailed summary of a discrete MRF model fitted via Bayesian MCMC sampling.
#' @param object a \code{dmrfit_bayes} object.
#' @param ... further arguments (currently unused).
#' @method summary dmrfit_bayes
#'
#' @return an object of class \code{summary.dmrfit_bayes}.
#'
#' @export
#'
summary.dmrfit_bayes <- function(object, ...) {
    if (!inherits(object, "dmrfit_bayes"))
        stop("object is not of class 'dmrfit_bayes'")

    P <- object$P
    n_categories <- object$n_categories
    n_thresholds <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2
    n_pars <- n_thresholds + n_interactions
    pars <- object$argument

    # posterior SDs from MCMC draws
    se <- apply(object$draws, 1, sd)
    names(se) <- names(pars)

    # posterior means from MCMC draws
    post_mean <- rowMeans(object$draws)
    names(post_mean) <- names(pars)

    # 95% credible intervals (they are not posterior HDI)
    ci <- apply(object$draws, 1, quantile, probs = c(0.025, 0.975))

    # build coefficient table
    coef_table <- cbind(
        `Post.Mode`  = pars,
        `Post.Mean`  = post_mean,
        `Post.SD`    = se,
        `2.5%`       = ci[1, ],
        `97.5%`      = ci[2, ]
    )
    rownames(coef_table) <- names(pars)

    inter_idx <- (n_thresholds + 1):n_pars

    out <- list(
        call             = object$call,
        coefficients     = coef_table,
        thresholds       = coef_table[1:n_thresholds, , drop = FALSE],
        interactions     = coef_table[inter_idx, , drop = FALSE],
        neg_pseudo_loglik = object$utils$value,
        P                = P,
        N                = object$N,
        n_categories     = n_categories,
        acceptance       = object$acceptance,
        seconds_elapsed  = object$seconds_elapsed,
        method           = object$method,
        scale            = object$scale,
        nsim             = ncol(object$draws),
        mess             = object$mess,
        savage_dickey    = object$savage_dickey
    )
    class(out) <- c("summary.dmrfit_bayes", "summary.dmrfit")
    out
}

# print.dmrfit_bayes
#' @title Print a \code{dmrfit_bayes} object
#' @rdname print.dmrfit_bayes
#' @description Prints a brief overview of a discrete MRF model fitted via Bayesian MCMC sampling.
#' @param x a \code{dmrfit_bayes} object.
#' @param ... further arguments passed to \code{print}.
#' @method print dmrfit_bayes
#'
#' @return the \code{dmrfit_bayes} object \code{x}, invisibly.
#'
#' @export
#'
print.dmrfit_bayes <- function(x, ...) {
    if (!inherits(x, "dmrfit_bayes"))
        stop("object is not of class 'dmrfit_bayes'")

    cat("\nCall:\n")
    print(x$call)

    P            <- x$P
    n_categories <- x$n_categories
    n_thresholds <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2
    n_pars <- n_thresholds + n_interactions
    pars <- x$argument

    cat("\nThresholds (posterior mode):\n")
    print(round(pars[1:n_thresholds], 4))

    inter_idx <- (n_thresholds + 1):n_pars
    cat("\nPairwise interactions (posterior mode):\n")
    print(round(pars[inter_idx], 4))

    cat("\nPosterior:", .method_label(x$method, x$scale), "\n")
    cat("MCMC samples:", ncol(x$draws),
        " Acceptance rate:", round(x$acceptance, 3),
        " Elapsed:", round(x$seconds_elapsed, 1), "sec\n")
    if (!is.null(x$mess) && !is.na(x$mess))
        cat("Multivariate effective sample size:", round(x$mess, 1), "\n")
    cat("\n")
    invisible(x)
}

# print.summary.dmrfit_bayes
#' @title Print a \code{summary.dmrfit_bayes} object
#' @rdname print.summary.dmrfit_bayes
#' @description Prints the detailed summary of a discrete MRF model fitted via Bayesian MCMC sampling.
#' @param x a \code{summary.dmrfit_bayes} object.
#' @param ... further arguments passed to \code{print}.
#' @method print summary.dmrfit_bayes
#'
#' @return the \code{summary.dmrfit_bayes} object \code{x}, invisibly.
#'
#' @export
#'
print.summary.dmrfit_bayes <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Discrete Markov Random Field")
    if (all(x$n_categories == 2)) cat(" (Ising)")
    cat("\n")
    cat("Estimation method: Bayesian (FisherMALA),", .method_label(x$method, x$scale), "\n")
    cat("Nodes:", x$P, " Observations:", x$N,
        " Parameters:", nrow(x$coefficients), "\n")
    cat("MCMC samples:", x$nsim,
        " Acceptance rate:", round(x$acceptance, 3),
        " Elapsed:", round(x$seconds_elapsed, 1), "sec\n")
    if (!is.null(x$mess) && !is.na(x$mess))
        cat("Multivariate effective sample size:", round(x$mess, 1), "\n")
    cat(paste0(rep("-", min(60, getOption("width"))), collapse = ""), "\n")

    cat("\nThresholds:\n")
    print(round(x$thresholds, 4))

    cat("\nPairwise interactions:\n")
    print(round(x$interactions, 4))

    cat("\nNegative pseudo-loglikelihood at posterior mode:", round(x$neg_pseudo_loglik, 4), "\n")

    if (!is.null(x$savage_dickey)) {
        sd <- x$savage_dickey
        cat("\nSavage-Dickey density ratio  [prior: Cauchy(", sd$interactions_location, ",",
            sd$interactions_scale, ")]\n")
        cat("H0: sigma = 0 for each pairwise interaction\n\n")
        tbl <- data.frame(
            `Post.Mode` = round(sd$estimate, 4),
            `Post.SD`   = round(sd$se, 4),
            `BF_01`     = formatC(sd$bf_01, format = "g", digits = 4),
            `Pr(=0|x)`  = formatC(sd$pr_null, format = "g", digits = 4),
            check.names = FALSE
        )
        floored <- if (is.null(sd$zero_beyond_draws)) rep(FALSE, nrow(tbl)) else sd$zero_beyond_draws
        tbl$BF_01[floored] <- paste0("< ", formatC(sd$bf_floor, format = "g", digits = 2))
        tbl$`Pr(=0|x)`[floored] <- paste0("< ", formatC(sd$bf_floor, format = "g", digits = 2))
        rownames(tbl) <- names(sd$estimate)
        print(tbl)
        if (any(floored))
            cat("\n'<': zero lies beyond all posterior draws, so BF_01 is only bounded (very strong evidence for an",
                "interaction).\n")
        ess_inter <- sd$ess[names(sd$estimate)]
        if (length(ess_inter) && all(!is.na(ess_inter))) {
            cat("\nMCMC samples:", sd$M, "  Effective sample size of the interactions: min",
                round(min(ess_inter), 1), ", median", round(stats::median(ess_inter), 1), "\n")
        } else {
            cat("\nMCMC samples:", sd$M, "\n")
        }
    }

    invisible(x)
}

#' ess_mcmc (internal)
#' @description Effective sample size of a single chain of draws, from its autocorrelations summed with Geyer's
#'   initial positive sequence (sums of pairs of consecutive autocorrelations, up to the first non-positive pair).
#'   The autocorrelations are computed through the fast Fourier transform.
#' @param x numeric vector, the draws of one parameter
#' @return the effective sample size (NA for fewer than 4 draws or constant draws)
#' @noRd
.ess_mcmc <- function(x) {

    n <- length(x)
    if (n < 4 || stats::var(x) == 0) return(NA_real_)

    # --- Autocorrelations at lags 0, ..., n - 1 (zero-padded to avoid wrap-around) ---
    xc <- x - mean(x)
    m <- 2^ceiling(log2(2 * n))
    f <- stats::fft(c(xc, rep(0, m - n)))
    acov <- Re(stats::fft(Mod(f)^2, inverse = TRUE))[seq_len(n)] / m
    rho <- acov / acov[1]

    # --- Integrated autocorrelation time (Geyer's initial positive sequence) ---
    tau <- -1
    for (k in seq(1, n - 1, by = 2)) {
        pair <- rho[k] + rho[k + 1]
        if (is.na(pair) || pair <= 0) break
        tau <- tau + 2 * pair
    }

    # --- Return the effective sample size ---
    return(n / max(tau, 1 / n))
}

#' mess_mcmc (internal)
#' @description Multivariate effective sample size of a chain of draws (Vats, Flegal and Jones, 2019, Biometrika),
#'   n * (det(Lambda) / det(Sigma))^(1/p), with Lambda the sample covariance of the draws and Sigma the asymptotic
#'   covariance of their mean. Sigma includes the autocorrelation of the chain, so the ratio shrinks when the chain
#'   mixes slowly. Sigma is estimated by lugsail batch means (Vats and Flegal, 2022), 2 * Sigma_b - Sigma_{b/3} with
#'   batch size b = floor(sqrt(n)), which corrects the underestimation of plain batch means for slowly mixing chains;
#'   plain batch means are used when the lugsail estimate is not positive definite.
#' @param draws matrix of size [p x n] with the draws (parameters by rows)
#' @return the multivariate effective sample size (NA when there are not more batches than parameters, i.e. Sigma is
#'   singular)
#' @noRd
.mess_mcmc <- function(draws) {

    p <- nrow(draws)
    n <- ncol(draws)
    b <- floor(sqrt(n))
    if (floor(n / b) <= p || floor(b / 3) < 1) return(NA_real_)
    x <- t(draws)

    # --- Batch-means estimate of Sigma with batch size b ---
    sigma_bm <- function(b) {
        a <- floor(n / b)
        batch_means <- rowsum(x[seq_len(a * b), , drop = FALSE], rep(seq_len(a), each = b), reorder = FALSE) / b
        b * stats::cov(batch_means)
    }

    # --- Log-determinant (NA if the matrix is not positive definite) ---
    log_det <- function(m) {
        d <- determinant(m, logarithm = TRUE)
        if (d$sign <= 0) NA_real_ else as.numeric(d$modulus)
    }

    # --- Lugsail estimate of Sigma, plain batch means as fallback ---
    Sigma_b <- sigma_bm(b)
    Sigma <- 2 * Sigma_b - sigma_bm(floor(b / 3))
    if (is.na(log_det(Sigma)) || any(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values <= 0)) Sigma <- Sigma_b

    # --- Return the multivariate effective sample size ---
    ratio <- (log_det(stats::cov(x)) - log_det(Sigma)) / p
    return(if (is.na(ratio)) NA_real_ else n * exp(ratio))
}

#' method_label (internal)
#' @description Label of the posterior targeted by a dmrfit_bayes fit, for print and summary.
#' @param method the sampler ("core", "adacore", "exact" or "dmh"; NULL for fits made before the method argument)
#' @param scale the rescaling covariance of method = "core" ("ghw", "mch" or "rm")
#' @return a character string
#' @noRd
.method_label <- function(method, scale) {
    if (is.null(method)) return("coordinate-rescaled pseudo-posterior (sandwich covariance)")
    switch(method,
           core = paste0("coordinate-rescaled pseudo-posterior (",
                         switch(scale, ghw = "sandwich covariance", mch = "Monte Carlo Hessian at the pseudo-posterior mode",
                                rm = "Monte Carlo Hessian at the Robbins-Monro estimate of the full-posterior mode"), ")"),
           adacore = "coordinate-rescaled pseudo-posterior (sandwich covariance adapted during burn-in)",
           exact = "full-likelihood posterior (exact normalizing constant)",
           dmh = "full-likelihood posterior (double Metropolis-Hastings)")
}
