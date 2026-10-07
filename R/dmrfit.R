#' @title Point estimation and Bayes factors for discrete Markov random fields
#'
#' @description Fits a discrete Markov random field (Ising or ordinal) through the pseudo-likelihood. The estimates
#' maximize the pseudo-likelihood, or the pseudo-posterior with \code{with_prior = TRUE}, and are obtained with a
#' trust region algorithm; the standard errors come from the Godambe-Huber-White (GHW, or sandwich) estimator. Optionally, Savage-Dickey
#' Bayes factors are computed for each pairwise interaction from the coordinate-rescaled pseudo-posterior (Arena and
#' Marsman, 2026).
#'
#' @param data data matrix, with rows as samples and columns as variables. Each variable should be rescaled to the range of 0 to m-1, where m is the number of categories for that variable. The baseline category is always the minimum value in the variable. The internal processing will check if the variables are rescaled and will rescale them if necessary. If there are any NAs in the data, they will be removed before optimization (listwise deletion).
#' @param parinit initial parameter values for the optimization, a vector of length equal to the number of parameters in the model, \code{sum(n_categories - 1) + P * (P - 1) / 2}, where \code{P} is the number of variables. If NULL (default), a vector of zeros.
#' @param structure network structure, a P x P symmetric matrix with 1 for an edge and 0 for no edge. If NULL (default), the network is fully connected.
#' @param with_prior logical, whether to include the prior in the optimization: a Beta-Prime prior on the thresholds and a Cauchy prior on the pairwise interactions. Default is FALSE.
#' @param savage_dickey logical, whether to compute the Savage-Dickey density ratio Bayes factor for each pairwise interaction via Bayesian Sampling Importance Resampling (BSIR) from the coordinate-rescaled pseudo-posterior: the pseudo-posterior rescaled around its mode to the Godambe-Huber-White (GHW) covariance, which corrects the underestimated posterior variability of the pseudo-likelihood. The proposal is a multivariate t with the GHW covariance. When zero lies beyond all resampled draws of an interaction, its density at zero cannot be estimated from the draws: the Bayes factor is then reported at the floor \code{1/(10 * M)} and flagged in \code{savage_dickey$zero_beyond_draws}. Only available when \code{with_prior = TRUE}. Default is FALSE.
#' @param M number of draws resampled in the BSIR step (default is 1000). Only used when \code{savage_dickey = TRUE}.
#' @param oversampling multiplier for the number of proposal draws in the BSIR step (default is 10). The total number of proposal draws is \code{M * oversampling}. Only used when \code{savage_dickey = TRUE}.
#' @param ncores number of cores used to compute the gradient, Hessian and pseudo-likelihood in parallel. It is capped at the number of available cores minus one, and has no effect when the package was built without OpenMP support. Default is 1.
#' @param thresholds_alpha alpha parameter for the Beta-Prime prior on thresholds (default is 0.5). Only used when \code{with_prior = TRUE}.
#' @param thresholds_beta beta parameter for the Beta-Prime prior on thresholds (default is 0.5). Only used when \code{with_prior = TRUE}.
#' @param interactions_location location parameter for the Cauchy prior on pairwise interactions (default is 0.0). Only used when \code{with_prior = TRUE}.
#' @param interactions_scale scale parameter for the Cauchy prior on pairwise interactions (default is 2.5). Only used when \code{with_prior = TRUE}.
#' @param proposal_df degrees of freedom for the multivariate t proposal distribution in the BSIR step (default is 5). Only used when \code{savage_dickey = TRUE}.
#' @param lrt_intervals whether to compute likelihood-ratio intervals (default is FALSE; see Details): \code{TRUE}
#'   for all free parameters, or a vector of parameter names (as in \code{names(fit$argument)}) for some of them. Each
#'   interval needs a series of constrained optimizations, which take a while for large networks. They are returned
#'   by \code{confint(fit, method = "lrt")}. When they are not computed and there are fewer than 10 observations per
#'   free parameter, a warning suggests them.
#' @param level confidence level of the intervals (default is 0.95): the Wald intervals printed by \code{summary()}
#'   and the likelihood-ratio intervals, so that both are at the same level. It is also the default level of
#'   \code{confint()}.
#' @param seed random seed for reproducibility of the BSIR step (default is 123). The caller's random number stream is restored on exit. Only used when \code{savage_dickey = TRUE}.
#'
#' @details The likelihood-ratio interval of a parameter \eqn{\theta_k} is based on its profile pseudo-likelihood:
#' for each value of \eqn{\theta_k}, the other parameters are set to the values that maximize the pseudo-likelihood
#' (or the pseudo-posterior with \code{with_prior = TRUE}). The interval collects the values where
#' \deqn{C_k \cdot 2 [\log PL(\hat\theta) - \max_{\theta_{-k}} \log PL(\theta_k, \theta_{-k})] \le \chi^2_{1, level}.}
#' For a full likelihood the statistic would be \eqn{\chi^2_1}-distributed (\eqn{C_k = 1}); for the
#' pseudo-likelihood it is \eqn{\lambda_k \chi^2_1}-distributed, with \eqn{\lambda_k = \Sigma_{kk} / (H^{-1})_{kk}},
#' because the curvature \eqn{H} of the pseudo-likelihood overstates the information in the data
#' (\eqn{\Sigma} is the GHW covariance). The scaling \eqn{C_k = 1 / \lambda_k} calibrates the statistic (Pace,
#' Salvan and Sartori, 2011; Varin, Reid and Firth, 2011). Unlike the Wald intervals, the likelihood-ratio intervals
#' can be asymmetric around the estimate. The bounds are found by \code{uniroot} on each side of the estimate. The
#' asymmetry matters when the sample size is small relative to the number of parameters (Arena and Marsman, 2026), so
#' with fewer than 10 observations per free parameter \code{dmrfit()} warns that the Wald intervals may be inaccurate.
#'
#' @return an object of class \code{dmrfit}, a list including the estimates (\code{argument}), the value, gradient,
#'   Hessian and GHW covariance of the objective at the estimates (\code{utils}), the matched call
#'   (\code{call}), the data dimensions (\code{P}, \code{N}, \code{n_categories}), and, with
#'   \code{savage_dickey = TRUE}, the Savage-Dickey Bayes factors with the importance-sampling effective sample
#'   size (\code{savage_dickey}).
#'
#' @references Arena, G. and Marsman, M. (2026). Bayesian inference for discrete Markov random fields through
#' coordinate rescaling. arXiv preprint. \doi{10.48550/arXiv.2601.17205}
#'
#' Pace, L., Salvan, A., and Sartori, N. (2011). Adjusting composite likelihood ratio statistics. \emph{Statistica
#' Sinica}, 21(1), 129-148.
#'
#' Skare, Ø., Bølviken, E., and Holden, L. (2003). Improved sampling-importance resampling and reduced bias importance
#' sampling. \emph{Scandinavian Journal of Statistics}, 30(4), 719-737.
#'
#' Varin, C., Reid, N., and Firth, D. (2011). An overview of composite likelihood methods. \emph{Statistica Sinica},
#' 21(1), 5-42.
#'
#' @seealso \code{\link{dmrfit_bayes}} for posterior sampling.
#'
#' @examples
#'
#' # binary responses of 1000 observations on 3 nodes (independent, for illustration)
#' set.seed(123)
#' n <- 1000
#' P <- 3
#' prob <- c(0.2, 0.5, 0.3)[rep(1:P, each = n)]
#' data <- matrix(rbinom(n * P, size = 1, prob = prob), nrow = n, ncol = P)
#' fit <- dmrfit(data, with_prior = TRUE, savage_dickey = TRUE)
#' fit
#' summary(fit)
#'
#' # constrained network: no edge between nodes 1 and 3
#' structure <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = P)
#' fit_constrained <- dmrfit(data, structure = structure)
#' summary(fit_constrained)
#'
#' @export
#'
dmrfit <- function(data, parinit = NULL, structure = NULL, with_prior = FALSE, savage_dickey = FALSE, M = 1000, oversampling = 10, ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5, interactions_location = 0.0, interactions_scale = 2.5, proposal_df = 5, lrt_intervals = FALSE, level = 0.95, seed = 123) {

    # save the matched call for print and summary methods
    cl <- match.call()

    if (!(is.character(lrt_intervals) || (is.logical(lrt_intervals) && length(lrt_intervals) == 1 && !is.na(lrt_intervals)))) {
        stop("lrt_intervals must be TRUE, FALSE or a vector of parameter names.")
    }
    if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
        stop("level must be a number between 0 and 1.")
    }

    # validate savage_dickey requires with_prior
    if (savage_dickey && !with_prior)
        stop("savage_dickey = TRUE requires with_prior = TRUE.")
    
    # --- Random seed for the BSIR step (the caller's random state is restored on exit) ---
    if(savage_dickey){
        had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit({
            if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
        }, add = TRUE)
        set.seed(seed)
    }

    # processing input arguments

    # check if data is a matrix, if not convert it to a matrix
    if(!is.matrix(data)) {
        data <- as.matrix(data)
    }

     # remove NAs from data if any exist
    data <- data[!is.na(rowSums(data)), , drop = FALSE]

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

    # if parinit is not provided, initialize it to a vector of zeros
    if(is.null(parinit)) {
        parinit <- rep(0.0, n_pars)
    }
    else if(length(parinit) != n_pars) {
        stop(paste("Length of parinit must be equal to the number of parameters in the model:", n_pars))
    }

    # if structure is provided, check if it is a square matrix with dimensions equal to the number of variables
    if(!is.null(structure)) {
        if(!is.matrix(structure) || nrow(structure) != P || ncol(structure) != P) {
                stop("structure must be a square matrix with dimensions equal to the number of variables (columns) in data.")
        }
        if(!isSymmetric(structure)) {
            stop("structure must be a symmetric matrix.")
        }
        if(any(structure[lower.tri(structure, diag = FALSE)] != 0 & structure[lower.tri(structure, diag = FALSE)] != 1)) {
            stop("structure must contain only 0s and 1s.")
        }   
    }

    # number of cores: the value given by the user, capped at the number of available cores minus one
    ncores <- max(1L, min(as.integer(ncores), parallel::detectCores() - 1L), na.rm = TRUE)

    # cross-product terms for the pairwise associations
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    if(is.null(structure)) {
        pmles <- suppressWarnings(tryCatch(expr = cpp_optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    } else {
        structure_input_optimize <- c(rep(1, n_thresholds), structure[lower.tri(structure, diag = FALSE)])
        pmles <- suppressWarnings(tryCatch(expr = cpp_optimize_with_structure(data = data, parinit = parinit, n_categories =  n_categories, P = P, structure = structure_input_optimize, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps) , n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    }

    if(is.null(pmles)) {
        warning("Optimization failed. Returning NULL.")
        return(NULL)
    }

    # metadata needed by print/summary
    pmles$call <- cl
    pmles$P <- P
    pmles$var_names <- if (is.null(colnames(data))) paste0("V", seq_len(P)) else colnames(data)[seq_len(P)] # node labels for plot()
    pmles$n_categories <- n_categories
    pmles$N <- nrow(data)
    pmles$with_prior <- with_prior
    pmles$structured <- !is.null(structure)
    pmles$structure <- structure
    pmles$ncores <- ncores
    pmles$level <- level # level of the Wald and likelihood-ratio intervals

    # label the parameter vector
    n_thresholds <- sum(n_categories - 1)
    thresh_names <- unlist(lapply(seq_len(P), function(p) {
        paste0("mu[", p, ",", seq_len(n_categories[p] - 1), "]")
    }))
    inter_names <- unlist(lapply(1:(P - 1), function(j) {
        lapply((j + 1):P, function(i) paste0("theta[", i, ",", j, "]"))
    }))
    names(pmles$argument) <- c(thresh_names, inter_names)

    # free parameters: all thresholds, and the interactions present in the structure
    free_idx <- if (is.null(structure)) seq_len(n_pars) else which(c(rep(TRUE, n_thresholds), structure[lower.tri(structure)] == 1))

    # --- GHW covariance of a constrained fit: from the free parameters only (the absent edges are not estimated),
    # with the same construction as dmrf_deriv(): the GHW covariance of the pseudo-likelihood, H^-1 J H^-1, and with a prior
    # the inverse of its inverse plus the prior precision ---
    if (!is.null(structure)) {
        lik <- if (with_prior) {
            dmrf_deriv(pars = pmles$argument, data = data, P = P, n_categories = n_categories, with_prior = FALSE,
                       ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                       interactions_location = interactions_location, interactions_scale = interactions_scale)
        } else pmles$utils
        H_free_inv <- solve(lik$hessian[free_idx, free_idx, drop = FALSE])
        HW_free <- H_free_inv %*% lik$FisherInfo[free_idx, free_idx, drop = FALSE] %*% H_free_inv
        if (with_prior) {
            prior_precision <- (pmles$utils$hessian - lik$hessian)[free_idx, free_idx, drop = FALSE] # diagonal
            HW_free <- solve(solve(HW_free) + prior_precision)
        }
        HW <- matrix(0, n_pars, n_pars)
        HW[free_idx, free_idx] <- HW_free
        pmles$utils$HW <- HW
    }

    # --- few observations per parameter: the Wald intervals assume a symmetric sampling distribution ---
    if (isFALSE(lrt_intervals) && nrow(data) < 10 * length(free_idx)) {
        warning("The sample size (", nrow(data), ") is small relative to the number of free parameters (", length(free_idx),
                "): the Wald intervals assume a symmetric sampling distribution and may be inaccurate. Likelihood-ratio ",
                "intervals are suggested: refit with lrt_intervals = TRUE (or a vector of parameter names).")
    }

    # --- likelihood-ratio intervals of the free parameters ---
    if (is.character(lrt_intervals) || isTRUE(lrt_intervals)) {
        lrt_idx <- if (is.character(lrt_intervals)) match(lrt_intervals, names(pmles$argument)) else free_idx
        if (anyNA(lrt_idx) || !all(lrt_idx %in% free_idx)) {
            stop("lrt_intervals: unknown parameter(s) or absent edges of the constrained fit: ",
                 paste(lrt_intervals[is.na(lrt_idx) | !(lrt_idx %in% free_idx)], collapse = ", "), ".")
        }
        pmles$lrt_intervals <- .lrt_intervals(pars = pmles$argument, nll_hat = pmles$utils$value, lrt_idx = lrt_idx,
                                              data = data, free_idx = free_idx, H = pmles$utils$hessian,
                                              Sigma = pmles$utils$HW, level = level, P = P, n_categories = n_categories,
                                              with_prior = with_prior, ncores = ncores, thresholds_alpha = thresholds_alpha,
                                              thresholds_beta = thresholds_beta, interactions_location = interactions_location,
                                              interactions_scale = interactions_scale)
    }

    # --- Savage-Dickey via BSIR ---
    if(savage_dickey){

        pars <- pmles$argument
        Sigma <- pmles$utils$HW
        se <- sqrt(diag(Sigma))

        # identify free parameters (all thresholds + interactions where structure == 1)
        if(!is.null(structure)){
            free_inter <- structure[lower.tri(structure, diag = FALSE)] == 1
            free_idx <- c(1:n_pars)[c(rep(TRUE, n_thresholds), free_inter)]
        } else {
            free_inter <- rep(TRUE, n_interactions)
            free_idx <- seq_len(n_pars)
        }
        n_pars_free <- length(free_idx)

        if (M <= n_pars_free)
            warning("M (", M, ") is not larger than the number of free parameters (", n_pars_free, "). Issues with SIR covariance and ESS estimation may arise. Consider increasing M.")

        inter_idx <- (n_thresholds + 1):n_pars
        free_inter_idx <- inter_idx[free_inter]

        # proposal only over free parameters (this is necessary if there are fixed interaction parameters = 0)
        pars_free  <- pars[free_idx]
        Sigma_free <- Sigma[free_idx, free_idx, drop = FALSE]

        # draw M samples from proposal
        M_importance <- M * oversampling
        Z <- cpp_mvnrnd_arma(mu = rep(0, n_pars_free), Sigma = Sigma_free, n = M_importance)
        V <- rchisq(n = M_importance, df = proposal_df)
        samples_free <- sweep(Z, 2, sqrt(proposal_df/V), "*")
        samples_free <- sweep(samples_free, 1, pars_free, "+")  # n_pars_free x M_importance

        # reconstruct full parameter vectors (fixed params stay at 0)
        samples <- matrix(0, n_pars, M_importance)
        samples[free_idx, ] <- samples_free

        # log-density of the proposal at each sample (multivariate t, free subspace)
        L <- chol(Sigma_free)
        log_det_Sigma <- 2 * sum(log(diag(L)))
        diff_mat <- sweep(samples_free, 1, pars_free, check.margin = FALSE)
        solve_L <- backsolve(L, diff_mat, transpose = TRUE)
        mahal <- colSums(solve_L^2)
        log_q <- lgamma((proposal_df + n_pars_free)/2) - lgamma(proposal_df/2) - (n_pars_free/2) * log(proposal_df*pi) - 0.5*log_det_Sigma - ((proposal_df + n_pars_free)/2) * log(1 + mahal/proposal_df)

        # --- Target: the coordinate-rescaled (CoRe) pseudo-posterior ---
        # Its covariance is the GHW covariance Sigma (the covariance of the proposal). Each proposal draw beta is
        # mapped back to the pseudo-posterior scale,
        # eta = A^{-1} (beta - pars) + pars with A^{-1} = R^{-1} L^{-T} (R'R = H, the curvature of the negative log
        # pseudo-posterior at the mode; L'L = Sigma), and the pseudo-posterior is evaluated there. The Jacobian of the
        # map is constant and cancels when the weights are normalized. Proposal and target thus share the same scale and
        # the weights only correct for the differences in shape (a proposal with the GHW covariance against the
        # pseudo-posterior itself is too wide in every direction, and its weights degenerate as the dimension grows).
        H <- pmles$utils$hessian
        H_free <- if (nrow(H) == n_pars) H[free_idx, free_idx, drop = FALSE] else H
        R_H <- chol((H_free + t(H_free)) / 2)
        eta_free <- backsolve(R_H, solve_L) + pars_free # solve_L = L^{-T} (beta - pars), computed for log_q above
        eta <- matrix(0, n_pars, M_importance)
        eta[free_idx, ] <- eta_free

        # --- Log pseudo-posterior at the mapped draws (all draws at once, over the unique response patterns) ---
        log_target <- -cpp_npseudologlik_draws(
            pars_draws = eta,
            data = data[, 1:P, drop = FALSE],
            n_categories = n_categories,
            with_prior = TRUE,
            thresholds_alpha = thresholds_alpha,
            thresholds_beta = thresholds_beta,
            interactions_location = interactions_location,
            interactions_scale = interactions_scale
        )

        # --- Importance weights on the log scale, shifted by their maximum (the largest weight is 1) ---
        log_w <- log_target - log_q
        log_w[!is.finite(log_w)] <- -Inf
        if (all(log_w == -Inf))
            stop("Savage-Dickey: all importance weights are zero (non-finite pseudo-posterior at every proposal draw).")
        w <- exp(log_w - max(log_w))
        w_sum <- sum(w)
        ess <- w_sum^2 / sum(w^2) # effective sample size of the importance sample (out of M * oversampling draws)
        if (ess < M / 10)
            warning("Savage-Dickey: the effective sample size of the importance sample (", round(ess, 1), ") is less than ",
                    "a tenth of M (", M, "); the Bayes factors may be unreliable. Consider increasing 'oversampling'.")

        # --- ISIR weights w_i / sum_{j != i} w_j (Skare et al., 2003, Scandinavian Journal of Statistics) ---
        # The denominator is bounded away from zero when a single draw carries almost all the weight
        w_isir <- w / pmax(w_sum - w, .Machine$double.xmin)

        # --- Resample the proposal draws (on the CoRe scale) according to the importance weights ---
        idx <- sample(x = seq_len(M_importance), size = M, replace = TRUE, prob = w_isir)
        sir_samples <- samples[, idx, drop = FALSE]

        log_prior_at_zero <- dcauchy(0, location = interactions_location, scale = interactions_scale, log = TRUE)

        bf_01 <- numeric(length(free_inter_idx))
        names(bf_01) <- names(pars)[free_inter_idx]

        # --- Savage-Dickey density ratio at zero for each free interaction ---
        # When zero lies beyond all resampled draws, the density at zero cannot be estimated from the draws (a kernel
        # estimate there only extrapolates the tail of the nearest kernel): BF_01 is then reported at a floor and flagged
        bf_floor <- 1 / (10 * M)
        zero_beyond_draws <- logical(length(free_inter_idx))
        names(zero_beyond_draws) <- names(bf_01)

        for (k in seq_along(free_inter_idx)){
            j <- free_inter_idx[k]
            sir_samples_j <- sir_samples[j, ]
            zero_beyond_draws[k] <- 0 < min(sir_samples_j) || 0 > max(sir_samples_j)
            if (zero_beyond_draws[k]) { bf_01[k] <- bf_floor; next }
            bf_01[k] <- tryCatch({
                # Gaussian kernel density estimate of the marginal posterior at zero, on the log scale (log-sum-exp)
                h <- stats::bw.nrd0(sir_samples_j)
                log_k <- stats::dnorm(0, mean = sir_samples_j, sd = h, log = TRUE)
                log_post_at_zero <- max(log_k) + log(mean(exp(log_k - max(log_k))))
                exp(log_post_at_zero - log_prior_at_zero)
            }, error = function(e) {
                warning("Savage-Dickey: density estimation failed for ", names(pars)[j],
                        ". The 'BF_01' is set to NA.")
                NA_real_
            })
        }

        pr_null <- bf_01 / (1 + bf_01) # This is the Pr(=0|x)

        pmles$savage_dickey <- list(
            bf_01 = bf_01,
            pr_null = pr_null,
            estimate = pars[free_inter_idx],
            se = se[free_inter_idx],
            interactions_location = interactions_location,
            interactions_scale = interactions_scale,
            ess = ess,
            M = M,
            n_proposals = M_importance,
            zero_beyond_draws = zero_beyond_draws,
            bf_floor = bf_floor
        )
    }

    return(structure(pmles, class = "dmrfit"))
}

# print.dmrfit
#' @title Print a \code{dmrfit} object
#' @rdname print.dmrfit
#' @description Prints a brief overview of a fitted discrete MRF model, analogous to \code{print.lm}.
#' @param x a \code{dmrfit} object.
#' @param ... further arguments passed to \code{print}.
#' @method print dmrfit
#'
#' @return the \code{dmrfit} object \code{x}, invisibly.
#'
#' @export
#'
print.dmrfit <- function(x, ...) {
    if (!inherits(x, "dmrfit"))
        stop("object is not of class 'dmrfit'")

    cat("\nCall:\n")
    print(x$call)

    P  <- x$P
    n_categories <- x$n_categories
    n_thresholds   <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2

    # determine which interactions are free
    if (x$structured && !is.null(x$structure)) {
        free_inter <- x$structure[lower.tri(x$structure, diag = FALSE)] == 1
    } else {
        free_inter <- rep(TRUE, n_interactions)
    }

    pars <- x$argument

    cat("\nThresholds:\n")
    print(round(pars[1:n_thresholds], 4))

    inter_idx <- (n_thresholds + 1):(n_thresholds + n_interactions)
    free_inter_idx <- inter_idx[free_inter]
    cat("\nPairwise interactions:\n")
    print(round(pars[free_inter_idx], 4))

    cat("\n")
    invisible(x)
}

# summary.dmrfit
#' @title Summary of a \code{dmrfit} object
#' @rdname summary.dmrfit
#' @description Produces a detailed summary of a fitted discrete MRF model, analogous to
#'   \code{summary.lm}. Standard errors are obtained from the Godambe-Huber-White (GHW) estimator.
#' @param object a \code{dmrfit} object.
#' @param ... further arguments (currently unused).
#' @method summary dmrfit
#'
#' @return an object of class \code{summary.dmrfit} containing:
#'   \item{call}{the matched call.}
#'   \item{coefficients}{a matrix with columns for the estimate, standard error, z-value and p-value.}
#'   \item{thresholds}{coefficient matrix for threshold parameters.}
#'   \item{interactions}{coefficient matrix for free pairwise interaction parameters.}
#'   \item{neg_pseudo_loglik}{the negative log pseudo-likelihood at convergence.}
#'   \item{P}{number of nodes.}
#'   \item{N}{number of observations.}
#'   \item{n_categories}{vector of category counts per node.}
#'   \item{with_prior}{logical, whether a prior was used.}
#'   \item{structured}{logical, whether the network structure was constrained.}
#'   \item{intervals}{a matrix with the Wald intervals of the free parameters, and their likelihood-ratio intervals
#'     when the fit has them (\code{lrt_intervals = TRUE}), at the level of the fit.}
#'   \item{level}{the level of the intervals.}
#'
#' @export
#'
summary.dmrfit <- function(object, ...) {
    if (!inherits(object, "dmrfit"))
        stop("object is not of class 'dmrfit'")

    P <- object$P
    n_categories <- object$n_categories
    n_thresholds <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2
    n_pars <- n_thresholds + n_interactions
    pars <- object$argument

    # standard errors from Godambe-Huber-White (GHW) estimator
    se <- sqrt(diag(object$utils$HW))
    names(se) <- names(pars)

    # build full coefficient table
    z_val <- pars / se
    p_val <- 2 * pnorm(-abs(z_val))
    coef_table <- cbind(Estimate = pars, `Std. Error` = se,
                        `z value` = z_val, `Pr(>|z|)` = p_val)
    rownames(coef_table) <- names(pars)
    if(object$with_prior){
        colnames(coef_table) <- c("Post.Mode", "Post.SD", "z value", "Pr(>|z|)")
    }
    else{
        colnames(coef_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }

    # determine free interactions
    if (object$structured && !is.null(object$structure)) {
        free_inter <- object$structure[lower.tri(object$structure, diag = FALSE)] == 1
    } else {
        free_inter <- rep(TRUE, n_interactions)
    }

    inter_idx <- (n_thresholds + 1):n_pars
    free_inter_idx <- inter_idx[free_inter]

    inter_table <- coef_table[free_inter_idx, , drop = FALSE]

    out <- list(
        call = object$call,
        coefficients = coef_table[c(1:n_thresholds, free_inter_idx), , drop = FALSE],
        thresholds = coef_table[1:n_thresholds, , drop = FALSE],
        interactions = inter_table,
        neg_pseudo_loglik = object$utils$value,
        P = P,
        N = object$N,
        n_categories = n_categories,
        with_prior = object$with_prior,
        structured = object$structured,
        savage_dickey = object$savage_dickey,
        intervals = .summary_intervals(object),
        level = .fit_level(object)
    )
    class(out) <- "summary.dmrfit"
    out
}

# print.summary.dmrfit
#' @title Print a \code{summary.dmrfit} object
#' @rdname print.summary.dmrfit
#' @description Prints the detailed summary of a fitted discrete MRF model,
#'   analogous to \code{print.summary.lm}.
#' @param x a \code{summary.dmrfit} object.
#' @param ... further arguments passed to \code{printCoefmat}.
#' @method print summary.dmrfit
#'
#' @return the \code{summary.dmrfit} object \code{x}, invisibly.
#'
#' @export
#'
print.summary.dmrfit <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Discrete Markov Random Field")
    if (all(x$n_categories == 2)) cat(" (Ising)")
    cat("\n")
    cat("Estimation method: maximum pseudo-likelihood")
    if (x$with_prior) cat(" (with prior)")
    cat("\n")
    if (x$structured) cat("Network structure: constrained\n")
    cat("Nodes:", x$P, " Observations:", x$N,
        " Free parameters:", nrow(x$coefficients), "\n")
    cat("Standard errors: Godambe-Huber-White (GHW) estimator\n")
    cat(paste0(rep("-", min(60, getOption("width"))), collapse = ""), "\n")

    cat("\nThresholds:\n")
    printCoefmat(x$thresholds, P.values = TRUE, has.Pvalue = TRUE,
                 signif.stars = FALSE, ...)

    cat("\nPairwise interactions:\n")
    printCoefmat(x$interactions, P.values = TRUE, has.Pvalue = TRUE,
                 signif.stars = TRUE, ...)

    cat("\nIntervals (", format(100 * x$level), "%): Wald (GHW standard errors)",
        if (ncol(x$intervals) == 4) " and likelihood-ratio", "\n", sep = "")
    print(round(x$intervals, 4))

    cat("\nNegative log pseudo-likelihood:", round(x$neg_pseudo_loglik, 4), "\n")

    if (!is.null(x$savage_dickey)) {
        sd <- x$savage_dickey
        cat("\nSavage-Dickey density ratio  [prior: Cauchy(",sd$interactions_location,",",
            sd$interactions_scale, ")]\n")
        cat("H0: theta = 0 for each pairwise interaction\n\n")
        tbl <- data.frame(Estimate = round(sd$estimate, 4), SE = round(sd$se, 4),
                          BF_01 = formatC(sd$bf_01, format = "g", digits = 4),
                          `Pr(=0|x)` = formatC(sd$pr_null, format = "g", digits = 4),
                          check.names = FALSE)
        floored <- if (is.null(sd$zero_beyond_draws)) rep(FALSE, nrow(tbl)) else sd$zero_beyond_draws
        tbl$BF_01[floored] <- paste0("< ", formatC(sd$bf_floor, format = "g", digits = 2))
        tbl$`Pr(=0|x)`[floored] <- paste0("< ", formatC(sd$bf_floor, format = "g", digits = 2))
        rownames(tbl) <- names(sd$estimate)
        print(tbl)
        if (any(floored))
            cat("\n'<': zero lies beyond all resampled draws, so BF_01 is only bounded (very strong evidence for an",
                "interaction).\n")
        cat("\nSIR samples:", sd$M,
            "  Effective sample size of the importance sample:", round(sd$ess, 1),
            if (!is.null(sd$n_proposals)) paste0("(of ", sd$n_proposals, " proposal draws)"), "\n")
    }

    invisible(x)
}


#' lrt_intervals (internal)
#' @description Adjusted profile likelihood-ratio intervals of parameters of a dmrfit() fit. For each parameter k, the
#'   profile negative log pseudo-likelihood fixes theta_k (and the absent edges of a constrained fit at 0) and
#'   minimizes over the other free parameters (cpp_optimize_profile, warm-started at the previous solution); the bounds
#'   are the values where C_k * 2 * (profile nll - nll_hat) = qchisq(level, 1), with C_k = (H^-1)_kk / Sigma_kk and H
#'   the negative Hessian of the free parameters (Pace, Salvan and Sartori, 2011). Each bound is bracketed by steps of
#'   sqrt(Sigma_kk) away from the estimate (doubling up to 2^10 steps) and found by uniroot.
#' @param pars named vector of the estimates
#' @param nll_hat negative log pseudo-likelihood (or pseudo-posterior) at the estimates
#' @param lrt_idx indices of the parameters that get an interval
#' @param data data matrix with the cross-product columns, as passed to cpp_optimize()
#' @param free_idx indices of the free parameters
#' @param H negative Hessian of the objective at the estimates
#' @param Sigma GHW covariance at the estimates
#' @param level confidence level
#' @param P,n_categories,with_prior,ncores,thresholds_alpha,thresholds_beta,interactions_location,interactions_scale as in
#'   dmrfit()
#' @return a data frame with one row per parameter in lrt_idx (row names: the parameter names) and the columns
#'   estimate, lower, upper and C; the level is stored in attr(, "level")
#' @noRd
.lrt_intervals <- function(pars, nll_hat, lrt_idx, data, free_idx, H, Sigma, level, P, n_categories, with_prior, ncores,
                           thresholds_alpha, thresholds_beta, interactions_location, interactions_scale) {
    absent <- setdiff(seq_along(pars), free_idx) # absent edges of a constrained fit, fixed at 0
    cutoff <- stats::qchisq(level, df = 1)
    H_inv <- matrix(0, length(pars), length(pars))
    H_inv[free_idx, free_idx] <- solve(H[free_idx, free_idx, drop = FALSE])

    # --- profile negative log pseudo-likelihood at theta_k = value, warm-started at the last solution (kept in an
    # environment, so that every evaluation can update it) ---
    warm <- new.env()
    warm$start <- pars
    profile_nll <- function(k, value) {
        fit <- cpp_optimize_profile(data = data, parinit = warm$start, which_parconstr = c(absent, k) - 1,
                                    parconstr = c(rep(0, length(absent)), value), n_categories = n_categories, P = P,
                                    f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100,
                                    rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores,
                                    thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                                    interactions_location = interactions_location, interactions_scale = interactions_scale)
        warm$start <- as.vector(fit$argument)
        return(fit$utils$value)
    }

    out <- data.frame(estimate = unname(pars[lrt_idx]), lower = NA_real_, upper = NA_real_,
                      C = diag(H_inv)[lrt_idx] / diag(Sigma)[lrt_idx], row.names = names(pars)[lrt_idx])
    for (m in seq_along(lrt_idx)) {
        k <- lrt_idx[m]
        # adjusted profile statistic minus the cutoff, as a function of theta_k
        g <- function(value) out$C[m] * 2 * (profile_nll(k, value) - nll_hat) - cutoff
        step <- sqrt(Sigma[k, k])
        for (side in c(-1, 1)) {
            warm$start <- pars
            # --- bracket the root: move away from the estimate until the statistic exceeds the cutoff ---
            near <- pars[k]
            g_near <- -cutoff # the statistic is 0 at the estimate
            far <- pars[k] + side * step
            g_far <- g(far)
            doublings <- 0
            while (is.finite(g_far) && g_far < 0 && doublings < 10) {
                doublings <- doublings + 1
                near <- far
                g_near <- g_far
                far <- pars[k] + side * step * 2^doublings
                g_far <- g(far)
            }
            # the signs at the ends are those found while bracketing (the profile optimizer stops at a tolerance, so a
            # new evaluation at the same point could differ slightly)
            ends <- if (side < 0) list(x = c(far, near), f = c(g_far, g_near)) else list(x = c(near, far), f = c(g_near, g_far))
            bound <- if (is.finite(g_far) && g_far >= 0) {
                tryCatch(stats::uniroot(g, interval = ends$x, f.lower = ends$f[1], f.upper = ends$f[2], tol = 1e-6 * step)$root,
                         error = function(e) NA_real_)
            } else NA_real_
            if (side < 0) out$lower[m] <- bound else out$upper[m] <- bound
        }
    }
    if (anyNA(out[, c("lower", "upper")])) {
        warning("Likelihood-ratio intervals: some bounds could not be found (", sum(is.na(out[, c("lower", "upper")])),
                " of ", 2 * nrow(out), "); they are NA.")
    }
    attr(out, "level") <- level
    return(out)
}


#' fit_level (internal)
#' @description Level of the intervals of a dmrfit fit (0.95 for fits made before the level argument).
#' @param object a dmrfit object
#' @return a number
#' @noRd
.fit_level <- function(object) {
    return(if (is.null(object$level)) 0.95 else object$level)
}


#' summary_intervals (internal)
#' @description Wald intervals of the free parameters of a dmrfit fit, with the likelihood-ratio intervals next to
#'   them when the fit has them, at the level of the fit.
#' @param object a dmrfit object
#' @return a matrix with the columns Wald lower and upper, and LRT lower and upper
#' @noRd
.summary_intervals <- function(object) {
    level <- .fit_level(object)
    out <- confint(object, level = level, method = "wald")
    colnames(out) <- c("Wald lower", "Wald upper")
    if (!is.null(object$lrt_intervals)) {
        lrt <- confint(object, level = level, method = "lrt")
        out <- cbind(out, `LRT lower` = lrt[, 1], `LRT upper` = lrt[, 2])
    }
    return(out)
}
