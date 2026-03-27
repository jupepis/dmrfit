#' @title dmrfit_bayes
#' 
#' @description Fit a discrete Markov Random Field model via pseudo-likelihood estimation. The optimization is performed using the Coordinate Rescaling sampler implemented by Arena and Marsman (2026). 
#' 
#' @param data data matrix, with rows as samples and columns as variables. Each variable should be rescaled to the range of 0 to m-1, where m is the number of categories for that variable. The baseline category is always the minimum value in the variable. The internal processing will check if the variables are rescaled and will rescale them if necessary. If there are any NAs in the data, they will be removed before optimization (listwise deletion).
#' @param parinit parinit, initial parameter estimates, a vector of length equal to the number of parameters in the model. If NULL, it will be initialized to a vector of zeros. The number of parameters in the model is calculated as sum(P * (n_categories - 1)) + P * (P - 1) / 2, where P is the number of variables.
#' @param nsim nsim, number of iterations for the optimization. Default is 1000.
#' @param burnin burnin, number of burnin iterations for the optimization. Default is 1000.
#' @param ncores ncores, number of cores to use for parallel computing of the gradient, hessian and likelihood values. If ncores is larger than the number of available cores, it will be set to the number of available cores minus one. Default is 1.
#' @param thresholds_alpha alpha parameter for the Beta-Prime prior on thresholds (default is 0.5). 
#' @param thresholds_beta beta parameter for the Beta-Prime prior on thresholds (default is 0.5). 
#' @param interactions_location location parameter for the Cauchy prior on pairwise interactions (default is 0.0). 
#' @param interactions_scale scale parameter for the Cauchy prior on pairwise interactions (default is 2.5). 
#' @param sigma2 initial value for the adaptive variance parameter in the FisherMALA sampler (default is 0.1).
#' @param seed random seed for reproducibility of the BSIR step (default is 123).
#' 
#' @return dmrfit S3 class object
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
#' @export 
#' 
dmrfit_bayes <- function(data, parinit = NULL, nsim = 1e03, burnin = 1e03, ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5, interactions_location = 0.0, interactions_scale = 2.5, sigma2 = 0.1, seed = 123) {

    # save the matched call for print and summary methods
    cl <- match.call()

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

    # if parinit is not provided, initialize it to a vector of zeros
    if(is.null(parinit)) {
        parinit <- rep(0.0, n_pars)
    }
    else if(length(parinit) != n_pars) {
        stop(paste("Length of parinit must be equal to the number of parameters in the model:", n_pars))
    }

    # if n_cores is less than 1, set it to 1 or the number of available cores minus one, or the number of cores specified by the user, whichever is smallest
    ncores <- parallel::detectCores()-1
    

    # cross-product terms for the pairwise associations
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    # finding the PMLEs via optimization of the pseudo-likelihood with trust region method
    pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = TRUE, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    if(is.null(pmles)) {
        warning("Optimization failed. Returning NULL.")
        return(NULL)
    }

    # transformation matrices
    current_sale <- chol(pmles$utils$hessian)
    new_scale <- t(chol(pmles$utils$HW))
    
    # run the core sampler 
    out <- suppressWarnings(tryCatch(expr = dmrfit:::omrf_core_sampler(
        data = t(data[,1:P, drop = FALSE]), # only the original data (without the cross-product terms) is needed for the core sampler, which computes the pseudo-likelihood and its gradient
        pars = pmles$argument,
        n_categories = n_categories,
        P = P,
        nsim = nsim,
        burnin = burnin,
        pmles = pmles$argument,
        current_scale = current_sale, # cholesky of hessian matrix of pseudo-likelihood calculate at the 'pmles' 
        new_scale = new_scale, # transposed cholesky of Godambe-Huber-White sandwich estimator
        adaptive_stage_n_iter = 500, # iterations for the initial adaptive stage (simple MALA), fixed to 500 for now
        sigma2 = sigma2, # for ordinal MRF, set this value to [[0.01]] or higher 0.1, by default this is set to 1.0 but with the ordinal MRF the algorithms runs fine with 0.1 or 0.01 (this sigma2 is adaptive and a good starting value is needed, however, it is allowed to vary over iterations)
        thresholds_alpha = thresholds_alpha,
        thresholds_beta = thresholds_beta,
        interactions_location = interactions_location,
        interactions_scale = interactions_scale
    ), error = function(e) {NULL}))
    if(is.null(out)) {
        warning("Core sampler failed. Returning NULL.")
        return(NULL)
    }

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

    # label the parameter vector
    n_thresholds <- sum(n_categories - 1)
    thresh_names <- unlist(lapply(seq_len(P), function(p) {
        paste0("mu[", p, ",", seq_len(n_categories[p] - 1), "]")
    }))
    inter_names <- unlist(lapply(1:(P - 1), function(j) {
        lapply((j + 1):P, function(i) paste0("sigma[", i, ",", j, "]"))
    }))
    par_names <- c(thresh_names, inter_names)
    names(pmles$argument) <- par_names

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

    for (k in seq_along(inter_idx)) {
        j <- inter_idx[k]
        draws_j <- out$draws[j, ]
        bf_01[k] <- tryCatch({
            range_j <- range(draws_j)
            d <- density(draws_j, n = 1024, from = range_j[1], to = range_j[2])
            log_post_at_zero <- log(approx(x = d$x, y = d$y, xout = 0.0)$y)
            log_post_at_zero <- ifelse(is.na(log_post_at_zero), log(.Machine$double.eps), log_post_at_zero)
            log_post_at_zero <- ifelse(log_post_at_zero == Inf, -log(.Machine$double.eps), log_post_at_zero)
            exp(log_post_at_zero - log_prior_at_zero)
        }, error = function(e) {
            warning("Savage-Dickey: density estimation failed for ", par_names[j],
                    ". The 'BF_01' is set to NA.")
            NA_real_
        })
    }

    pr_null <- bf_01 / (1 + bf_01) # this is the Pr(=0|x)

    # ESS via determinant ratio
    draws_cov <- cov(t(out$draws))  # n_pars x n_pars covariance of posterior draws
    Sigma <- pmles$utils$HW         # reference covariance (Huber-White)

    logdetX <- as.numeric(determinant(draws_cov, logarithm = TRUE)$modulus)
    logdetZ <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
    ess <- nsim * exp((logdetX - logdetZ) / n_pars)

    pmles$savage_dickey <- list(
        bf_01 = bf_01,
        pr_null = pr_null,
        estimate = pars[inter_idx],
        se = se[inter_idx],
        interactions_location = interactions_location,
        interactions_scale = interactions_scale,
        ess = ess,  # rough ESS approximation
        M = nsim
    )

    # store posterior draws
    pmles$draws <- out$draws

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
        nsim             = ncol(object$draws),
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

    cat("\nMCMC samples:", ncol(x$draws),
        " Acceptance rate:", round(x$acceptance, 3),
        " Elapsed:", round(x$seconds_elapsed, 1), "sec\n")
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
    cat("Estimation method: Bayesian (FisherMALA)\n")
    cat("Nodes:", x$P, " Observations:", x$N,
        " Parameters:", nrow(x$coefficients), "\n")
    cat("MCMC samples:", x$nsim,
        " Acceptance rate:", round(x$acceptance, 3),
        " Elapsed:", round(x$seconds_elapsed, 1), "sec\n")
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
        tbl <- cbind(
            `Post.Mode` = sd$estimate,
            `Post.SD`   = sd$se,
            `BF_01`     = sd$bf_01,
            `Pr(=0|x)`  = sd$pr_null
        )
        print(round(tbl, 4))
        cat("\nMCMC samples:", sd$M,
            "  Effective sample size:", round(sd$ess, 1), "\n")
    }

    invisible(x)
}