#' @title dmrfit
#' 
#' @description Fit a discrete Markov Random Field model via pseudo-likelihood estimation. The optimization is performed using the trust region algorithm implemented by 
#' 
#' @param data data matrix, with rows as samples and columns as variables. Each variable should be rescaled to the range of 0 to m-1, where m is the number of categories for that variable. The baseline category is always the minimum value in the variable. The internal processing will check if the variables are rescaled and will rescale them if necessary. If there are any NAs in the data, they will be removed before optimization (listwise deletion).
#' @param parinit parinit, initial parameter estimates, a vector of length equal to the number of parameters in the model. If NULL, it will be initialized to a vector of zeros. The number of parameters in the model is calculated as sum(P * (n_categories - 1)) + P * (P - 1) / 2, where P is the number of variables.
#' @param structure network structure, P x P matrix, with 0 for no edge and 1 for edge. If NULL, then fully connected graph is assumed. Default is NULL
#' @param with_prior with_prior, logical, whether to include the prior in the optimization. If TRUE, a Beta-Prime and a Cauchy prior are applied to thresholds and pairwise associations respectively. If FALSE, no prior is applied. Default is FALSE
#' @param savage_dickey logical, whether to compute the Savage-Dickey density ratio Bayes factor for each pairwise interaction via Bayesian Sampling Importance Resampling (BSIR). Only available when \code{with_prior = TRUE}. Default is FALSE.
#' @param M number of importance samples for the BSIR step (default is 1000). Only used when \code{savage_dickey = TRUE}.
#' @param oversampling multiplier for the number of proposal draws in the BSIR step (default is 10). The total number of proposal samples is \code{M * oversampling}. Only used when \code{savage_dickey = TRUE}.
#' @param ncores ncores, number of cores to use for parallel computing of the gradient, hessian and likelihood values. If ncores is larger than the number of available cores, it will be set to the number of available cores minus one. Default is 1.
#' @param thresholds_alpha alpha parameter for the Beta-Prime prior on thresholds (default is 0.5). Only used when \code{with_prior = TRUE}.
#' @param thresholds_beta beta parameter for the Beta-Prime prior on thresholds (default is 0.5). Only used when \code{with_prior = TRUE}.
#' @param interactions_location location parameter for the Cauchy prior on pairwise interactions (default is 0.0). Only used when \code{with_prior = TRUE}.
#' @param interactions_scale scale parameter for the Cauchy prior on pairwise interactions (default is 2.5). Only used when \code{with_prior = TRUE}.
#' @param proposal_df degrees of freedom for the multivariate t proposal distribution in the BSIR step (default is 5). Only used when \code{savage_dickey = TRUE}. 
#' @param seed random seed for reproducibility of the BSIR step (default is 123). Only used when \code{savage_dickey = TRUE}.
#' 
#' @return dmrfit S3 class object
#' 
#' @example 
#' 
#' # simulate data from a 3-node Ising model with 2 categories per node
#' set.seed(123)
#' n <- 1000
#' P <- 3
#' structure <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = P)
#' data <- matrix(rbinom(n * P, size = 1, prob = c(0.2, 0.5, 0.3)[rep(1:P, each = n)]), nrow = n, ncol = P)
#' fit <- dmrfit(data, structure = structure, with_prior = TRUE, savage_dickey = TRUE, M = 10000, ncores = 2)
#' summary(fit)
#'
#' @export 
#' 
dmrfit <- function(data, parinit = NULL, structure = NULL, with_prior = FALSE, savage_dickey = FALSE, M = 1000, oversampling = 10, ncores = 1, thresholds_alpha = 0.5, thresholds_beta = 0.5, interactions_location = 0.0, interactions_scale = 2.5, proposal_df = 5, seed = 123) {

    # save the matched call for print and summary methods
    cl <- match.call()

    # validate savage_dickey requires with_prior
    if (savage_dickey && !with_prior)
        stop("savage_dickey = TRUE requires with_prior = TRUE.")

    # processing input arguments

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

    # if n_categories is a single number, replicate it for all variables
    if(length(n_categories) == 1) {
        n_categories <- rep(n_categories, P)
    }
    else if(length(n_categories) != P) {
        stop("Length of n_categories must be either 1 or equal to the number of variables (columns) in data.")
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

    # if n_cores is less than 1, set it to 1 or the number of available cores minus one, or the number of cores specified by the user, whichever is smallest
    if(ncores < 1) {
        ncores <- min(c(ncores,parallel::detectCores()-1, 1))
    }

    # cross-product terms for the pairwise associations
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    if(is.null(structure)) {
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    } else {
        structure_input_optimize <- c(rep(1, n_thresholds), structure[lower.tri(structure, diag = FALSE)])
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize_with_structure(data = data, parinit = parinit, n_categories =  n_categories, P = P, structure = structure_input_optimize, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps) , n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale), error = function(e) {NULL}))
    }

    if(is.null(pmles)) {
        warning("Optimization failed. Returning NULL.")
        return(NULL)
    }

    # metadata needed by print/summary
    pmles$call <- cl
    pmles$P <- P
    pmles$n_categories <- n_categories
    pmles$N <- nrow(data)
    pmles$with_prior <- with_prior
    pmles$structured <- !is.null(structure)
    pmles$structure <- structure
    pmles$ncores <- ncores

    # label the parameter vector
    n_thresholds <- sum(n_categories - 1)
    thresh_names <- unlist(lapply(seq_len(P), function(p) {
        paste0("mu[", p, ",", seq_len(n_categories[p] - 1), "]")
    }))
    inter_names <- unlist(lapply(1:(P - 1), function(j) {
        lapply((j + 1):P, function(i) paste0("sigma[", i, ",", j, "]"))
    }))
    names(pmles$argument) <- c(thresh_names, inter_names)

    # --- Savage-Dickey via BSIR ---
    if(savage_dickey){

        pars <- pmles$argument
        Sigma <- pmles$utils$HW
        se <- sqrt(diag(Sigma))

        # draw M samples from proposal
        M_importance <- M * oversampling  # draw more samples than needed to ensure enough effective samples after weighting
        # draw random multivariatenormal samples (n_pars x M_importance)
        Z <- dmrfit:::mvnrnd_arma(mu = rep(0, n_pars), Sigma = Sigma, n = M_importance)
        # draw chi-squared scaling factors
        V <- rchisq(n = M_importance, df = proposal_df)

        # scale each column: Z / sqrt(V/nu), then shift by mu
        samples <- sweep(Z, 2, sqrt(proposal_df/V), "*")
        samples <- sweep(samples, 1, pars, "+")  # n_pars x M_importance

        # log-density of the proposal at each sample (multivariate log t density)
        L <- chol(Sigma)
        log_det_Sigma <- 2 * sum(log(diag(L))) # this is log determinant of Sigma calculated as 2 * sum of log-diagonal-elements of the Cholesky factor
        diff_mat <- sweep(samples, 1, pars, check.margin = FALSE)  # n_pars x M_importance
        solve_L <- backsolve(L, diff_mat, transpose = TRUE)
        mahal <- colSums(solve_L^2)
        log_q <- lgamma((proposal_df + n_pars)/2) - lgamma(proposal_df/2) - (n_pars/2) * log(proposal_df*pi) - 0.5*log_det_Sigma - ((proposal_df + n_pars) / 2) * log(1 + mahal / proposal_df)

        log_target <- numeric(M_importance)
        for(m in seq_len(M_importance)){
            log_target[m] <- -dmrfit:::npseudologlik(
                pars = samples[, m],
                data = data,
                P = P,
                n_categories = n_categories,
                with_prior = TRUE,
                ncores = ncores,
                thresholds_alpha = thresholds_alpha,
                thresholds_beta = thresholds_beta,
                interactions_location = interactions_location,
                interactions_scale = interactions_scale
            )
        }

        irlw <- (log_target - max(log_target) + min(log_q)) - log_q  # stabilize weights by subtracting max log_target and adding min log_q
        s_irlw <- log(sum(exp(irlw)) - exp(irlw)) # ISIR (Skare et al., 2003 Scandinavian Journal of Statistics)
        irlw <- irlw - s_irlw
        w <- exp(irlw) # importance resampling weights (irw)


        idx <- sample(x = 1:M_importance, size = M, replace = TRUE, prob = w)
        sir_samples <- samples[, idx, drop = FALSE]

        sir_var <- cov(t(sir_samples))
        logdetX <- determinant(sir_var, logarithm = TRUE)
        logdetZ <- determinant(Sigma, logarithm = TRUE)
        logX <- as.numeric(logdetX$modulus)
        logZ <- as.numeric(logdetZ$modulus)
        ess <- M * exp((logX - logZ) / n_pars) # * (logdetX$sign / logdetZ$sign)

        if (!is.null(structure)) {
            free_inter <- structure[lower.tri(structure, diag = FALSE)] == 1
        } else {
            free_inter <- rep(TRUE, n_interactions)
        }
        inter_idx <- (n_thresholds + 1):n_pars
        free_inter_idx <- inter_idx[free_inter]

        log_prior_at_zero <- dcauchy(0, location = interactions_location, scale = interactions_scale, log = TRUE)

        bf_01 <- numeric(length(free_inter_idx))
        names(bf_01) <- names(pars)[free_inter_idx]

        for (k in seq_along(free_inter_idx)){
            j <- free_inter_idx[k]
            sir_samples_j <- sir_samples[j, ]
            range_j <- range(sir_samples_j)
            d <- density(sir_samples_j, n = 1024, from = min(range_j), to = max(range_j))
            log_post_at_zero <- log(approx(x = d$x, y = d$y, xout = 0.0)$y)
            log_post_at_zero <- ifelse(is.na(log_post_at_zero), log(.Machine$double.eps), log_post_at_zero) #if NA, it means 0.0 is not in range and therefore it is -Inf
            log_post_at_zero <- ifelse(log_post_at_zero == Inf, -log(.Machine$double.eps), log_post_at_zero)
            bf_01[k] <- exp(log_post_at_zero - log_prior_at_zero)
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
            M = M
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
#'   \code{summary.lm}. Standard errors are obtained from the Huber-White sandwich estimator.
#' @param object a \code{dmrfit} object.
#' @param ... further arguments (currently unused).
#' @method summary dmrfit
#'
#' @return an object of class \code{summary.dmrfit} containing:
#'   \item{call}{the matched call.}
#'   \item{coefficients}{a matrix with columns for the estimate, standard error, z-value and p-value.}
#'   \item{thresholds}{coefficient matrix for threshold parameters.}
#'   \item{interactions}{coefficient matrix for free pairwise interaction parameters.}
#'   \item{neg_pseudo_loglik}{the negative pseudo-loglikelihood at convergence.}
#'   \item{P}{number of nodes.}
#'   \item{N}{number of observations.}
#'   \item{n_categories}{vector of category counts per node.}
#'   \item{with_prior}{logical, whether a prior was used.}
#'   \item{structured}{logical, whether the network structure was constrained.}
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

    # standard errors from Huber-White sandwich estimator
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
        savage_dickey = object$savage_dickey
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
    cat("Standard errors: Huber-White sandwich estimator\n")
    cat(paste0(rep("-", min(60, getOption("width"))), collapse = ""), "\n")

    cat("\nThresholds:\n")
    printCoefmat(x$thresholds, P.values = TRUE, has.Pvalue = TRUE,
                 signif.stars = FALSE, ...)

    cat("\nPairwise interactions:\n")
    printCoefmat(x$interactions, P.values = TRUE, has.Pvalue = TRUE,
                 signif.stars = TRUE, ...)

    cat("\nNegative pseudo-loglikelihood:", round(x$neg_pseudo_loglik, 4), "\n")

    if (!is.null(x$savage_dickey)) {
        sd <- x$savage_dickey
        cat("\nSavage-Dickey density ratio  [prior: Cauchy(",sd$interactions_location,",",
            sd$interactions_scale, ")]\n")
        cat("H0: sigma = 0 for each pairwise interaction\n\n")
        tbl <- cbind(Estimate = sd$estimate, SE = sd$se,
                     BF_01 = sd$bf_01, `Pr(=0|x)` = sd$pr_null)
        print(round(tbl, 4))
        cat("\nSIR samples:", sd$M,
            "  Effective sample size:", round(sd$ess, 1), "\n")
    }

    invisible(x)
}


