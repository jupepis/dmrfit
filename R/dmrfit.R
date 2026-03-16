#' @title dmrfit
#' 
#' @description Fit a discrete Markov Random Field model via pseudo-likelihood estimation. The optimization is performed using the trust region algorithm implemented by 
#' 
#' @param data data matrix, with rows as samples and columns as variables. Each variable should be rescaled to the range of 0 to m-1, where m is the number of categories for that variable. The baseline category is always the minimum value in the variable. The internal processing will check if the variables are rescaled and will rescale them if necessary. If there are any NAs in the data, they will be removed before optimization (listwise deletion).
#' @param n_categories n_categories, a vector of length equal to the number of variables, where each element specifies the number of categories for the corresponding variable. If a single number is provided, it will be replicated for all variables.
#' @param parinit parinit, initial parameter estimates, a vector of length equal to the number of parameters in the model. If NULL, it will be initialized to a vector of zeros. The number of parameters in the model is calculated as sum(P * (n_categories - 1)) + P * (P - 1) / 2, where P is the number of variables.
#' @param structure network structure, P x P matrix, with 0 for no edge and 1 for edge. If NULL, then fully connected graph is assumed. Default is NULL
#' @param with_prior with_prior, logical, whether to include the prior in the optimization. If TRUE, a Beta-Prime and a Cauchy prior are applied to thresholds and pairwise associations respectively. If FALSE, no prior is applied. Default is FALSE
#' @param ncores ncores, number of cores to use for parallel computing of the gradient, hessian and likelihood values. If ncores is larger than the number of available cores, it will be set to the number of available cores minus one. Default is 1.
#' @return dmrfit S3 class object
#' @export 
#' 
dmrfit <- function(data, n_categories, parinit = NULL, structure = NULL, with_prior = FALSE, ncores = 1) {

    cl <- match.call()

    # processing input arguments

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

    # check if data is a matrix, if not convert it to a matrix
    if(!is.matrix(data)) {
        data <- as.matrix(data)
    }

    # remove NAs from data if any exist
    data <- data[!is.na(rowSums(data)), ]

    # check if variables are rescaled from 0 to m-1
    for(i in 1:ncol(data)) {
        if(any(data[, i] < 0) || any(data[, i] >= n_categories[i])) {
            stop(paste("Variable", i, "contains values outside the range of 0 to n_categories - 1. Please rescale the variable accordingly."))
        }
    }

    # rescale the data to the range of 0 to m-1 if necessary
    for(i in 1:ncol(data)) {
        if(any(data[, i] < 0) || any(data[, i] >= n_categories[i])) {
            data[, i] <- data[, i] - min(data[, i]) # baseline category is always the minimum value in the variable
        }
    }

    # cross-product terms for the pairwise associations
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    if(is.null(structure)) {
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores), error = function(e) {NULL}))
    } else {
        structure_input_optimize <- c(rep(1, n_thresholds), structure[lower.tri(structure, diag = FALSE)])
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize_with_structure(data = data, parinit = parinit, n_categories =  n_categories, P = P, structure = structure_input_optimize, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps) , n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores), error = function(e) {NULL}))
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

    # label the parameter vector
    n_thresholds <- sum(n_categories - 1)
    thresh_names <- unlist(lapply(seq_len(P), function(p) {
        paste0("mu[", p, ",", seq_len(n_categories[p] - 1), "]")
    }))
    inter_names <- unlist(lapply(1:(P - 1), function(j) {
        lapply((j + 1):P, function(i) paste0("sigma[", i, ",", j, "]"))
    }))
    names(pmles$argument) <- c(thresh_names, inter_names)

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

    P  <- object$P
    n_categories <- object$n_categories
    n_thresholds   <- sum(n_categories - 1)
    n_interactions <- P * (P - 1) / 2
    n_pars <- n_thresholds + n_interactions
    pars   <- object$argument

    # standard errors from Huber-White sandwich estimator
    se <- sqrt(diag(object$utils$HW))
    names(se) <- names(pars)

    # build full coefficient table
    z_val  <- pars / se
    p_val  <- 2 * pnorm(-abs(z_val))
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

    # Build interactions table; add Savage-Dickey column when prior was used
    inter_table <- coef_table[free_inter_idx, , drop = FALSE]
    if (object$with_prior) {
        # Savage-Dickey density ratio: prior Cauchy(0, 2.5), posterior ~ Normal(est, se)
        inter_est <- pars[free_inter_idx]
        inter_se  <- se[free_inter_idx]
        prior_at_zero <- dcauchy(0, location = 0, scale = 2.5)
        post_at_zero  <- dnorm(0, mean = inter_est, sd = inter_se)
        bf_01 <- post_at_zero / prior_at_zero
        pr_null <- bf_01 / (1 + bf_01)
        inter_table <- cbind(inter_table, `Pr(=0|x)` = pr_null)
    }

    out <- list(
        call            = object$call,
        coefficients    = coef_table[c(1:n_thresholds, free_inter_idx), , drop = FALSE],
        thresholds      = coef_table[1:n_thresholds, , drop = FALSE],
        interactions    = inter_table,
        neg_pseudo_loglik = object$utils$value,
        P               = P,
        N               = object$N,
        n_categories    = n_categories,
        with_prior      = object$with_prior,
        structured      = object$structured
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
    inter <- x$interactions
    if (x$with_prior) {
        # print first 4 cols via printCoefmat, then the Pr(=0|x) column
        printCoefmat(inter[, 1:4, drop = FALSE], P.values = TRUE, has.Pvalue = TRUE,
                     signif.stars = TRUE, ...)
        cat("\nSavage-Dickey Pr(=0|x)  [prior: Cauchy(0, 2.5)]:\n")
        print(round(inter[, "Pr(=0|x)", drop = FALSE], 4))
    } else {
        printCoefmat(inter, P.values = TRUE, has.Pvalue = TRUE,
                     signif.stars = TRUE, ...)
    }

    cat("\nNegative pseudo-loglikelihood:", round(x$neg_pseudo_loglik, 4), "\n")
    invisible(x)
}

