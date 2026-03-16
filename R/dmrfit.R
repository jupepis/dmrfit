# process_output <- function(pmles, P, n_categories, structure) {

# }


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
dmrfit <- function(data, n_categories, parinit, structure = NULL, with_prior = FALSE, ncores = 1) {

    # processing input arguments
    
    # if parinit is not provided, initialize it to a vector of zeros
    if(is.null(parinit)) {
        parinit <- rep(0.0, sum(P * (n_categories - 1)) + P * (P - 1) / 2)
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

    # number of variables
    P <- ncol(data)

    # cross-product terms for the pairwise associations
    cross_product_stats <- t(apply(data,1,function(x) {
        S <- x%*%t(x)
        S[lower.tri(S,diag=FALSE)]
    }))
    data <- cbind(data, 2.0 * cross_product_stats)

    if(is.null(structure)) {
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize(data = data, parinit = parinit, n_categories =  n_categories, P = P, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps), n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores), error = function(e) {NULL}))
    } else {
        pmles <- suppressWarnings(tryCatch(expr = dmrfit:::optimize_with_structure(data = data, parinit = parinit, n_categories =  n_categories, P = P, structure = structure, f_term = sqrt(.Machine$double.eps), m_term = sqrt(.Machine$double.eps) , n_iter_max = 100, rinit = 1.0, rmax = 10.0, with_prior = with_prior, epsilon = 1e-06, ncores = ncores), error = function(e) {NULL}))
    }

    if(is.null(pmles)) {
        warning("Optimization failed. Returning NULL.")
        return(NULL)
    }

    # metadata needed by print/summary
    pmles$P <- P
    pmles$n_categories <- n_categories
    pmles$N <- nrow(data)
    pmles$with_prior <- with_prior
    pmles$structured <- !is.null(structure)

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

#######################################################################################
#######################################################################################
