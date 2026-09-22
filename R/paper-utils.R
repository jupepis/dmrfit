#' generate_smallworld_graph (internal)
#' @description A function that generates a graph structure following the smallworld type structure, Watts-Strogatz model
#' @param p [integer], number of variables in the graph
#' @return A randomly generated smallworld structure
#' @noRd
.generate_smallworld_graph <- function(p) {
    prob    <- 0.2
    S       <- 0
    max_nei <- floor((p - 1) / 2)
    nei     <- sample(x = 1:max_nei, size = 1) # random nei (from igraph documentation: "the neighborhood within which the vertices of the lattice will be connected")
    if (p < 3) nei <- 1 # [deprecated because we only generate for p > 3] for small p's, set nei to 1 to avoid errors in igraph::sample_smallworld
    net <- NULL
    while (S < 1 && prob <= 1) {
        net <- igraph::sample_smallworld(dim = 1, size = p, nei = nei, p = prob)
        C   <- igraph::transitivity(net, type = "average")
        L   <- igraph::mean_distance(graph = net, unconnected = TRUE)
        Cr_vec <- rep(0, 100)
        Lr_vec <- rep(0, 100)

        for (i in 1:100) {
            rand_net    <- igraph::sample_gnm(igraph::vcount(net), igraph::ecount(net))
            Cr_vec[i]   <- igraph::transitivity(rand_net, type = "average")
            Lr_vec[i]   <- igraph::mean_distance(graph = rand_net, unconnected = TRUE)
        }

        Cr      <- mean(Cr_vec)
        Lr      <- mean(Lr_vec)
        S       <- if (Cr == 0 || Lr == 0) 0 else (C / Cr) / (L / Lr) 
        prob    <- prob + 0.1
    }

    if (S < 1) {
        warning(sprintf("generate_smallworld_graph: no graph reached small-worldness S >= 1 (last S = %.3f); returning closest candidate.", S), call. = FALSE)
    }

    G <- as.matrix(igraph::as_adjacency_matrix(net, sparse = FALSE))
    return(G)
}

#' generate_random_graph (internal)
#' @description A function that generates a graph structure following the Erdős-Rényi model
#' @param p [integer], number of variables in the graph
#' @return A randomly generated random structure
#' @noRd
.generate_random_graph <- function(p) {
    prob    <- stats::runif(1, 0, 1) # random inlcusion probability
    net     <- igraph::sample_gnp(n = p, p = prob)
    G       <- as.matrix(igraph::as_adjacency_matrix(net, sparse = FALSE))
    return(G)
}

#' initialize_structure (internal)
#' @description This function initializes the structure for simulation
#' @param structure structure name
#' @param p number of variables in the graph
#' @param n_categories number of categories per node
#' @return  vector o parameters to use for design condition: it includes the thresholds and only the PRESENT pairwise associations in the generated graph
#' @noRd
.initialize_structure <- function(structure, p, n_categories){
    # 1. Generate structure (0/1 of the interactions)
    if(structure == "smallworld"){
        str_coeff <- .generate_smallworld_graph(p = p) 
        str_coeff <- str_coeff[lower.tri(str_coeff, diag = FALSE)]
    }
    else if(structure == "random"){
        str_coeff <- .generate_random_graph(p = p) 
        str_coeff <- str_coeff[lower.tri(str_coeff, diag = FALSE)]
    }
    else if(structure == "full"){
        str_coeff <- rep(1.0, p * (p-1) / 2) 
    }
    # 2. Add thresholds (1.0 for each threshold)
    str_coeff <- c(rep(1.0, sum(n_categories-1)), str_coeff)
    return(str_coeff)
}


#' simulate_sample (internal)
#' @description This function simulate a single sample 
#' @param X starting value for random sample generation (matrix of size N x P)
#' @param pars vector of pars to use as true parameters to simulate networks
#' @param N sample size
#' @param P number of nodes
#' @param n_categories number of categories for each node in the network
#' @param iter number of iterations
#' @return  matrix of sampled networks
#' @noRd
.simulate_sample <- function(X, pars, N, P, n_categories, iter = 5e03){

    n_thresholds <- sum(n_categories - 1) # total number of threshold parameters
    H <- max(n_categories - 1) # maximum number of thresholds across variables

    # arrange thresholds in a [P x H] matrix
    thresholds <- matrix(0.0, nrow = P, ncol = H) 
    i <- 1
    for(p in 1:P){
        for(l in 1:(n_categories[p]-1)){
            thresholds[p,l] <- pars[i]
            i <- i+1
        }
    }

    # arrange interactions in a [P x P] matrix
    interactions <- matrix(0, nrow = P, ncol = P)
    interactions[lower.tri(interactions,diag=FALSE)] <- pars[c((n_thresholds+1):length(pars))]
    interactions[upper.tri(interactions,diag=FALSE)] <- t(interactions)[upper.tri(interactions,diag=FALSE)]

    # generate sample from Ordinal MRF using Gibbs sampler (using X as starting value)
    rnd_sample <- cpp_gibbs_sampler_omrf(mu = thresholds, sigma = interactions, 
                                        n_categories = n_categories, N = N, P = P, 
                                        iter = iter, X_start = X, save_iter = FALSE) # or generate without starting value via cpp_gibbs_sampler_omrf(mu = thresholds_mat, sigma = interactions_mat, n_categories = n_categories, N = N, P = P, iter = iter, X = NULL, save_iter = FALSE)

    return(rnd_sample$X)
}

#' process_data (internal)
#' @description This function processes the simulated data
#' @param data the input data matrix
#' @param suff_stats a logical indicating whether to compute sufficient statistics (default is FALSE)
#' @return the processed data matrix
#' @noRd
.process_data <- function(data, suff_stats = FALSE){
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
    colnames(data) <- paste0("X", 1:P) # name the columns 
    # add sufficient statistics to the data if suff_stats is TRUE
    if (suff_stats) {
        cross_product_stats <- t(apply(data,1,function(x) {  # cross-product terms for the pairwise associations
            S <- x%*%t(x)
            S[lower.tri(S,diag=FALSE)]
        }))
        data <- cbind(data, 2.0 * cross_product_stats)
    }

    return(list(
        data = data,
        n_categories = n_categories,
        n_thresholds = n_thresholds,
        n_interactions = n_interactions,
        n_pars = n_pars
    ))
}

#' validate_generated (internal)
#' @param X integer matrix, N x P, generated data (categories 0..m-1 per column)
#' @param truth one element of truth_info: must carry $n_categories (per-node cat counts)
#' @return list(collapsed, min_coverage, pair_stat_range)
#' @noRd
.validate_generated <- function(X, truth) {
    X <- as.matrix(X)
    P <- ncol(X)
    # realized distinct categories per node
    uniq <- vapply(seq_len(P), function(j) length(unique(X[, j])), integer(1L))
    # aggregate pairwise sufficient statistics: sum_i x_ip x_iq
    Sxx  <- crossprod(X)                              # P x P, off-diagonals are the pair sums
    return(list(
        collapsed       = any(uniq == 1L),            # any node frozen in ONE category
        min_coverage    = min(uniq / truth$n_categories),  # realized / possible, worst node
        pair_stat_range = range(Sxx[lower.tri(Sxx)])  # sanity range of pairwise stats
    ))
}

#' save_and_move_object (internal)
#' @description This function saves an object and moves it to a specified location
#' @param obj the object to be saved
#' @param folder the top-level project folder (e.g. "omrf/")
#' @param subfolder the phase/method-specific path segment, e.g.
#'   "simulate_data", "reference/condition_5", "PPH", "CoRe/condition_5"
#' @param results_name the file name to save under (no path components)
#' @return NULL
#' @noRd
.save_and_move_object <- function(obj, folder, subfolder, results_name){ 
    save_path <- paste0(folder, subfolder, "/")
    cp_from      <- paste0('cp "$TMPDIR"/', save_path)
    cp_to        <- paste0(' "$HOME"/', save_path)
    saveRDS(obj,file = paste0(save_path, results_name)) 
    system(command = paste(cp_from, results_name, cp_to,sep="")) 
}

#' build_conditions
#' @description This function builds the full set of simulation conditions
#' @param P integer vector of network sizes (number of variables)
#' @param N sample sizes; a vector applied to every P, or a function of P returning the sample sizes for that P
#' @param structures character vector of graph structures ("full", "random", "smallworld")
#' @return a data frame of conditions with columns index, P, structure, N, reference ("exact" or "dmh"),
#'   and dmh_bridge (TRUE only for P = 12, where the "reference" run also produces a DMH-vs-exact
#'   comparison as a byproduct -- no separate method = "DMH" run is needed or allowed for these rows)
#' @export
build_conditions <- function(P = c(6, 9, 12, 24), N = c(500, 1000, 2000, 3000), structures = c("full", "random", "smallworld")) {
  base <- expand.grid(P = P, structure = structures, N = N,
                      stringsAsFactors = FALSE)

    # This function selects the reference method for a given number of nodes in the network (network size)
    .reference_for_p <- function(p) {
        if (p <= 12) "exact"   # exact is generated for P <= 12; for P = 12 this run also yields DMH as a byproduct (see dmh_bridge)
        else         "dmh"     # p = 24 (and any p > 12)
    }

  conditions <- cbind(base, reference = vapply(base$P, .reference_for_p, character(1L)),
                       row.names = NULL)
  conditions$dmh_bridge <- conditions$P == 12

  # number the runnable rows -- one per P/N/structure combination, no fan-out
  conditions <- cbind(index = seq_len(nrow(conditions)), conditions)
  return(conditions)
}


#' select_bivariate_pairs (internal)
#' @description Randomly selects one representative parameter pair from each
#'   of three categories (interaction-interaction, threshold-interaction,
#'   threshold-threshold) for bivariate overlap computation. For
#'   interaction-interaction and threshold-interaction, restricts selection to
#'   pairs sharing a variable (the pairs most likely to show real posterior
#'   correlation, and most directly testing CoRe's covariance correction) and
#'   to interactions actually present in the true structure.
#' @param structure the 0/1 mask of present interactions (from .initialize_structure,
#'   length P*(P-1)/2, lower-triangle order)
#' @param P number of variables
#' @param n_thresholds number of threshold parameters
#' @return a named list: $interaction_interaction, $threshold_interaction,
#'   $threshold_threshold -- each a length-2 integer vector of parameter
#'   indices (1-based, into the full n_pars-length parameter vector), or NULL
#'   if no valid pair exists for that category
.select_bivariate_pairs <- function(structure, P, n_thresholds) {

    # Generate the lower-triangle indices for the interactions
     lower_indices <- { # same as the function sigma_lower_tri_indices() in paper-utils.cpp
        idx <- matrix(0L, nrow = 2, ncol = P * (P - 1) / 2)
        l <- 1
        for (j in 1:(P - 1)) {
            for (i in (j + 1):P) {
                idx[1, l] <- i
                idx[2, l] <- j
                l <- l + 1
            }
        }
        idx
    }
 

    present <- which(structure == 1)   # indices (into the interaction block) of present edges

    out <- list(interaction_interaction = NULL,
                threshold_interaction   = NULL,
                threshold_threshold     = NULL)

    if (length(present) >= 2) {
        # --- interaction-interaction: two present edges sharing a node ---
        found <- FALSE
        candidates <- present
        for (a in seq_along(candidates)) {
            for (b in seq_along(candidates)) {
                if (a >= b) next
                ij_a <- lower_indices[, candidates[a]]
                ij_b <- lower_indices[, candidates[b]]
                if (length(intersect(ij_a, ij_b)) > 0) {   # share a node
                    out$interaction_interaction <- n_thresholds + c(candidates[a], candidates[b])
                    found <- TRUE
                    break
                }
            }
            if (found) break
        }
        # The fallback: no shared-node pair exists (sparse/disconnected structure) -> pick any two present edges
        if (!found && length(present) >= 2) {
            picked <- sample(present, 2)
            out$interaction_interaction <- n_thresholds + picked
        }
    }

    if (length(present) >= 1 && n_thresholds >= 1) {
        # --- threshold-interaction: a present edge's threshold + that edge ---
        picked_edge <- present[sample(length(present), 1)]
        ij <- lower_indices[, picked_edge]
        node <- ij[sample(2, 1)]                 # pick one of the edge's two endpoints
        thr_idx <- node                            # HOOK: adapt to your actual threshold indexing
        # (if multiple thresholds per variable -- ordinal case -- pick one, e.g. the first)
        out$threshold_interaction <- c(thr_idx, n_thresholds + picked_edge)
    }

    if (n_thresholds >= 2) {
        # --- threshold-threshold: any two distinct thresholds ---
        picked <- sample(1:n_thresholds, 2)
        out$threshold_threshold <- picked
    }

    out
}

#' @title marginal_overlapping_index (overlapping index, only on marginal distributions)
#'
#' @description computes the distribution-free overlapping index (cite here [..])
#'
#' @param d1 vector of posterior draws for one model parameter
#' @param d2 vector of posterior draws for one model parameter for the gold standard method
#' @param range_d2 range of values for d2
#' 
#' @return  a value quantifying the divergence, between 0 and 1
#' 
#' @export
#' 
marginal_overlapping_index <- function(d1, d2, range_d2){

    # calculate upper and lower bounds for the joint support over the two distributions, d1 and d2
    bounds <- range(c(d1, range_d2))
    from <- bounds[1]
    to <- bounds[2]

    # calculate density of d2
    density_d2 <- density(x = d2, n = 1024, from = from, to = to)
    density_d2$y[which(density_d2$y < .Machine$double.xmin)] <- .Machine$double.xmin
    density_d2$y <- density_d2$y / sum(density_d2$y)

    # calculate density of d1
    density_d1 <- density(x = d1, n = 1024, from = from, to = to)  
    density_d1$y[which(density_d1$y < .Machine$double.xmin)] <- .Machine$double.xmin
    density_d1$y <- density_d1$y / sum(density_d1$y)
     
    integral_min <- sum(pmin(density_d1$y, density_d2$y))
    integral_max <- sum(pmax(density_d1$y, density_d2$y))
  
    return(integral_min / integral_max)
}

#' @title bivariate_overlapping_index (overlapping index on the joint 2D distribution)
#'
#' @description computes the distribution-free overlapping index on the JOINT
#'   density of a pair of parameters.
#'
#' @param d1 matrix/data.frame, n x 2, posterior draws for the parameter pair
#' @param d2 matrix/data.frame, n x 2, posterior draws for the parameter pair
#'   (gold standard / reference method)
#' @param n_grid number of grid points per dimension for the 2D KDE (default 128;
#'   kde2d is O(n_grid^2), so keep this modest -- 1024 (matching the marginal
#'   function) would be 1024^2 grid cells, likely too slow to run per pair per replicate)
#'
#' @return a value quantifying the joint overlap, between 0 and 1
#'
#' @export
bivariate_overlapping_index <- function(d1, d2, n_grid = 128) {

    # joint bounding box over both samples, so both KDEs are evaluated on the same grid
    x_range <- range(c(d1[, 1], d2[, 1]))
    y_range <- range(c(d1[, 2], d2[, 2]))

    kde1 <- MASS::kde2d(d1[, 1], d1[, 2], n = n_grid, lims = c(x_range, y_range))
    kde2 <- MASS::kde2d(d2[, 1], d2[, 2], n = n_grid, lims = c(x_range, y_range))

    z1 <- kde1$z
    z1[z1 < .Machine$double.xmin] <- .Machine$double.xmin
    z1 <- z1 / sum(z1)

    z2 <- kde2$z
    z2[z2 < .Machine$double.xmin] <- .Machine$double.xmin
    z2 <- z2 / sum(z2)

    # same min/max-overlap logic as the marginal version, applied cell-by-cell
    integral_min <- sum(pmin(z1, z2))
    integral_max <- sum(pmax(z1, z2))

    return(integral_min / integral_max)
}

.log_dcauchy <- function(x, scale, location) {
    # log-density of the Cauchy distribution at x
    log_density <- log(1 / (pi * scale * (1 + ((x - location) / scale)^2)))
    return(log_density)
}

# .compute_metrics_utils (internal) only works one replicate (z) because it will be run within the foreach() for each method
.compute_metrics_utils <- function(x, reference, method, n_pars, n_thresholds, structure, P, interactions_location, interactions_scale) {

    # compute_metrics() is called within each replicate (inside foreach):
    # 1. posterior quantities (mean, median, sd, map)
    # 2. acceptance rate (if applicable)
    # 3. time (seconds elapsed)
    # 4. counter update (if AdaCoRe)
    # 5. ESS
    # 6. Savage-Dickey density ratio
    # 7. metric: marginal and bivariate overlapping index
    # 8. Return as a list of results for replicate z, to be combined later acrtoss replicates (z = 1:Q)

    # --- Posterior quantities (mean, median, sd, map) ---
    post <- matrix(0.0, nrow = n_pars, ncol = 4)
    post[, 1] <- apply(x$draws, 2, mean)
    post[, 2] <- apply(x$draws, 2, median)
    post[, 3] <- apply(x$draws, 2, sd)
    post[, 4] <- apply(x$draws, 2, function(z) {
        d_x <- density(z)
        return(d_x$x[which.max(d_x$y)])
    })
    colnames(post) <- c("mean", "median", "sd", "map")

    # --- Acceptance rate (applicable to all methods except for Post-Hoc's) ---
    acceptance <- if(!(method %in% c("RM", "GHW", "MCH"))) {
        x$acceptance
    } else {
        NA_real_
    }
    
    # --- Time (seconds elapsed) ---
    time <- x$seconds_elapsed

    # --- Counter update (only for AdaCoRe) ---
    counter_update <- if(method == "AdaCoRe") {
        x$counter_update
    } else {
        NA_integer_
    }

    # --- ESS (Effective Sample Size, using coda::effectiveSize) ---
    ESS <- rep(0.0,n_pars)
    for(j in 1:n_pars){
        y_j <- x$draws[,j]
        ESS[j] <- coda::effectiveSize(x = y_j)
    }

    # --- Savage-Dickey density ratio (only for interaction parameters) ---
    log_sd_density_ratio <- rep(0.0, n_pars - n_thresholds)
    log_prior_at_zero <- .log_dcauchy(x = 0.0, scale = interactions_scale, location = interactions_location)
    for(j in (n_thresholds + 1):n_pars){
        # Extract the draws for the j-th parameter
        y_j <- x$draws[,j]

        # Define the range for density estimation
        range_y_j <- range(y_j)

        # Find density estimate for the j-th parameter
        d_j <- density(x = y_j, n = 1024, from = min(range_y_j), to = max(range_y_j))

        # Define approximation function for the density
        d_j_approxfun <- approxfun(x = d_j$x, y = d_j$y)

        # Calculate the log posterior probability at 0.0 using the approximation function, and correction for too large or too low numbers
        log_post_at_zero <- log(d_j_approxfun(0.0))
        log_post_at_zero <- ifelse(is.na(log_post_at_zero) | log_post_at_zero == -Inf, log(.Machine$double.eps), log_post_at_zero)

        # [[OLD CODE: refers to the two lines below]]
        # post_pr <- ifelse(is.na(post_pr), log(.Machine$double.eps), post_pr) # if NA or Inf, it means 0.0 is not in range and therefore it is -Inf
        # post_pr <- ifelse(post_pr == Inf, - log(.Machine$double.eps), post_pr)

        # Compute logarithm of the Savage-Dickey Density ratio
        log_sd_density_ratio[j - n_thresholds] <- log_post_at_zero - log_prior_at_zero
    }

    # --- Metric: marginal and bivariate overlapping index (only for methods that are not exact) ---
    metric_marginal <- NULL
    metric_bivariate <- NULL
    if(method != "exact"){ # we need to load the exact draws from the rds files

        # --- Marginal metric: for each parameter, compute the marginal overlapping index between the posterior draws and the exact draws
        metric_marginal <- rep(0.0,n_pars)
        for(j in 1:n_pars){
            y_j <- reference$draws[,j]
            range_y_j <- range(y_j)
            x_m_j <- x$draws[,j]
            range_x_m_j <- range(x_m_j)
            metric_marginal[j] <- marginal_overlapping_index(d1 = x_m_j, d2 = y_j, range_d2 = range_y_j)
        }  

        # --- Bivariate metric: for a random but motivated selection of pairs of parameters (not all of them because expensive; see '.select_bivariate_pairs' for more details on the selection): computes the bivariate overlapping index between the posterior draws and the exact draws
        pairs <- .select_bivariate_pairs(structure, P, n_thresholds)  # computed once per condition/truth, passed in
        metric_bivariate <- lapply(pairs, function(pr) {
            if (is.null(pr)) return(NA_real_)
            bivariate_overlapping_index(d1 = x$draws[, pr], d2 = reference$draws[, pr])
        })
        names(metric_bivariate) <- names(pairs)  # "interaction_interaction", "threshold_interaction", "threshold_threshold"
    }

    # --- Return the results as a list ---
    return(list(
        post = post,
        acceptance = acceptance,
        time = time,
        counter_update = counter_update,
        ESS = ESS,
        log_sd_density_ratio = log_sd_density_ratio,
        metric_marginal = metric_marginal,
        metric_bivariate = metric_bivariate
    ))
}

# .compute_metrics
.compute_metrics <- function(x, reference, method, n_pars, n_thresholds, structure, P, interactions_location, interactions_scale) {
    # compute_metrics() is called within each replicate (inside foreach):
    # 1. posterior quantities (mean, median, sd, map)
    # 2. acceptance rate (if applicable)
    # 3. time (seconds elapsed)
    # 4. counter update (if AdaCoRe)
    # 5. ESS
    # 6. Savage-Dickey density ratio
    # 7. metric: marginal and bivariate overlapping index
    # Return as a list of results for replicate z, to be combined later across replicates (z = 1:Q)

    metrics <- NULL

    if(all(names(reference) %in% c("exact", "dmh")) && P == 12) { # if P = 12
        # reference is a list of two elements (exact and dmh)
        metrics <- lapply(reference, function(ref) {
            .compute_metrics_utils(x = x, reference = ref,
                                   method = method, n_pars = n_pars,
                                   n_thresholds = n_thresholds,
                                   structure = structure,
                                   P = P,
                                   interactions_location = interactions_location,
                                   interactions_scale = interactions_scale)
        })
    } else { # if P != 12
        metrics <- .compute_metrics_utils(x = x, reference = reference,
                                            method = method, n_pars = n_pars,
                                            n_thresholds = n_thresholds,
                                            structure = structure,
                                            P = P,
                                            interactions_location = interactions_location,
                                            interactions_scale = interactions_scale)
    }

    return(metrics)
}

# .process_metrics (internal) --  to be written as post processing function when results are in

