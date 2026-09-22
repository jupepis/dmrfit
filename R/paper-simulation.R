.HYPER <- list(
    # --- Sampler hyperparameters ---
    nsim = 2e04, # number of posterior draws
    burnin = 4500, # number of burnin draws
    adaptive_stage_n_iter = 500, # number of iterations for the adaptive stage of the exact sampler (for P <= 12)
    sigma2 = 1.0, # default value for sigma2 (meant to work well for any case; however, for ordinal casa with >2 levels, we set it to 0.001)

    # --- Priors hyperparameters ---
    thresholds_alpha = 0.5, # shape parameter for the gamma prior on thresholds
    thresholds_beta = 0.5, # rate parameter for the gamma prior on thresholds
    interactions_location = 0.0, # location parameter for the normal prior on interactions
    interactions_scale = 2.5, # scale parameter for the normal prior on interactions

    # --- Gibbs sampler hyperparameters ---
    L = 1e05, # number of samples used by RM and MCH routines
    L_DMH = 25000, # number of samples to approximate gradient (used inside the DMH) 
    new_data_sample_inner_gibbs_iter = 5000, # inner Gibbs iterations for generating new data samples (using random start value)
    inner_gibbs_iter = 5, # inner Gibbs iterations for approximating gradient (using warm start value)

    # --- Robbins-Monro hyperparameters ---
    rm_step_thresholds = 0.001, # default value is set to 0.001 according to Bouranis et al.
    rm_step_interactions = 0.001, # default value is set to 0.001 according to Bouranis et al.
    rm_max_iter = 200, # max number of iterations
    tolerance = 0.0001, # tolerance for convergence

    # --- Output control ---
    verbose = FALSE, # whether to print verbose output during the simulation
    progress = FALSE # whether to print progress output during the simulation
)

#' @title simulation
#' @description compute posterior distribution according to a specific 'method' (repeating <simulate>|<reference>|<pseudo> as this saves storage space and time for saving rds of exact and posterior draws and also for memory serialization
#' @param data raw dataset with only complete cases and X1,...,XP coded in 0,...,M-1 levels (IMPORTANT!)
#' @param condition row from the build_conditions() output, which contains all the simulation conditions. It contains the following columns: N, P, structure, index, and reference method (exact or dmh)
#' @param master_seed master seed to use with RNG for reproducibility of parallelized simulations (IF WE GENERATE EVERYTIME FOR ONE METHOD TO RUN, WE NEED TO USE THE SAME MASTER SEED FOR ALL METHODS TO ENSURE THAT THE SAME RANDOM STRUCTURES ARE GENERATED)
#' @param Q number of random structures to generate (default is 10)
#' @param Q_nsim number of simulation per random structure to generate (default is 10)
#' @param nthreads number of threads to use for parallelization, default is 100
#' @param method one of the following methods: "reference" (exact likelihood for P <= 12 or DMH for P = 24), "PPH" (pseudo posterior including all Post Hoc methods, RM, GHW, and MCH; one run, four methods), "DMH" (double Metropolis-Hastings), "CoRe" (Godambe-Huber-White correction), "CoRe-RM" (Robbins-Monro correction), "CoRe-MCH" (Monte Carlo Hessian correction), "AdaCoRe" (adaptive covariance structure correction)
#' @param folder the folder where the results will be saved
#' @return  depending on "method", a summary list of results or list of draws, or helpers is saved
#' @export
simulation <- function(data, condition, master_seed, folder = "pl_project/", Q = 10, Q_nsim = 10, nthreads = 100L, method = c("reference","PPH","DMH","CoRe","CoRe-RM","CoRe-MCH","AdaCoRe")){

    # ---- Process method ----
    method  <- match.arg(arg = method, choices = c("reference","PPH","DMH","CoRe","CoRe-RM","CoRe-MCH","AdaCoRe"), several.ok = FALSE)

    # ---- Generate master seed if not provided ----
    if (is.null(master_seed)) {
        master_seed <- as.integer(Sys.time()) + as.integer(Sys.getpid())
        
    }

    # ---- Initialize a parallel cluster ----
    #cl <- makeCluster(nthreads, type = "FORK") # FORK because we want to share large objects in memory (e.g., the partition function matrix X for P <= 12)
    cl <- makeCluster(nthreads, type = "FORK", outfile = file.path(Sys.getenv("HOME"), "worker.log"))
    registerDoParallel(cl)

    # --- Set the random seed for reproducibility across parallel threads ---
    set.seed(master_seed)
    registerDoRNG(seed = master_seed) # ensure each thread has a unique seed

    # ---- Unpack HYPER parameters ----

    # ---- Sampler parameters ----
    nsim <- .HYPER$nsim
    burnin <- .HYPER$burnin
    adaptive_stage_n_iter <- .HYPER$adaptive_stage_n_iter
    sigma2 <- .HYPER$sigma2

    # ---- Priors parameters ----
    thresholds_alpha <- .HYPER$thresholds_alpha
    thresholds_beta <- .HYPER$thresholds_beta
    interactions_location <- .HYPER$interactions_location
    interactions_scale <- .HYPER$interactions_scale

    # ---- Gibbs sampler parameters ----
    L <- .HYPER$L
    L_DMH <- .HYPER$L_DMH
    new_data_sample_inner_gibbs_iter <- .HYPER$new_data_sample_inner_gibbs_iter
    inner_gibbs_iter <- .HYPER$inner_gibbs_iter

    # --- Robbins-Monro parameters ----
    rm_step_thresholds <- .HYPER$rm_step_thresholds
    rm_step_interactions <- .HYPER$rm_step_interactions
    rm_max_iter <- .HYPER$rm_max_iter
    tolerance <- .HYPER$tolerance

    # ---- Other parameters ----
    verbose <- .HYPER$verbose
    progress <- .HYPER$progress

    # ---- Initialize variables and process data ----
    N                       <- as.numeric(condition[["N"]])                  # sample size
    P                       <- as.numeric(condition[["P"]])                  # number of nodes (variables)
    S                       <- as.character(condition[["structure"]])        # structure type
    index_sim               <- as.numeric(condition[["index"]])              # index of the simulation condition
    reference               <- as.character(condition[["reference"]])        # reference condition (dmh or exact; both only for P = 12)
    design_condition        <- list(index = index_sim, N = N, P = P, S = S, reference = reference) # list to store the design condition
    proc_data               <- .process_data(data = data, suff_stats = FALSE)   # process data to include sufficient statistics
    data                    <- proc_data$data
    n_categories            <- proc_data$n_categories
    n_thresholds            <- proc_data$n_thresholds
    .equal_categories_check <- all(n_categories == n_categories[1])             # check if all variables have the same number of categories

    # ---- Print simulation start message ----
    cat(paste(Sys.time()," - simulation started. \n",sep=""))

    # ---- Phase: <simulate + reference draws> ---- note: saving .rds objects of reference draws and simulated datasets only if method is "reference"
    if(method == "reference"){
  
        # ---- Step 1: Q truths from the FULL dataset ----
        truth_info <- vector("list", Q)
        for (i in seq_len(Q)) {
            fit <- NULL
            while (is.null(fit)) {
                select_items    <- sort(sample(seq_len(ncol(data)), size = P, replace = FALSE)) # randomly select P variables (without replacement)
                ncat_i          <- n_categories[select_items]              
                n_thr_i         <- sum(ncat_i - 1L)
                n_int_i         <- P * (P - 1) / 2
                n_pars_i        <- n_thr_i + n_int_i
                str_i           <- .initialize_structure(structure = S, p = P, n_categories = ncat_i)
                X_i <- NULL
                if(!.equal_categories_check && P <= 12){ # if equal categories case is not satisfied and only for P<=12 where exact likelihood is feasible [[CHECK IF on .equal_categories_check, is it right?]]
                    permutations_i <- expand.grid(sapply(1:P, function(x) list(0:(ncat_i[x]-1)))) 
                    permutations_i <- as.matrix(permutations_i) 
                    X_i <- dmrfit:::cpp_build_permutations_stats(permutations = permutations_i, n_pars = n_pars_i, n_thresholds = n_thr_i, n_categories = ncat_i)
                }

                data_est <- data[, select_items, drop = FALSE]     # Select all rows, samples columns (P variables)
                Sxx <- t(apply(data_est, 1, function(x){ M <- x %*% t(x); M[lower.tri(M)] }))
                data_est <- cbind(data_est, 2.0 * Sxx) # multiply by 2 (sufficient statistics for pseudo-likelihood function)

                fit <- suppressWarnings(tryCatch(cpp_optimize_with_structure(data = data_est, 
                                                                            parinit = rep(0.0, n_pars_i), 
                                                                            structure = str_i, 
                                                                            n_categories = ncat_i, 
                                                                            P = P, 
                                                                            f_term = sqrt(.Machine$double.eps), 
                                                                            m_term = sqrt(.Machine$double.eps), 
                                                                            n_iter_max = 100, 
                                                                            rinit = 1.0, 
                                                                            rmax = 10.0, 
                                                                            with_prior = TRUE, 
                                                                            epsilon = 1e-06, 
                                                                            ncores = 1L, 
                                                                            thresholds_alpha = thresholds_alpha, 
                                                                            thresholds_beta = thresholds_beta, 
                                                                            interactions_location = interactions_location, 
                                                                            interactions_scale = interactions_scale),
                                                                            error = function(e) NULL)) # [[ CHECK / CHANGE FUNCTION ]]
            }
            true_pars <- fit$argument
            truth_info[[i]] <- list(index_sim = index_sim, true_pars = true_pars, select_items = select_items,
                                n_categories = ncat_i, n_thresholds = n_thr_i, n_pars = n_pars_i,
                                structure = str_i,
                                density = sum(abs(true_pars[-(seq_len(n_thr_i))]) > 0) / n_int_i, X = X_i)
        }

        # ---- Step 2: Q_nsim datasets per each generated structure, each with N observations ----
        ref_str    <- rep(seq_len(Q), each = Q_nsim) # create a reference vector to identify which random structure each dataset corresponds to
        n_datasets <- Q * Q_nsim # total number of datasets to simulate
        which_dataset_sigma2 <- sample.int(n = n_datasets, size = 5)
        samples_ls <- foreach(q = seq_len(n_datasets), .packages = "dmrfit") %dopar% {
            tr <- truth_info[[ref_str[q]]]
            max_regen <- 20L
            max_solve <- 20L
            regen  <- 0L
            pmles  <- list(argument = NULL)
            X_suff <- NULL

            repeat {
                X   <- .simulate_sample(X = NULL, pars = tr$true_pars, N = N, P = P,
                                        n_categories = tr$n_categories,
                                        iter = new_data_sample_inner_gibbs_iter)
                chk <- .validate_generated(X, tr)

                if (chk$collapsed && regen < max_regen) {   # data is degenerate --> regenerate
                    regen <- regen + 1L
                    next
                }

                X_suff <- .process_data(data = X, suff_stats = TRUE)$data
                for (s in seq_len(max_solve)) {             # data is valid --> retry solver
                    init  <- if (s == 1) rep(0.0, tr$n_pars) else stats::rnorm(tr$n_pars, 0, 0.1 * s)
                    pmles <- suppressWarnings(tryCatch(cpp_optimize(
                                data                  = X_suff,
                                parinit               = init,
                                n_categories          = tr$n_categories,
                                P                     = P,
                                f_term                = sqrt(.Machine$double.eps),
                                m_term                = sqrt(.Machine$double.eps),
                                n_iter_max            = 100,
                                rinit                 = 1.0,
                                rmax                  = 10.0,
                                with_prior            = TRUE,
                                epsilon               = 1e-06,
                                ncores                = 1L,
                                thresholds_alpha      = thresholds_alpha,
                                thresholds_beta       = thresholds_beta,
                                interactions_location = interactions_location,
                                interactions_scale    = interactions_scale),
                            error = function(e) list(argument = NULL))) # [[ CHECK / CHANGE FUNCTION ]]
                    if (!is.null(pmles$argument)) break
                }
                if (!is.null(pmles$argument) || regen >= max_regen) break
                regen <- regen + 1L
            }

            solved <- !is.null(pmles$argument)
            capped <- regen >= max_regen
            forced <- capped && !solved

            list(sample     = X_suff,
                pmles      = pmles$argument,      # may be NULL when forced == TRUE
                n_regen    = regen,
                collapsed  = chk$collapsed,
                coverage   = chk$min_coverage,
                solved     = solved,              # TRUE = normal converged pmle
                forced     = forced,              # TRUE = hit cap without a usable pmle
                ref_str     = ref_str[q],
                save_sigma2 = q %in% which_dataset_sigma2,
                cap_reason = if (!capped) NA_character_
                            else if (chk$collapsed) "collapsed" else "solver")
        }
        cat(paste(Sys.time()," - data simulation done. \n",sep=""))

        # ---- Save .rds objects of simulated datasets and truth_info ----

        # ... sample_ls 
        .save_and_move_object(obj = samples_ls, folder = folder, subfolder = "simulate_data", results_name = paste0("simulate_data_condition_", index_sim, ".rds"))

        # ... info data condition (true parameters, structure and other info)
        light_truth_info <- lapply(truth_info, function(x) {x$X <- NULL; return(x)}) # remove the X matrix from truth_info to save space
        light_truth_info$design <- list(design_condition = design_condition, n_categories = n_categories, n_thresholds = n_thresholds, Q = Q, Q_nsim = Q_nsim, nthreads = nthreads, index_sim = index_sim, master_seed = master_seed) # attach info on design to the light_truth_info object
        .save_and_move_object(obj = light_truth_info, folder = folder, subfolder = "simulate_data", results_name = paste0("info_data_condition_", index_sim, ".rds"))
 
        # ---- Step 3: <reference draws> ---- note: exact posterior (for P <= 12) or DMH (for P = 24) 
        if(P <= 12){ # use exact sampler for P <= 12           
            X <- NULL
            if(.equal_categories_check){ # if equal categories case is satisfied and only for P<=12 where exact likelihood is feasible
                ncat_P <- rep(n_categories[1], P)
                permutations <- expand.grid(sapply(1:P, function(x) list(0:(ncat_P[1] - 1)))) 
                permutations <- as.matrix(permutations) 
                n_thresholds_equal <- sum(ncat_P - 1)
                n_pars_equal <- n_thresholds_equal + P * (P - 1) / 2
                X <- dmrfit:::cpp_build_permutations_stats(permutations = permutations, n_pars = n_pars_equal, n_thresholds = n_thresholds_equal, n_categories = ncat_P) 
                rm(permutations, n_pars_equal, n_thresholds_equal, ncat_P)
            }
            draws <- foreach(z = 1:n_datasets, .packages = "dmrfit") %dopar% {
                # Unpack objects from replicate z
                sample_z <- samples_ls[[z]] # select replicate z
                data_z <- sample_z$sample[,1:P]
                pars_z <- sample_z$pmles
                ref_str_z <- sample_z$ref_str # reference structure index
                X_partition_z <- if(.equal_categories_check) X else truth_info[[ref_str_z]]$X
                n_categories_z <- truth_info[[ref_str_z]]$n_categories
                n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
                n_pars_z <- truth_info[[ref_str_z]]$n_pars
                structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
                save_sigma2_z <- sample_z$save_sigma2

                reference_draws <- list()

                # Reference draws using exact sampler
                reference_draws$exact <- dmrfit:::cpp_exact_sampler(data = data_z, pars = pars_z, n_categories = n_categories_z,
                                                                    X = X_partition_z, nsim = nsim, burnin = burnin, adaptive_stage_n_iter = adaptive_stage_n_iter,
                                                                    sigma2 = sigma2, thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta,
                                                                    interactions_location = interactions_location, interactions_scale = interactions_scale,
                                                                    verbose = verbose, progress = progress)
                if(save_sigma2_z){
                    .save_and_move_object(obj = reference_draws$exact$sigma2, folder = folder, subfolder = "Exact", results_name = paste0("exact_sigma2_condition_", index_sim, "_sample_", z, ".rds"))
                }
                # Reference draws using DMH sampler --> only for P = 12 we also add DMH and return a list with two sets of samples
                if(P == 12) {
                    reference_draws$dmh <- dmrfit:::cpp_dmh_sampler(data = data_z, pars = pars_z, n_categories = n_categories_z, 
                                                                    nsim = nsim, burnin = burnin, L = L_DMH, inner_sampler_n_iter = inner_gibbs_iter,
                                                                    adaptive_stage_n_iter = adaptive_stage_n_iter, sigma2 = sigma2, 
                                                                    thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, 
                                                                    interactions_location = interactions_location, interactions_scale = interactions_scale,
                                                                    verbose = verbose, progress = progress)
                    dmh_metrics <- .compute_metrics(x = reference_draws$dmh, 
                                                    reference = reference_draws$exact, 
                                                    method = "DMH", 
                                                    n_pars = n_pars_z, 
                                                    n_thresholds = n_thresholds_z,
                                                    structure = structure_z, 
                                                    P = P, 
                                                    interactions_location = interactions_location, 
                                                    interactions_scale = interactions_scale)
                    # -- Save Metrics (reference to exact) for DMH draws (only for P = 12) --
                    .save_and_move_object(obj = dmh_metrics, folder = folder, subfolder = "DMH", results_name = paste0("summary_condition_", index_sim, ".rds"))
                    if(save_sigma2_z){
                        .save_and_move_object(obj = reference_draws$dmh$sigma2, folder = folder, subfolder = "DMH", results_name = paste0("dmh_sigma2_condition_", index_sim, "_sample_", z, ".rds"))
                    }
                }
                # ---- Reference draws ----
                if(P == 12) {
                    .save_and_move_object(obj = reference_draws, folder = folder, subfolder = paste0(method, "/condition_", index_sim), results_name = paste0("reference_draws_sample_", z, ".rds"))
                } else {
                .save_and_move_object(obj = reference_draws$exact, folder = folder, subfolder = paste0(method, "/condition_", index_sim), results_name = paste0("reference_draws_sample_", z, ".rds"))  
                }    
            }
        } else { # otherwise use DMH for P = 24
            # reference draws using DMH 
            draws <- foreach(z = 1:n_datasets, .packages = "dmrfit") %dopar% {
                # Unpack objects from replicate z
                sample_z <- samples_ls[[z]] # select replicate z
                data_z <- sample_z$sample[,1:P]
                pars_z <- sample_z$pmles
                ref_str_z <- sample_z$ref_str # reference structure index
                save_sigma2_z <- sample_z$save_sigma2
                n_categories_z <- truth_info[[ref_str_z]]$n_categories
                # Compute reference draws using DMH
                reference_draws <- dmrfit:::cpp_dmh_sampler(data = data_z, pars = pars_z, n_categories = n_categories_z, 
                                                                nsim = nsim, burnin = burnin, L = L_DMH, inner_sampler_n_iter = inner_gibbs_iter,
                                                                adaptive_stage_n_iter = adaptive_stage_n_iter, sigma2 = sigma2, 
                                                                thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, 
                                                                interactions_location = interactions_location, interactions_scale = interactions_scale,
                                                                verbose = verbose, progress = progress)
                # ---- Reference draws ----
                .save_and_move_object(obj = reference_draws, folder = folder, subfolder = paste0(method, "/condition_", index_sim), results_name = paste0("reference_draws_sample_", z, ".rds"))
                if(save_sigma2_z){
                    .save_and_move_object(obj = reference_draws$sigma2, folder = folder, subfolder = "DMH", results_name = paste0("dmh_sigma2_condition_", index_sim, "_sample_", z, ".rds"))
                }
            }
        }

        return(invisible()) # function ends here if method is "reference"
    }
    else { # load simulated data and ground truth (reference draws are loaded within the foreach loop for each sample under a specific method)
        samples_ls <- readRDS(paste0(folder, "simulate_data/simulate_data_condition_", index_sim, ".rds"))
        truth_info  <- readRDS(paste0(folder, "simulate_data/info_data_condition_", index_sim, ".rds"))
    }

    # ---- Second phase : <calibration methods> ----

    # empty list where to store posterior distribution for metrics calculation (x) and empty list for posterior estimates from each method (post)
    summary_ls <- list() # list of one list, object named after the method's name containing: posterior estimates, overlapping metrics, Savage-Dickey density ratio, ESS, elapsed time

    # ---- Pseudo + Post-Hoc methods (RM - GHW - MCH) ----
    if(method == "PPH"){
        # ---- Pseudo-posterior draws first, then Post hoc methods ----
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {

            metrics <- list() # initialize metrics list for each replicate z

            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample[,1:P] # select the data for replicate z
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Draws for pseudo-posterior ---
            draws_pp <- dmrfit:::cpp_pseudo_sampler(data = data_z, pars = pars_z, n_categories = n_categories_z, nsim = nsim, burnin = burnin, 
                                                    adaptive_stage_n_iter = adaptive_stage_n_iter, sigma2 = sigma2, thresholds_alpha = thresholds_alpha,
                                                    thresholds_beta = thresholds_beta, interactions_location = interactions_location, interactions_scale = interactions_scale, 
                                                    verbose = verbose, progress = progress)
            if(save_sigma2_z){
                .save_and_move_object(obj = draws_pp$sigma2, folder = folder, subfolder = method, results_name = paste0("pseudo_sigma2_condition_", index_sim, "_sample_", z, ".rds"))
            }
            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Metrics for pseudo-posterior draws ---
            metrics$Pseudo <- .compute_metrics(x = draws_pp, 
                                            reference = ref_z, 
                                            method = "Pseudo", 
                                            n_pars = n_pars_z, 
                                            n_thresholds = n_thresholds_z, 
                                            structure = structure_z, 
                                            P = P, 
                                            interactions_location = interactions_location, 
                                            interactions_scale = interactions_scale)

            # --- Helpers ---

            # --- Pseudo-posterior estimates (means) ---
            pmles_z <- apply(draws_pp$draws, 2, mean)

            # --- Utils for pseudo-posterior draws ---
            utils_z <- dmrfit:::dmrf_deriv(pars = pmles_z, 
                                            data = sample_z$sample, # include sufficient statistics in the data for pseudo-likelihood function
                                            P = P, 
                                            n_categories = n_categories_z, 
                                            with_prior = TRUE,
                                            ncores = 1L,
                                            thresholds_alpha = thresholds_alpha,
                                            thresholds_beta = thresholds_beta,
                                            interactions_location = interactions_location,
                                            interactions_scale = interactions_scale)

            # --- Cholesky of Hessian matrix ---
            chol_hessian_z <- chol(x = utils_z$hessian)

            # --- Post hoc methods: RM, GHW, MCH ---

            # ---- RM (Robbins-Monro) ----
            draws_rm <- dmrfit:::cpp_rm_correction(pars = pmles_z, draws = draws_pp$draws, data = data_z, n_categories = n_categories_z,
                                                        chol_hessian = chol_hessian_z, L = L, sampler_n_iter = inner_gibbs_iter, 
                                                        thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, 
                                                        interactions_location = interactions_location, interactions_scale = interactions_scale, 
                                                        rm_step_thresholds = rm_step_thresholds, rm_step_interactions = rm_step_interactions, 
                                                        rm_max_iter = rm_max_iter, tolerance = tolerance)

            # Adjust elapsed time estimate (including RM correction time)
            draws_rm$seconds_elapsed <- draws_pp$seconds_elapsed + draws_rm$seconds_elapsed

            # Metrics for RM draws
            metrics$RM <- .compute_metrics(x = draws_rm, 
                                            reference = ref_z, 
                                            method = "RM", 
                                            n_pars = n_pars_z, 
                                            n_thresholds = n_thresholds_z,  
                                            structure = structure_z, 
                                            P = P, 
                                            interactions_location = interactions_location, 
                                            interactions_scale = interactions_scale)
            
            # Free space
            rm(draws_rm)

            # ---- GHW (Godambe-Huber-White correction) ----
            draws_ghw <- dmrfit:::cpp_ghw_correction(pars = pmles_z ,draws = draws_pp$draws, chol_hessian = chol_hessian_z, GHW = utils_z$HW) # [[utils_z$HW may be utils_z$GHW]]

            # Adjust elapsed time estimate (including GHW correction time)
            draws_ghw$seconds_elapsed <- draws_pp$seconds_elapsed + draws_ghw$seconds_elapsed

            # Metrics for GHW draws
            metrics$GHW <- .compute_metrics(x = draws_ghw, 
                                            reference = ref_z, 
                                            method = "GHW", 
                                            n_pars = n_pars_z, 
                                            n_thresholds = n_thresholds_z, 
                                            structure = structure_z, 
                                            P = P, 
                                            interactions_location = interactions_location, 
                                            interactions_scale = interactions_scale)

            # Free space
            rm(draws_ghw)

            # ---- MCH (Monte Carlo Hessian) ----
            draws_mch <- dmrfit:::cpp_mc_hessian_correction(pars = pmles_z, draws = draws_pp$draws, data = data_z, chol_hessian = chol_hessian_z, 
                                                            n_categories = n_categories_z, L = L, sampler_n_iter = inner_gibbs_iter, 
                                                            thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, 
                                                            interactions_location = interactions_location, interactions_scale = interactions_scale)

            # Adjust elapsed time estimate (including MCH correction time)
            draws_mch$seconds_elapsed <- draws_pp$seconds_elapsed + draws_mch$seconds_elapsed

            # Metrics for MCH draws
            metrics$MCH <- .compute_metrics(x = draws_mch, 
                                            reference = ref_z, 
                                            method = "MCH", 
                                            n_pars = n_pars_z, 
                                            n_thresholds = n_thresholds_z, 
                                            structure = structure_z, 
                                            P = P, 
                                            interactions_location = interactions_location, 
                                            interactions_scale = interactions_scale)
           
            # Free space
            rm(draws_mch)

            return(metrics)
        }

        # --- Metrics for pseudo-posterior draws and post hoc methods altogether done ---
        cat(paste(Sys.time()," - Pseudo-posterior draws and helpers done. \n",sep=""))
 
    }

    # ---- CoRe (with Godambe-Huber-White correction) ----
    if(method == "CoRe"){

        # --- CoRe ---
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {
 
            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample[,1:P] # select the data for replicate z
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Compute utils for the current replicate z ---
            utils_z <- dmrfit:::dmrf_deriv(pars = pars_z, 
                                            data = sample_z$sample, # include sufficient statistics in the data for pseudo-likelihood function
                                            P = P, 
                                            n_categories = n_categories_z, 
                                            with_prior = TRUE,
                                            ncores = 1L,
                                            thresholds_alpha = thresholds_alpha,
                                            thresholds_beta = thresholds_beta,
                                            interactions_location = interactions_location,
                                            interactions_scale = interactions_scale)

            # --- Current scale = Cholesky of Hessian matrix ---
            current_scale_z <- chol(x = utils_z$hessian)

            # --- New scale = Godambe-Huber-White matrix ---
            new_scale_z <- t(chol(utils_z$HW))

            # --- Remove utils_z to free space ---
            rm(utils_z)

            # --- Compute CoRe draws ---
            draws_core <- dmrfit:::cpp_core_sampler(data = data_z,
                                                    pars = pars_z, 
                                                    n_categories = n_categories_z,
                                                    pmles = pars_z, 
                                                    current_scale = current_scale_z,   
                                                    new_scale = new_scale_z, 
                                                    nsim = nsim, 
                                                    burnin = burnin, 
                                                    adaptive_stage_n_iter = adaptive_stage_n_iter, 
                                                    sigma2 = sigma2,
                                                    thresholds_alpha = thresholds_alpha, 
                                                    thresholds_beta = thresholds_beta, 
                                                    interactions_location = interactions_location, 
                                                    interactions_scale = interactions_scale,
                                                    verbose = verbose,
                                                    progress = progress)
            
            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Compute metrics ---
            metrics <- .compute_metrics(x = draws_core, 
                                        reference = ref_z, 
                                        method = "CoRe", 
                                        n_pars = n_pars_z, 
                                        n_thresholds = n_thresholds_z, 
                                        structure = structure_z, 
                                        P = P, 
                                        interactions_location = interactions_location, 
                                        interactions_scale = interactions_scale)

            return(metrics)
        }
        # --- Metrics for CoRe draws done ---
        cat(paste(Sys.time()," - CoRe draws computed. \n",sep=""))
    }

    # --- CoRe-RM ---
    if(method == "CoRe-RM"){

        # --- CoRe-RM ---
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {
 
            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample[,1:P] # select the data for replicate z
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Compute utils for the current replicate z ---
            utils_z <- dmrfit:::dmrf_deriv(pars = pars_z, 
                                            data = sample_z$sample, # include sufficient statistics in the data for pseudo-likelihood function
                                            P = P, 
                                            n_categories = n_categories_z, 
                                            with_prior = TRUE,
                                            ncores = 1L,
                                            thresholds_alpha = thresholds_alpha,
                                            thresholds_beta = thresholds_beta,
                                            interactions_location = interactions_location,
                                            interactions_scale = interactions_scale)

            # --- Current scale = Cholesky of Hessian matrix ---
            current_scale_z <- chol(x = utils_z$hessian)

            # --- Remove utils_z to free space ---
            rm(utils_z)

            # --- New scale = Robbins-Monro-based covariance matrix ---
            robbins_monro_ls <- dmrfit:::cpp_compute_robbins_monro(data = data_z,
                                                                    pars_init = pars_z,
                                                                    n_categories = n_categories_z,
                                                                    thresholds_alpha = thresholds_alpha, 
                                                                    thresholds_beta = thresholds_beta, 
                                                                    interactions_location = interactions_location, 
                                                                    interactions_scale = interactions_scale,
                                                                    rm_step_thresholds = rm_step_thresholds, 
                                                                    rm_step_interactions = rm_step_interactions,
                                                                    L = L, 
                                                                    sampler_n_iter = inner_gibbs_iter,
                                                                    rm_max_iter = rm_max_iter, 
                                                                    tolerance = tolerance) 
            new_scale_z <- t(chol(qr.solve(-robbins_monro_ls$hessian)))

            # --- Compute CoRe-RM draws ---
            draws_core_rm <- dmrfit:::cpp_core_sampler(data = data_z,
                                                        pars = pars_z, 
                                                        n_categories = n_categories_z,
                                                        pmles = pars_z, 
                                                        current_scale = current_scale_z,   
                                                        new_scale = new_scale_z, 
                                                        nsim = nsim, 
                                                        burnin = burnin, 
                                                        adaptive_stage_n_iter = adaptive_stage_n_iter, 
                                                        sigma2 = sigma2,
                                                        thresholds_alpha = thresholds_alpha, 
                                                        thresholds_beta = thresholds_beta, 
                                                        interactions_location = interactions_location, 
                                                        interactions_scale = interactions_scale,
                                                        verbose = verbose,
                                                        progress = progress)
            
            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Compute metrics ---
            metrics <- .compute_metrics(x = draws_core_rm, 
                                        reference = ref_z, 
                                        method = "CoRe-RM", 
                                        n_pars = n_pars_z, 
                                        n_thresholds = n_thresholds_z,
                                        structure = structure_z, 
                                        P = P, 
                                        interactions_location = interactions_location, 
                                        interactions_scale = interactions_scale)
            return(metrics)
        }

        # --- Metrics for CoRe-RM draws done ---
        cat(paste(Sys.time()," -  CoRe-RM draws done. \n",sep=""))
    }

    # --- CoRe-MCH --- 
    if(method == "CoRe-MCH"){
        # --- CoRe-MCH ---
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {
 
            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample[,1:P] # select the data for replicate z
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Compute utils for the current replicate z ---
            utils_z <- dmrfit:::dmrf_deriv(pars = pars_z, 
                                            data = sample_z$sample, # include sufficient statistics in the data for pseudo-likelihood function
                                            P = P, 
                                            n_categories = n_categories_z, 
                                            with_prior = TRUE,
                                            ncores = 1L,
                                            thresholds_alpha = thresholds_alpha,
                                            thresholds_beta = thresholds_beta,
                                            interactions_location = interactions_location,
                                            interactions_scale = interactions_scale)

            # --- Current scale = Cholesky of Hessian matrix ---
            current_scale_z <- chol(x = utils_z$hessian)

            # --- Remove utils_z to free space ---
            rm(utils_z)

            # --- New scale = Covariance matrix based on the Monte-Carlo-estimated Hessian ---
            mc_hessian_z <- dmrfit:::cpp_compute_mc_hessian(data = data_z,
                                                            pars = pars_z,
                                                            n_categories = n_categories_z,
                                                            L = L,
                                                            sampler_n_iter = inner_gibbs_iter,
                                                            thresholds_alpha = thresholds_alpha, 
                                                            thresholds_beta = thresholds_beta, 
                                                            interactions_location = interactions_location, 
                                                            interactions_scale = interactions_scale)                            
            new_scale_z <- t(chol(qr.solve(-mc_hessian_z)))

            # --- Compute CoRe-MCH draws ---
            draws_core_mch <- dmrfit:::cpp_core_sampler(data = data_z,
                                                        pars = pars_z, 
                                                        n_categories = n_categories_z,
                                                        pmles = pars_z, 
                                                        current_scale = current_scale_z,   
                                                        new_scale = new_scale_z, 
                                                        nsim = nsim, 
                                                        burnin = burnin, 
                                                        adaptive_stage_n_iter = adaptive_stage_n_iter, 
                                                        sigma2 = sigma2,
                                                        thresholds_alpha = thresholds_alpha, 
                                                        thresholds_beta = thresholds_beta, 
                                                        interactions_location = interactions_location, 
                                                        interactions_scale = interactions_scale,
                                                        verbose = verbose,
                                                        progress = progress)
            
            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Compute metrics ---
            metrics <- .compute_metrics(x = draws_core_mch, 
                                        reference = ref_z, 
                                        method = "CoRe-MCH", 
                                        n_pars = n_pars_z, 
                                        n_thresholds = n_thresholds_z, 
                                        structure = structure_z, 
                                        P = P, 
                                        interactions_location = interactions_location, 
                                        interactions_scale = interactions_scale)

            return(metrics)
        }

        # --- Metrics for CoRe-MCH draws done ---
        cat(paste(Sys.time()," - CoRe-MCH draws done. \n",sep=""))
    }

    # --- AdaCoRe ----
    if(method == "AdaCoRe"){

        # --- AdaCoRe ---
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {
 
            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample # select the data for replicate z (including P*(P-1)/2 crossproducts columns)
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Compute AdaCoRe draws ---
            draws_adacore <- dmrfit:::cpp_adacore_sampler(data = data_z, # columns: the P variables and the P*(P-1)/2 crossproducts 
                                                            pars = pars_z, 
                                                            n_categories = n_categories_z,
                                                            pmles = pars_z, 
                                                            nsim = nsim, 
                                                            burnin = burnin, 
                                                            adaptive_stage_n_iter = adaptive_stage_n_iter, 
                                                            sigma2 = sigma2,
                                                            thresholds_alpha = thresholds_alpha, 
                                                            thresholds_beta = thresholds_beta, 
                                                            interactions_location = interactions_location, 
                                                            interactions_scale = interactions_scale,
                                                            verbose = verbose,
                                                            progress = progress)

            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Compute metrics ---
            metrics <- .compute_metrics(x = draws_adacore, 
                                        reference = ref_z, 
                                        method = "AdaCoRe", 
                                        n_pars = n_pars_z, 
                                        n_thresholds = n_thresholds_z, 
                                        structure = structure_z, 
                                        P = P, 
                                        interactions_location = interactions_location, 
                                        interactions_scale = interactions_scale)           

            return(metrics)
        }
        # --- Metrics for AdaCoRe draws done ---
        cat(paste(Sys.time()," - AdaCoRe draws done. \n",sep=""))
    }

    # DMH only for P != 12
    if(method == "DMH" && P != 12){

        # --- Compute DMH draws ---
        summary_ls <- foreach(z = seq_along(samples_ls), .packages = "dmrfit") %dopar% {

            # --- Unpack objects from replicate z ---
            sample_z <- samples_ls[[z]] # select replicate z
            data_z <- sample_z$sample[,1:P] # select the data for replicate z
            pars_z <- sample_z$pmles # pmles
            ref_str_z <- sample_z$ref_str # reference structure index
            n_categories_z <- truth_info[[ref_str_z]]$n_categories
            n_thresholds_z <- truth_info[[ref_str_z]]$n_thresholds
            n_pars_z <- truth_info[[ref_str_z]]$n_pars
            structure_z <- truth_info[[ref_str_z]]$structure[-c(1:n_thresholds_z)] # only the 0/1 mask of the interaction parameters
            save_sigma2_z <- sample_z$save_sigma2

            # --- Compute DMH draws ---
            draws_dmh <- dmrfit:::cpp_dmh_sampler(data = data_z, pars = pars_z, n_categories = n_categories_z, 
                                                            nsim = nsim, burnin = burnin, L = L_DMH, inner_sampler_n_iter = inner_gibbs_iter,
                                                            adaptive_stage_n_iter = adaptive_stage_n_iter, sigma2 = sigma2, 
                                                            thresholds_alpha = thresholds_alpha, thresholds_beta = thresholds_beta, 
                                                            interactions_location = interactions_location, interactions_scale = interactions_scale,
                                                            verbose = verbose, progress = progress)

            # --- Load reference draws for replicate z ---
            ref_z <- readRDS(paste0(folder, "reference/condition_", index_sim, "/reference_draws_sample_", z, ".rds")) # load reference draws for replicate z

            # --- Compute metrics ---
            metrics <- .compute_metrics(x = draws_dmh, 
                                        reference = ref_z, 
                                        method = "DMH", 
                                        n_pars = n_pars_z, 
                                        n_thresholds = n_thresholds_z, 
                                        structure = structure_z, 
                                        P = P, 
                                        interactions_location = interactions_location, 
                                        interactions_scale = interactions_scale)               

            return(metrics)
        }

        # --- Metrics for DMH draws done ---
        cat(paste(Sys.time()," - DMH draws done. \n",sep=""))
        
    } else if(method == "DMH" && P == 12){
        stop("For P = 12, use method = 'reference' to obtain exact and DMH draws. This is a current design choice that reflects the paper's narrative.")
    }
 
    # save results object [[NOTE HERE: probably will have to check if .save_and_move_object() works fine in the line below]]
    summary_out <- list(summary = summary_ls, master_seed = master_seed)
    .save_and_move_object(obj = summary_out, folder = folder, subfolder = method, results_name = paste0("summary_condition_",index_sim,".rds"))
    cat(paste(Sys.time()," - Results saved. \n",sep=""))

    # stop the parallel cluster
    stopCluster(cl)
}