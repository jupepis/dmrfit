#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"
#include "paper-utils.h"
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// CoRe posterior sampler
// @param data matrix of size [N x P] with the observed data
// @param pars vector of size [n_pars] with the initial values of the parameters, in order [thresholds,interactions]
// @param n_categories vector of size [P] with the number of categories for each node
// @param pmles Pseudo-maximum likelihood estimates (pmles) 
// @param current_scale cholesky of variance and covariance matrix of pseudo-likelihood calculate ad the 'pmles' 
// @param new_scale transposed cholesky of inverse variance and covariance matrix - MC-hessian, Godambe-Huber-White or Robbins-Monro's
// @param nsim number of iterations after burnin
// @param burnin number of burnin iterations
// @param adaptive_stage_n_iter number of iterations for the initial adaptive stage (simple MALA)
// @param sigma2 initial value for the step size (this is adaptive and a good starting value is needed (defailt is 1.0), however, it is allowed to vary over iterations) --> for OMRF is set to 0.01 (still needed?)
// @param thresholds_alpha hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param thresholds_beta hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param interactions_location location parameter for the cauchy prior on the interactions (default is 0.0)
// @param interactions_scale scale parameter for the cauchy prior on the interactions (default is 2.5)
// @param verbose whether to print out the progress of each iteration (default is false). Also, verbose = true only if progress = false
// @param progress whether to show a progress bar (default is true)
// @return a list with the draws
// [[Rcpp::export]]
Rcpp::List cpp_core_sampler(const arma::mat &data,
                                const arma::vec &pars, 
                                const arma::uvec &n_categories,
                                const arma::vec &pmles, 
                                const arma::mat &current_scale,   
                                const arma::mat &new_scale, 
                                arma::uword nsim, 
                                arma::uword burnin, 
                                arma::uword adaptive_stage_n_iter = 500, 
                                double sigma2 = 1.0,
                                double thresholds_alpha = 0.5, 
                                double thresholds_beta = 0.5, 
                                double interactions_location = 0.0, 
                                double interactions_scale = 2.5,
                                bool verbose = false,
                                bool progress = true) {

    // --- Number of observations, variables, categories, and parameters ---
    arma::uword P = data.n_cols;
    arma::uword H = arma::max(n_categories - 1);
    arma::uword n_thresholds = arma::accu(n_categories - 1);
    arma::uword n_pars = pars.n_elem;

    // --- Number of total iterations (including burnin and adaptive stage) ---
    nsim += burnin + adaptive_stage_n_iter;
    arma::uword past_burnin = burnin + adaptive_stage_n_iter;
    arma::uword n_keep = nsim - past_burnin;

    // --- Precomputing matrices for the coordinate rescaling step ---
    arma::mat A = new_scale * current_scale;
    arma::mat Ainv = arma::inv(A);
    arma::mat AinvT = Ainv.t();

    // --- Compute sufficient statistics for the observed data ---
    arma::vec obs_stats;
    sufficient_statistics_omrf(obs_stats, data, n_categories, P, n_thresholds, n_pars, 2.0);

    // --- Find unique rows of the data and their frequencies ---
    arma::mat unique_data;
    arma::vec frequency;
    get_data_unique(unique_data, frequency, data);
    arma::uword N_unique = unique_data.n_rows; // For the pseudo posterior we need this sample size which is the number of unique rows in the data

    // --- Preallocate draws and sigma2_vec ---
    arma::mat draws(n_pars, n_keep, arma::fill::zeros);
    arma::vec sigma2_vec(nsim, arma::fill::zeros);

    // --- Parameters of adaptive proposal (FisherMALA) ---
    arma::mat R_n(n_pars, n_pars, arma::fill::eye);
    const arma::mat I_d(n_pars, n_pars, arma::fill::eye);   // diagonal matrix used in the adaptation of R_n
    double target_ar        = 0.574;   // optimal acceptance rate for MALA-type proposals (vs. 0.234 for RWM)
    double lambda           = 10.0;    // default value, Algorithm 1 in Fisher-adaptive MALA
    double learning_rate    = 0.015;   // default adaptation learning rate
    double sigma2_R         = sigma2; // this is  sigma2 / ((1/n_pars) * trace(R_n * R_n.t())), so sigma2_R = sigma2 at initial values

    // --- Initialize counters ---
    arma::uword s           = 1,
                rejected    = 0,
                accepted    = 0;
                                    
    // --- Setting starting values ---
    arma::vec beta_current = (A * (pars - pmles)) + pmles; // CoRe space 
    arma::vec eta_current = pars; // Original space
    sigma2_vec(s-1) = sigma2;

    // --- Initialize structures used by compute_pseudo_gradient() ---

    // Lower triangular matrix indices excluding diagonal elements
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                           
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // Utility vector indicating the offsets for each category in the thresholds vector
    arma::uvec category_offsets;
    category_offset_map(category_offsets, n_categories, P);

    // --- Calculate pseudo gradient and normalizing constant at current parameters, then trasnform to CoRe space gradient ---
    PseudoGradient core_current = compute_pseudo_gradient(obs_stats, unique_data, frequency, eta_current, eta_current, 
                                                            n_categories, which_stats, lower_indices, category_offsets, P, N_unique, H, 
                                                            n_pars, n_thresholds, thresholds_alpha, thresholds_beta, 
                                                            interactions_location, interactions_scale);
    core_current.gradient = AinvT * core_current.gradient; // transform to CoRe space gradient

    // --- Start timer for the chain ---
    arma::wall_clock timer_chain;
    timer_chain.tic();

    // --- Progress bar setup (only if progress == true) ---
    arma::uword print_every = 0;
    Progress p(nsim, progress);
    if (progress) {
        print_every = std::max<arma::uword>(1, nsim / 100);
    }

    // --- Start sampling ---
    while(s < nsim){

        if (progress && Progress::check_abort()) {
            Rcpp::stop("Sampler aborted by user.");
        }

        // (1) Propose new parameters and transform to original space
        arma::vec beta_proposed = fishermala_propose(beta_current, core_current.gradient, R_n, sigma2_R); // Core space
        arma::vec eta_proposed = (Ainv * (beta_proposed - pmles)) + pmles; // Original space

        // (2) Pseudo gradient and logZ at new parameters
        PseudoGradient core_proposed = compute_pseudo_gradient(obs_stats, unique_data, frequency, eta_proposed, eta_current,
                                                                n_categories, which_stats, lower_indices, category_offsets, P, N_unique, H, 
                                                                n_pars, n_thresholds, thresholds_alpha, thresholds_beta, 
                                                                interactions_location, interactions_scale);
        core_proposed.gradient = AinvT * core_proposed.gradient; // transform to CoRe space gradient

        // (3) Compute acceptance ratio
        double log_a = fishermala_core_log_acceptance_ratio(beta_current, beta_proposed, eta_current, eta_proposed, 
                                                            core_current.gradient, core_proposed.gradient, R_n, obs_stats, 
                                                            core_proposed.logZ_ratio, sigma2_R, n_thresholds, n_pars, thresholds_alpha, 
                                                            thresholds_beta, interactions_location, interactions_scale);

        // (4)  Adaptive step for R_n (FisherMALA)
        fishermala_update_preconditioner(R_n, log_a, core_current.gradient, core_proposed.gradient, I_d, lambda, s, adaptive_stage_n_iter);

        // (5) Adapt and normalize step size sigma2 and sigma2_R (FisherMALA)
        fishermala_update_step_size(sigma2, sigma2_R, R_n, n_pars,log_a, target_ar, learning_rate);
        // save adaptive sigma2 parameter
        sigma2_vec(s) = sigma2;   

        // (6) MH step: accept or reject the proposed parameters
        double u = R::runif(0.0, 1.0);

        // if(u < std::exp(log_a)) then set beta_current = beta_proposed update counter s++
        if(u < std::exp(log_a)){
            beta_current                        = beta_proposed;
            eta_current                         = eta_proposed;
            core_current.gradient               = core_proposed.gradient;
            if(s >= past_burnin){
                draws.col(s - past_burnin)  = beta_proposed;
                accepted++;
            }
        }
        else{
             
            if(s >= past_burnin){
                draws.col(s - past_burnin) = beta_current;
                rejected++;
            }
        }

        // --- Print progress if verbose == true ---
        if (verbose && !progress) {
            Rcpp::Rcout << "(" << accepted << " of " << (accepted + rejected) << ") || iteration " << s
                        << " || acceptance = " << std::exp(log_a)
                        << " || sigma2 (" << sigma2 << ") with (de)increment = "
                        << learning_rate * (std::exp(log_a) - target_ar) << "\n";
        }

        if (progress && s % print_every == 0) p.increment(print_every); // Update progress bar

        s++; // Increment counter
    }

    if (progress) {
        Rcpp::Rcout << "\n";
    }

    double time = timer_chain.toc(); // time elapsed (in seconds)

    // --- Calculate acceptance probability --- s / (s + rejected)
    double acceptance_probability = static_cast<double>(accepted)/(static_cast<double>(accepted) + static_cast<double>(rejected));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("draws") = draws.t(),
        Rcpp::Named("acceptance") = acceptance_probability,
        Rcpp::Named("seconds_elapsed") = time,
        Rcpp::Named("sigma2") = sigma2_vec
    );

    return out;
}

