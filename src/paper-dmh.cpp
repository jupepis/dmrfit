#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"
#include "paper-utils.h"
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// FisherMALA DMH log acceptance ratio
double fishermala_dmh_log_acceptance_ratio(const arma::vec& current_pars,
                                        const arma::vec& proposed_pars,
                                        const arma::vec& approximate_gradient_current,
                                        const arma::vec& approximate_gradient_proposed,
                                        const arma::mat& R_n,
                                        const arma::vec& obs_stats,
                                        double logZ_ratio,
                                        double sigma2_R,
                                        arma::uword n_thresholds,
                                        arma::uword n_pars,
                                        double thresholds_alpha,
                                        double thresholds_beta,
                                        double interactions_location,
                                        double interactions_scale) {

    // --- Calculate log-proposal --- (fast way to calculate without inverting vcov matrix - Proposition 1 from Michalis Titsias)
    arma::vec step_current = current_pars - proposed_pars - ((sigma2_R / 4.0) * (R_n * (R_n.t() * approximate_gradient_proposed)));
    double h_current_pars = arma::dot(step_current, approximate_gradient_proposed);
    arma::vec step_proposed = proposed_pars - current_pars - ((sigma2_R / 4.0) * (R_n * (R_n.t() * approximate_gradient_current)));
    double h_proposed_pars = arma::dot(step_proposed, approximate_gradient_current);
    double log_proposal = 0.5 * (h_current_pars - h_proposed_pars);

    // --- Calculate log-likelihood difference between current and proposed parameters ---
    double loglik_data = arma::accu((proposed_pars - current_pars) % obs_stats) + logZ_ratio; 

    // --- Calculate log-prior ---
    double log_prior = arma::accu(
                        log_beta_prime(proposed_pars(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta) - 
                        log_beta_prime(current_pars(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta)
                        ) + 
                        arma::accu(
                            log_dcauchy(proposed_pars(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale) - 
                            log_dcauchy(current_pars(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale)
                        );

    // --- Calculate log-acceptance ratio ---
    double log_a = log_prior + loglik_data + log_proposal;

    if(log_a > 0.0){ 
        log_a = 0.0; // This is equivalent to log min{1,acceptance_ratio}
    }

    return log_a;
}

// DMH posterior sampler
// @param data matrix of size [N x P] with the observed data
// @param pars vector of size [n_pars] with the initial values of the parameters, in order [thresholds,interactions]
// @param n_categories vector of size [P] with the number of categories for each node
// @param nsim number of iterations after burnin
// @param burnin number of burnin iterations
// @param L number of samples from the state space to approximate gradient and normalizing constant (default is 25,000)
// @param inner_sampler_n_iter number of inner Gibbs sampler iterations (default is 5)
// @param adaptive_stage_n_iter number of iterations for the initial adaptive stage (simple MALA)
// @param sigma2 initial value for the step size (this is adaptive and a good starting value is needed (defailt is 1.0), however, it is allowed to vary over iterations) --> for OMRF is set to 0.01 (still needed?)
// @param thresholds_alpha hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param thresholds_beta hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param interactions_location location parameter for the cauchy prior on the interactions (default is 0.0)
// @param interactions_scale scale parameter for the cauchy prior on the interactions (default is 2.5)
// @param verbose whether to print out the progress of each iteration (default is false). Also, verbose = true only if progress = false
// @param progress whether to display a progress bar (default is true)
// @return a list with the draws
// [[Rcpp::export]]
Rcpp::List cpp_dmh_sampler(const arma::mat &data,
                                const arma::vec &pars, 
                                const arma::uvec &n_categories, 
                                arma::uword nsim, 
                                arma::uword burnin, 
                                arma::uword L = 25000,
                                arma::uword inner_sampler_n_iter = 5,
                                arma::uword adaptive_stage_n_iter = 500, 
                                double sigma2 = 1.0,
                                double thresholds_alpha = 0.5, 
                                double thresholds_beta = 0.5, 
                                double interactions_location = 0.0, 
                                double interactions_scale = 2.5,
                                bool verbose = false,
                                bool progress = true) {

    // --- Number of observations, variables, max number of categories, categories, and parameters ---
    arma::uword N = data.n_rows;
    arma::uword P = data.n_cols;
    arma::uword H = arma::max(n_categories - 1);
    arma::uword n_thresholds = arma::accu(n_categories - 1);
    arma::uword n_pars = pars.n_elem;

    // --- Number of total iterations (including burnin and adaptive stage) ---
    nsim += burnin + adaptive_stage_n_iter;
    arma::uword past_burnin = burnin + adaptive_stage_n_iter;
    arma::uword n_keep = nsim - past_burnin;

    // --- Compute sufficient statistics for the observed data ---
    arma::vec obs_stats;
    sufficient_statistics_omrf(obs_stats, data, n_categories, P, n_thresholds, n_pars, 1.0);

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
    arma::vec current_pars = pars;
    sigma2_vec(s-1) = sigma2;

    // --- Initialize structures used by compute_approximate_gradient() ---

    // Lower triangular matrix indices excluding diagonal elements
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                               
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // --- Calculate approximate gradient and normalizing constant at current parameters ---
    ApproximateGradient approximate_current = compute_approximate_gradient(obs_stats, data, current_pars, current_pars, n_categories,
                                                                            which_stats, lower_indices, category_stats, L, P, N, H, n_pars,
                                                                            n_thresholds, inner_sampler_n_iter, thresholds_alpha, thresholds_beta,
                                                                            interactions_location, interactions_scale);

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

        // (1) Propose new parameters 
        arma::vec proposed_pars = fishermala_propose(current_pars, approximate_current.gradient, R_n, sigma2_R);

        // (2) Approximate gradient and logZ at new parameters
        ApproximateGradient approximate_proposed = compute_approximate_gradient(obs_stats, data, proposed_pars, current_pars, n_categories,
                                                                                which_stats, lower_indices, category_stats, L, P, N, H, n_pars,
                                                                                n_thresholds, inner_sampler_n_iter, thresholds_alpha, thresholds_beta,
                                                                                interactions_location, interactions_scale);

        // (3) Compute acceptance ratio
        double log_a = fishermala_dmh_log_acceptance_ratio(current_pars, proposed_pars, approximate_current.gradient, approximate_proposed.gradient,
                                                            R_n, obs_stats, approximate_proposed.logZ_ratio, sigma2_R, n_thresholds, n_pars,
                                                            thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);

        // (4)  Adaptive step for R_n (FisherMALA)
        fishermala_update_preconditioner(R_n, log_a, approximate_current.gradient, approximate_proposed.gradient, I_d, lambda, s, adaptive_stage_n_iter);

        // (5) Adapt and normalize step size sigma2 and sigma2_R (FisherMALA)
        fishermala_update_step_size(sigma2, sigma2_R, R_n, n_pars,log_a, target_ar, learning_rate);
        // save adaptive sigma2 parameter
        sigma2_vec(s) = sigma2;   

        // (6) MH step: accept or reject the proposed parameters
        double u = R::runif(0.0, 1.0);

        // if(u < std::exp(log_a)) then set current_pars = proposed_pars update counter s++
        if(u < std::exp(log_a)){
            current_pars                        = proposed_pars; 
            approximate_current.gradient        = approximate_proposed.gradient;
            if(s >= past_burnin){
                draws.col(s - past_burnin) = proposed_pars;
                accepted++;
            }
        }
        else{ 
            if(s >= past_burnin){
                draws.col(s - past_burnin) = current_pars;
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
    double acceptance_probability = static_cast<double>(accepted)/(static_cast<double>(accepted) + static_cast<double>(rejected)); //static_cast<double>(n_iter)/(static_cast<double>(n_iter)+static_cast<double>(rejected));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("draws") = draws.t(),
        Rcpp::Named("acceptance") = acceptance_probability,
        Rcpp::Named("seconds_elapsed") = time,
        Rcpp::Named("sigma2") = sigma2_vec
    );

    return out;
}