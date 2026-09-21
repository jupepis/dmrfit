#include <string>
#include <RcppArmadillo.h>
#include "priors.h"
#include "paper-utils.h"

// Draw a random row index in [0, N) using R's RNG
arma::uword sample_index(arma::uword N) {
    double u = R::runif(0.0, 1.0);
    arma::uword v = static_cast<arma::uword>(std::floor(u * static_cast<double>(N)));
    if (v >= N) v = N - 1;   // guard the u==1 edge case
    return v;
}

// Generate a sample of size N from an Ordinal MRF on P variables with parameters (mu, sigma)
// [[Rcpp::export]]
Rcpp::List cpp_gibbs_sampler_omrf(
                        const arma::mat& mu,  // matrix of thresholds (size P x (m-1))
                        const arma::mat& sigma,  // matrix of interactions (size P x P)
                        const arma::uvec& n_categories, // vector of number of categories for each variable (size P)
                        const arma::uword& N, // number of samples to generate
                        const arma::uword& P, // number of variables
                        const arma::uword& iter, // number of iterations for the Gibbs sampler
                        Rcpp::Nullable<arma::imat> X_start, // must be a matrix of size N x P
                        bool save_iter) // if true, save the state of X at each iteration
{
    // initialize X matrix
    arma::imat initial_X(N, P);
    if(X_start.isNull()){
        for(arma::uword p = 0; p < P; p++){
            arma::uword k = n_categories(p);              // categories 0..k-1
            for(arma::uword n = 0; n < N; n++){
                double u = R::runif(0.0, 1.0);
                arma::uword c = static_cast<arma::uword>(std::floor(u * k));
                if(c >= k) c = k - 1;                       // guard u==1 edge
                initial_X(n, p) = static_cast<int>(c);
            }
        }
    } else {
        initial_X = Rcpp::as<arma::imat>(X_start);
    }
    arma::mat X = arma::conv_to<arma::mat>::from(initial_X);

    // initialize output cube if save_iter is true
    arma::cube X_out;
    if(save_iter){
        X_out = arma::cube(N, P, iter, arma::fill::zeros);
    }

    for(arma::uword s = 0; s < iter; s++){
        if(save_iter) X_out.slice(s) = X; // save the current state of X if save_iter = true
        for(arma::uword p = 0; p < P; p++){
            arma::uword n_categories_p = n_categories(p);
            for(arma::uword n = 0; n < N; n++){
                arma::rowvec X_n = X.row(n);
                double score_p = arma::accu(X_n * sigma.col(p));
                // (1) create a vector of probabiity of the Xp for h = 1, ..., m-1
                arma::vec pr_p(n_categories_p,arma::fill::zeros);
                pr_p(0) = 1.0; // setting the first numerator of the probability vector to 1.0 (for Xp = 0)
                double denom = 1.0;
                arma::uword j;
                for(j = 1; j < n_categories_p; j++){
                    // numerator
                    pr_p(j) = std::exp(mu(p,j-1) + score_p*static_cast<double>(j)); 
                    // denominator 
                    denom += std::exp(mu(p,j-1) + score_p*static_cast<double>(j));
                }
                pr_p /= denom;

                // (2) generate the category for the p-th variable (use while loop)
                j = 0;
                arma::vec cdf_p = arma::cumsum(pr_p);
                double u = R::runif(0.0, 1.0); // arma::randu(); we use R to ensure reproducibility with R's RNG
                while(u > cdf_p(j)){
                   j++;
                }
                X(n,p) = static_cast<double>(j);
            }
        }
    }

    // return the output
    if(save_iter){
        return Rcpp::List::create(
            Rcpp::Named("X") = X,
            Rcpp::Named("X_iter") = X_out
        );
    }
    return Rcpp::List::create(
        Rcpp::Named("X") = X
    );
}

// Generate a single draw from an Ordinal MRF on P variables with parameters (mu, sigma) using Gibbs sampling
arma::rowvec gibbs_single_draw(const arma::mat& mu,
                                      const arma::mat& sigma,
                                      const arma::uvec& n_categories,
                                      arma::rowvec x_start,  
                                      arma::uword P,
                                      arma::uword iter) {
    arma::rowvec x = x_start;

    for (arma::uword s = 0; s < iter; s++) {
        for (arma::uword p = 0; p < P; p++) {
            arma::uword n_categories_p = n_categories(p);
            double score_p = arma::accu(x * sigma.col(p));
            // (1) create a vector of probabiity of the Xp for h = 1, ..., m-1
            arma::vec pr_p(n_categories_p, arma::fill::zeros);
            pr_p(0) = 1.0; // setting the first numerator of the probability vector to 1.0 (for Xp = 0)
            double denom = 1.0;
            arma::uword j;
            for (j = 1; j < n_categories_p; j++) {
                // numerator 
                pr_p(j) = std::exp(mu(p, j - 1) + score_p * static_cast<double>(j));
                // denominator
                denom  += std::exp(mu(p, j - 1) + score_p * static_cast<double>(j));
            }
            pr_p /= denom;

            // (2) generate the category for the p-th variable (use while loop)
            j = 0;
            arma::vec cdf_p = arma::cumsum(pr_p);
            double u = R::runif(0.0, 1.0);   // we use R to ensure reproducibility with R's RNG
            while (u > cdf_p(j)) {
                j++;
            }
            x(p) = static_cast<double>(j);
        }
    }
    return x;
}

// Convert thresholds vector to matrix
void thresholds_vec_to_mat(arma::mat &thresholds_mat,  // pre-sized P x H by caller, reused across calls
                                    const arma::vec &pars, 
                                    const arma::uvec &n_categories, 
                                    arma::uword P, 
                                    arma::uword H) {
    thresholds_mat.zeros();                                    
    arma::uword i = 0;
    for(arma::uword p = 0; p < P; p++){
        for(arma::uword l = 0; l < (n_categories(p)-1); l++){
            thresholds_mat(p,l) = pars(i); 
            i++;
        }
    }
}

// Convert interactions vector to matrix
void interactions_vec_to_mat(arma::mat &interactions_mat,
                                    const arma::vec &pars,
                                    const arma::uvec &lower_indices,
                                    arma::uword P,
                                    arma::uword n_thresholds,
                                    arma::uword n_pars) {
    interactions_mat.zeros();
    interactions_mat(lower_indices) = pars(arma::span(n_thresholds, n_pars - 1));
    for (arma::uword j = 0; j < P - 1; ++j)
        for (arma::uword i = j + 1; i < P; ++i)
            interactions_mat(j, i) = interactions_mat(i, j);   // mirror upper from lower, no new allocation
}

// Define utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects 
void sigma_lower_tri_indices(arma::umat& matrix_indices_sigma, arma::uword P) {
    matrix_indices_sigma.zeros(2, P * (P - 1) / 2);
    arma::uword l = 0; // resent index operator l
    for(arma::uword j = 0; j < (P-1); j++){ // column j
        for(arma::uword i = (j+1); i < P; i++){ // row i
            matrix_indices_sigma(0,l) = i;
            matrix_indices_sigma(1,l) = j;
            l++;
        }
    }
}


// Define utility [matrix 2 x (P*(P-1)/2 + P)] lower triangular matrix including diagonal elements
void sigma_lower_tri_indices_with_diagonal(arma::umat& matrix_indices_sigma, arma::uword P) {
    matrix_indices_sigma.zeros(2, P * (P - 1) / 2 + P);
    arma::uword l = 0; // resent index operator l
    for(arma::uword j = 0; j < P; j++){ // column j
        for(arma::uword i = j; i < P; i++){ // row i
            matrix_indices_sigma(0,l) = i;
            matrix_indices_sigma(1,l) = j;
            l++;
        }
    }
}

// Define utility vectors indicating which stats and which category (of length n_thresholds - used in the hessian computation)
void threshold_stat_index_map(arma::uvec& which_stats,
                                     arma::vec& category_stats,
                                     const arma::uvec &n_categories,
                                     arma::uword n_thresholds,
                                     arma::uword P) {
    which_stats.zeros(n_thresholds);
    category_stats.zeros(n_thresholds);

    arma::uword l = 0;
    for (arma::uword p = 0; p < P; ++p){
        for (arma::uword h = 1; h < n_categories(p); ++h) {
            which_stats(l)    = p;
            category_stats(l) = static_cast<double>(h);
            l++;
        }
    }

}

// Compute the cumulative offset into the threshold parameter block for each variable (i.e. category_offsets(p) = number of threshold parameters for variables 0..p (inclusive). 
// Used to locate variable p's threshold slots within the flat n_thresholds-length block (offset by category_offsets(p-1) for p > 0).
void category_offset_map(arma::uvec& category_offsets,
                         const arma::uvec& n_categories,
                         arma::uword P) {
    category_offsets.set_size(P);
    category_offsets(0) = n_categories(0) - 1;
    for (arma::uword p = 1; p < P; p++) {
        category_offsets(p) = category_offsets(p - 1) + n_categories(p) - 1;
    }
}

// Compute sufficient statistics for the exact likelihood (used in the cpp_exact_sampler function)
void sufficient_statistics_omrf(arma::vec& obs_stats,
                                const arma::mat& data,
                                const arma::uvec &n_categories,
                                arma::uword P,
                                arma::uword n_thresholds,
                                arma::uword n_pars,
                                double interaction_scale) {

    // Utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects 
    arma::umat matrix_indices_sigma;
    sigma_lower_tri_indices(matrix_indices_sigma, P);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                               
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // Calculate sufficient statistics for the observed data
    obs_stats.set_size(n_pars);
    for(arma::uword i = 0; i < n_pars; i++){
        if(i < n_thresholds){
            // calculate threshold statistic for the auxiliary data
            arma::uword which_stats_idx = which_stats(i);
            double category_stats_i = category_stats(i);
            obs_stats(i) = arma::accu(data.col(which_stats_idx) == category_stats_i);
        }
        else{
            // calculate interaction statistic for the auxiliary data
            arma::uvec ij = matrix_indices_sigma.col(i - n_thresholds);
            obs_stats(i) = interaction_scale * arma::dot(data.col(ij(0)), data.col(ij(1))); // interaction_scale = 1.0 for exact and dmh, interaction_scale = 2.0 for pseudo-likelihood
        }
    }
}

// FisherMALA propose 
arma::vec fishermala_propose(const arma::vec& pars,
                                const arma::vec& gradient,
                                const arma::mat& R_n,
                                double sigma2_R) {
    // (we use Rcpp for reproducibility, using same R's RNG)                    
    Rcpp::NumericVector eta_n_r = Rcpp::rnorm(pars.n_elem, 0.0, 1.0);
    arma::vec eta_n = Rcpp::as<arma::vec>(eta_n_r);
    return pars + (sigma2_R / 2.0) * (R_n * (R_n.t() * gradient)) + std::sqrt(sigma2_R) * (R_n * eta_n);                   
}

// FisherMALA update preconditioner R_n
void fishermala_update_preconditioner(arma::mat& R_n, 
                                                double log_a, 
                                                const arma::vec& current_gradient, 
                                                const arma::vec& proposed_gradient, 
                                                const arma::mat& I_d,
                                                double lambda, 
                                                arma::uword s, 
                                                arma::uword adaptive_stage_n_iter) {

        if((s - adaptive_stage_n_iter) == 0){
            double r_proposed = 1.0 / (1.0 + std::sqrt(lambda / (lambda + arma::dot(proposed_gradient,proposed_gradient))));
            R_n = (1.0 / std::sqrt(lambda)) * (I_d - r_proposed * (proposed_gradient * proposed_gradient.t())/(lambda + (arma::dot(proposed_gradient,proposed_gradient)))); // Proposition 4 Formula (12)
        }
        else if(s > adaptive_stage_n_iter){
            // Compute adaptation signal
            arma::vec adapt_signal = std::sqrt(std::exp(log_a)) * (proposed_gradient - current_gradient);
            arma::vec phi_proposed = R_n.t() * adapt_signal;
            double r_proposed =  1.0 / (1.0 + std::sqrt(1.0 / (1.0 + arma::dot(phi_proposed,phi_proposed))));
            R_n -= r_proposed * (((R_n * phi_proposed) * phi_proposed.t())/(1.0 + (arma::dot(phi_proposed,phi_proposed)))); // Proposition 4 Formula (13)
        }
}

// FisherMALA update step size sigma2 and sigma2_R
void fishermala_update_step_size(double& sigma2, 
                                        double& sigma2_R, 
                                        const arma::mat& R_n, 
                                        arma::uword n_pars, 
                                        double log_a, 
                                        double target_ar, 
                                        double learning_rate) {
    // Adapt step size ..
    sigma2 = sigma2 * (1.0 + learning_rate * (std::exp(log_a) - target_ar));
    // .. and normalize it
    sigma2_R = sigma2 / ((1.0 / static_cast<double>(n_pars)) * arma::accu(arma::square(R_n))); // tr(R_n * R_n.t()) = sum(square(R_n))
}

PseudoGradient compute_pseudo_gradient(
                                    const arma::vec &stats,
                                    const arma::mat &data, // dataset with unique rows of the original data
                                    const arma::vec &frequency, // frequency of each unique row in the data
                                    const arma::vec &pars,
                                    const arma::vec &current_pars,
                                    const arma::uvec &n_categories,
                                    const arma::uvec &which_stats,
                                    const arma::uvec &lower_indices,
                                    const arma::uvec &category_offsets,
                                    arma::uword P,
                                    arma::uword N, // dimension of the decoupled data
                                    arma::uword H,
                                    arma::uword n_pars,
                                    arma::uword n_thresholds,
                                    double thresholds_alpha, 
                                    double thresholds_beta, 
                                    double interactions_location, 
                                    double interactions_scale)
{

    PseudoGradient out;
    arma::vec gradient(n_pars,arma::fill::zeros);
    arma::vec gradient_interactions_n(P*(P-1)/2,arma::fill::zeros);
    double logZ_ratio = 0.0;

    // --- Save proposed and current thresholds from pars and current_pars respectively ---
    arma::vec thresholds = pars(arma::span(0, n_thresholds - 1));  // vector of thresholds parameters
    arma::vec thresholds_current = current_pars(arma::span(0, n_thresholds - 1)); 

    // --- Convert proposed and current interactions from vector to matrix form ---
    arma::mat interactions(P,P);
    interactions_vec_to_mat(interactions, pars, lower_indices, P, n_thresholds, n_pars);
    arma::mat interactions_current(P,P);
    interactions_vec_to_mat(interactions_current, current_pars, lower_indices, P, n_thresholds, n_pars);

    // --- Compute log-prior at current parameters ---
    arma::vec log_prior(n_pars,arma::fill::zeros);
    log_prior(arma::span(0, n_thresholds - 1)) = log_beta_prime_first_derivative(pars(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta);
    log_prior(arma::span(n_thresholds, n_pars - 1)) = log_dcauchy_first_derivative(pars(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale);

    // --- Add prior contribution and sufficient statistics to gradient ---
    gradient += log_prior + stats; 

    // --- Compute the pseudo gradient and logZ_ratio (looping over N samples, the observations in the data)  ---
    for(arma::uword n = 0; n < N; n++){
        // --- Select n-th person statistics ---
        arma::vec stats_n = data.row(n).t(); // stats {X1,X2,...,XP} for the n-th person 
        // --- For each p we have to compute the support (normalizing constant, denominator) for each item, that is ln[1+sum_h{exp(mu_h+h*sum_{j!=p}{x_jsigma_pj})}] ---
        double logZ_ratio_n = 0.0;
        arma::vec probs_n(n_thresholds,arma::fill::zeros); // vector of probabilities  P(Xp = h)
        for(arma::uword p = 0; p < P; p++){
            double denom_p = 1.0;
            double denom_p_current_pars = 1.0;
            arma::vec stats_excl_p = stats_n;
            stats_excl_p(p) = 0.0; // this could be avoided because the interaction matrix has 0.0 in the diagonal
            double xixj_sigma = arma::dot(stats_excl_p, interactions.col(p));
            double xixj_sigma_current = arma::dot(stats_excl_p, interactions_current.col(p));
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h - 1;
                if(p > 0){
                    index_threshold_Xp += category_offsets(p - 1); // [[IT WAS -->]] arma::accu(n_categories(arma::span(0, p - 1)) - 1);
                }
                // --- At the proposed pars: pars ---
                double success_event = arma::accu(thresholds(index_threshold_Xp) + static_cast<double>(h) * xixj_sigma);

                // --- Calculating denom_p, and probs_n numerator for proposed pars ---
                denom_p += std::exp(success_event);
                probs_n(index_threshold_Xp) = std::exp(success_event); // success probability for p-th variable

                // --- At the current pars: current_pars ---
                double success_event_current_pars = arma::accu(thresholds_current(index_threshold_Xp) + static_cast<double>(h) * xixj_sigma_current); 

                // --- Calculating denom_p for current pars ---
                denom_p_current_pars += std::exp(success_event_current_pars);

            }
            // --- Calculate pr(Xp = h) --- 
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h - 1;
                if(p > 0){
                    index_threshold_Xp += category_offsets(p - 1); // [[IT WAS -->]] arma::accu(n_categories(arma::span(0, p - 1)) - 1);
                }
                probs_n(index_threshold_Xp) /= denom_p; // divide probs_n for the calculated denom_p
            }
            logZ_ratio_n += (std::log(denom_p_current_pars) - std::log(denom_p)); // this is the difference between the logZ_current and logZ_proposed for the n-th person cumulated over all the P variables
        }
        logZ_ratio += frequency(n) * logZ_ratio_n; // this is the logZ_ration_n (for the n-th person) multiplied by the frequency of the n-th person in the data (because we are working with unique rows of the data)

        // --- Gradient for thresholds ---
        gradient(arma::span(0, n_thresholds - 1)) -= frequency(n) * probs_n; // subtract expected value of thresholds

        // --- Gradient for interactions ---

        // First compute expected value for each variable 
        arma::vec expected_X(P, arma::fill::zeros);
        for(arma::uword p = 0; p < P; p++){
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword idx = h - 1;
                if(p > 0) {
                    idx += category_offsets(p-1);
                }
                expected_X(p)    += static_cast<double>(h) * probs_n(idx);
            }
        }
        // Then calculate gradient
        arma::uword index_ij = 0;
        gradient_interactions_n.zeros(); // reset the gradient for interactions for the n-th person
        for(arma::uword j = 0; j < (P-1); j++){
            for(arma::uword i = (j+1); i < P; i++){
                gradient_interactions_n(index_ij) -= (stats_n(j) * expected_X(i) + stats_n(i) * expected_X(j));
                index_ij++;
            }
        }
        gradient(arma::span(n_thresholds, n_pars - 1)) += frequency(n) * gradient_interactions_n;
        
    }


    out.gradient = gradient;
    out.logZ_ratio = logZ_ratio;

    return out;
}


// FisherMALA CoRe log acceptance ratio
double fishermala_core_log_acceptance_ratio(const arma::vec& beta_current,
                                        const arma::vec& beta_proposed,
                                        const arma::vec& eta_current,
                                        const arma::vec& eta_proposed,
                                        const arma::vec& core_gradient_current,
                                        const arma::vec& core_gradient_proposed,
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

    // --- Calculate log-proposal (on CoRe scale) --- (fast way to calculate without inverting vcov matrix - Proposition 1 from Michalis Titsias)
    arma::vec step_current = beta_current - beta_proposed - ((sigma2_R / 4.0) * (R_n * (R_n.t() * core_gradient_proposed)));
    double h_eta_current = arma::dot(step_current, core_gradient_proposed);
    arma::vec step_proposed = beta_proposed - beta_current - ((sigma2_R / 4.0) * (R_n * (R_n.t() * core_gradient_current)));
    double h_eta_proposed = arma::dot(step_proposed, core_gradient_current);
    double log_proposal = 0.5 * (h_eta_current - h_eta_proposed);

    // --- Calculate log-likelihood difference between current and proposed parameters (on original scale) ---
    double loglik_data = arma::accu((eta_proposed - eta_current) % obs_stats) + logZ_ratio; // logZ_ratio is calcualted on the original scale (eta)

    // --- Calculate log-prior (on original scale) ---
    double log_prior = arma::accu(
                        log_beta_prime(eta_proposed(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta) - 
                        log_beta_prime(eta_current(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta)
                        ) + 
                    arma::accu(
                        log_dcauchy(eta_proposed(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale) - 
                        log_dcauchy(eta_current(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale)
                    );

    // --- Calculate log-acceptance ratio ---
    double log_a = log_prior + loglik_data + log_proposal;

    if(log_a > 0.0){ 
        log_a = 0.0; // This is equivalent to log min{1,acceptance_ratio}
    }

    return log_a;
}


ApproximateGradient compute_approximate_gradient(
                                    const arma::vec &stats,
                                    const arma::mat &data,
                                    const arma::vec &pars,
                                    const arma::vec &current_pars,
                                    const arma::uvec &n_categories,
                                    const arma::uvec &which_stats,
                                    const arma::uvec &lower_indices,
                                    const arma::vec &category_stats,
                                    arma::uword L,
                                    arma::uword P,
                                    arma::uword N,
                                    arma::uword H,
                                    arma::uword n_pars,
                                    arma::uword n_thresholds,
                                    arma::uword sampler_n_iter,
                                    double thresholds_alpha, 
                                    double thresholds_beta, 
                                    double interactions_location, 
                                    double interactions_scale)
{

    ApproximateGradient out;
    arma::rowvec sim_l;
    arma::vec expected_stats(n_pars,arma::fill::zeros);
    arma::vec stats_all_l(n_pars);
    double expected_Z_ratio = 0.0;

    // --- Convert thresholds and interactions from vector to matrix form ---
    arma::mat thresholds_mat(P,H);
    thresholds_vec_to_mat(thresholds_mat, pars, n_categories, P, H);
    arma::mat interactions_mat(P,P);
    interactions_vec_to_mat(interactions_mat, pars, lower_indices, P, n_thresholds, n_pars);

    // --- Compute log-prior at current parameters ---
    arma::vec log_prior(n_pars,arma::fill::zeros);
    log_prior(arma::span(0, n_thresholds - 1)) = log_beta_prime_first_derivative(pars(arma::span(0, n_thresholds - 1)), thresholds_alpha, thresholds_beta);
    log_prior(arma::span(n_thresholds, n_pars - 1)) = log_dcauchy_first_derivative(pars(arma::span(n_thresholds, n_pars - 1)), interactions_location, interactions_scale);

    // --- Compute the approximate gradient (expectation using Monte Carlo on only one sample) ---
    for(arma::uword l = 0; l < L ; l++){
        // --- Randomly select a row from the observed data matrix to use as the starting point for the Gibbs sampler ---
        arma::uword v = sample_index(N); // uses R's RNG via Rcpp functions

        // --- Run the Gibbs sampler to generate a sample from the state space at pars ---
        sim_l = gibbs_single_draw(thresholds_mat, interactions_mat, n_categories, data.row(v), P, sampler_n_iter); // cpp_gibbs_sampler_omrf() needs to be exported using a header; the function returns an Rcpp::List with a named "X" matrix [1 x P]

        // --- Empty vector for sufficient statistics ---
        stats_all_l.zeros();

        // --- Sufficient statistic for thresholds: sum of I(Xi == h) ---
        arma::uword i = 0;
        for (arma::uword p = 0; p < P; p++) {
            for (arma::uword h = 1; h < n_categories(p); h++) {
                arma::uword which_stats_i = which_stats(i);
                double category_stats_p = category_stats(i);
                if (sim_l(which_stats_i) == category_stats_p) {
                    stats_all_l(i) += 1.0;
                }
                i++;
            }
        }

        // --- Sufficient statistic for interactions: sum of crossproducts XiXj ---
        arma::uword idx = n_thresholds;
        for (arma::uword j = 0; j < (P - 1); j++) {
            for (arma::uword i = j + 1; i < P; i++) {
                stats_all_l(idx) = sim_l(i) * sim_l(j);
                idx++;
            }
        }

        expected_stats      += stats_all_l/static_cast<double>(L); 
        expected_Z_ratio    += std::exp(arma::accu(stats_all_l % (current_pars - pars)))/static_cast<double>(L);

    }

    out.gradient = (stats - (expected_stats*static_cast<double>(N))) + log_prior;
    out.logZ_ratio = static_cast<double>(N)*std::log(expected_Z_ratio);

    return out;
}



// Compute Hessian matrix using Monte Carlo approximation (single-sample-based)
// [[Rcpp::export]]
arma::mat cpp_compute_mc_hessian(const arma::mat &data,
                            const arma::vec &pars,
                            const arma::uvec &n_categories,
                            arma::uword L,
                            arma::uword sampler_n_iter,
                            double thresholds_alpha, 
                            double thresholds_beta, 
                            double interactions_location, 
                            double interactions_scale)
{

    arma::uword n_thresholds = arma::accu(n_categories-1);
    arma::uword H = arma::max(n_categories-1);
    arma::uword N = data.n_rows;
    arma::uword P = n_categories.n_elem;
    arma::uword n_pars = pars.n_elem;

    // --- Initialize structures ---

    // Lower triangular matrix indices excluding diagonal elements
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    // Lower triangular matrix indices including diagonal elements
    arma::uvec lower_indices_including_diagonal = arma::trimatl_ind(arma::size(P, P), 0);

    // Lower triangular matrix indices including diagonal elements (i,j) pairs
    arma::umat matrix_indices_lower_diag;
    sigma_lower_tri_indices_with_diagonal(matrix_indices_lower_diag,P);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                               
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // --- Convert thresholds and interactions from vector to matrix form ---
    arma::mat thresholds_mat(P,H);
    thresholds_vec_to_mat(thresholds_mat, pars, n_categories, P, H);
    arma::mat interactions_mat(P,P);
    interactions_vec_to_mat(interactions_mat, pars, lower_indices, P, n_thresholds, n_pars);

    // --- Initialize vectors and matrices ---
    arma::rowvec sim_l;
    arma::vec stats_all_l(n_pars);
    arma::mat hessian(n_pars,n_pars,arma::fill::zeros);
    arma::vec E_stats(n_pars,arma::fill::zeros);
    arma::vec E_Xi_XiXj(n_pars * (n_pars - 1) / 2 + n_pars,arma::fill::zeros);
    arma::vec lower_triangular_Xi_XiXj(n_pars * (n_pars - 1) / 2 + n_pars);  

    for(arma::uword l = 0; l < L; l++){
       // --- Randomly select a row from the observed data matrix to use as the starting point for the Gibbs sampler ---
        arma::uword v = sample_index(N); // uses R's RNG via Rcpp functions

        // --- Run the Gibbs sampler to generate a sample from the state space at pars ---
        sim_l = gibbs_single_draw(thresholds_mat, interactions_mat, n_categories, data.row(v), P, sampler_n_iter); 

        // --- Empty vector for sufficient statistics ---
        stats_all_l.zeros();

        // --- Sufficient statistic for thresholds: sum of I(Xi == h) ---
        arma::uword i = 0;
        for (arma::uword p = 0; p < P; p++) {
            for (arma::uword h = 1; h < n_categories(p); h++) {
                arma::uword which_stats_i = which_stats(i);
                double category_stats_p = category_stats(i);
                if (sim_l(which_stats_i) == category_stats_p) {
                    stats_all_l(i) += 1.0;
                }
                i++;
            }
        }

        // --- Sufficient statistic for interactions: sum of crossproducts XiXj ---
        arma::uword idx = n_thresholds;
        for (arma::uword j = 0; j < (P - 1); j++) {
            for (arma::uword i = j + 1; i < P; i++) {
                stats_all_l(idx) = sim_l(i) * sim_l(j);
                idx++;
            }
        }

        // update expectation vector
        E_stats += stats_all_l / static_cast<double>(L);

        // update expectation of the product of stats vector
        lower_triangular_Xi_XiXj.zeros(); 
        for (arma::uword k = 0; k < lower_triangular_Xi_XiXj.n_elem; ++k) {
            arma::uword i = matrix_indices_lower_diag(0,k);   
            arma::uword j = matrix_indices_lower_diag(1,k);
            lower_triangular_Xi_XiXj(k) = stats_all_l(i) * stats_all_l(j);
        }
        E_Xi_XiXj += lower_triangular_Xi_XiXj / static_cast<double>(L);
    }

    arma::mat E_Xi_E_XiXj = E_stats * E_stats.t();
    arma::vec lower_triangular_E_Xi_E_XiXj_next = E_Xi_E_XiXj(lower_indices_including_diagonal);
    arma::vec lower_triangular_H = - ((E_Xi_XiXj - lower_triangular_E_Xi_E_XiXj_next)*static_cast<double>(N)); 
    hessian(lower_indices_including_diagonal) = lower_triangular_H;
    hessian = arma::symmatl(hessian); // make hessian symmetric from its lower triangular

    // adding priors to hessian
    arma::vec hessian_prior = arma::join_cols(log_beta_prime_second_derivative(pars(arma::span(0,n_thresholds-1)), thresholds_alpha, thresholds_beta),
                                            log_dcauchy_second_derivative(pars(arma::span(n_thresholds,n_pars-1)), interactions_location, interactions_scale));
    arma::mat diag_prior = arma::diagmat(hessian_prior);
    hessian += diag_prior;
    
   return hessian;
}

// Robbins-Monro Ordinal MRF
// [[Rcpp::export]]
Rcpp::List cpp_compute_robbins_monro(const arma::mat &data,
                                    const arma::vec &pars_init,
                                    const arma::uvec &n_categories,
                                    double thresholds_alpha, 
                                    double thresholds_beta, 
                                    double interactions_location, 
                                    double interactions_scale,
                                    double rm_step_thresholds, // default value is set to 0.001 according to Bouranis et al.
                                    double rm_step_interactions, // default value is set to 0.001 according to Bouranis et al.
                                    arma::uword L, // number of simulated networks at each iteration of the RM-algorithm (used to approximate gradient and hessian from the correct model) --> matches DMH's value
                                    arma::uword sampler_n_iter, // number of iterations of the binary MRF Gibbs sampler (suggested number of iteration is (#nodes)^2)  
                                    arma::uword rm_max_iter, // max number of iterations
                                    double tolerance) { // tolerance value used by the stopping rule

    arma::uword n_thresholds = arma::accu(n_categories-1);
    arma::uword H = arma::max(n_categories-1);
    arma::uword N = data.n_rows;
    arma::uword P = n_categories.n_elem;
    arma::uword n_pars = pars_init.n_elem;
    Rcpp::List pars_iterations = Rcpp::List::create();
    arma::vec gradient(n_pars);
    arma::uword i = 0;

    // --- Initialize structures ---

    // Lower triangular matrix indices excluding diagonal elements
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                               
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // --- Compute sufficient statistics for the observed data ---
    arma::vec obs_stats;
    sufficient_statistics_omrf(obs_stats, data, n_categories, P, n_thresholds, n_pars, 1.0);


    double calculated_tolerance = 1000.0;
    arma::vec pars = pars_init;
    arma::vec old_pars = pars;
    while((calculated_tolerance > tolerance) && (i < rm_max_iter)){
        // (0) saving current pars 
        pars_iterations.push_back(pars);

        // (1) approximate gradient at current iteration pars 
        ApproximateGradient approximate_gradient = compute_approximate_gradient(obs_stats, data, pars, pars, n_categories,
                                                                            which_stats, lower_indices, category_stats, L, P, N, H, n_pars,
                                                                            n_thresholds, sampler_n_iter, thresholds_alpha, thresholds_beta,
                                                                            interactions_location, interactions_scale);
        gradient = approximate_gradient.gradient; // extract only the gradient
        double eps_i_thresholds = rm_step_thresholds / static_cast<double>(i + 1); // we sum +1 because i starts from 0
        double eps_i_interactions = rm_step_interactions / static_cast<double>(i + 1);

        // (2) next iteration parameters
        pars(arma::span(0,n_thresholds-1)) += eps_i_thresholds * gradient(arma::span(0, n_thresholds - 1));
        pars(arma::span(n_thresholds,n_pars-1)) += eps_i_interactions * gradient(arma::span(n_thresholds, n_pars - 1));

        // (3) Update iteration counter
        i++;

        // (4) Update tolerance and old_pars
        calculated_tolerance = arma::max(arma::abs(old_pars-pars)); 
        old_pars = pars;
    }

    // --- Approximate Hessian via MCMC at the converged parameters (pars) ---
    arma::mat hessian = cpp_compute_mc_hessian(data, pars, n_categories, L, sampler_n_iter, thresholds_alpha, 
                                            thresholds_beta, interactions_location, interactions_scale); 
    
   // create output list and return it
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("pars") = pars,
        Rcpp::Named("pars_iterations") = pars_iterations,
        Rcpp::Named("hessian") = hessian
        );

   return out;
}

// [[Rcpp::export]]
arma::mat cpp_build_permutations_stats(const arma::mat &permutations,
                                    arma::uword n_pars,
                                    arma::uword n_thresholds,
                                    arma::uvec n_categories) {

    arma::uword P = n_categories.n_elem; // number of variables
    arma::uword n_interactions = P * (P - 1) / 2;

    // --- Initialize structures used by compute_approximate_gradient() ---

    // Utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects
    arma::umat matrix_indices_sigma;
    sigma_lower_tri_indices(matrix_indices_sigma, P);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                               
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // --- Allocate matrix to store sufficient statistics for all permutations ---
    arma::mat X(permutations.n_rows,n_pars,arma::fill::zeros);

    // --- Loop over all permutations and compute sufficient statistics ---
    for(arma::uword r = 0; r < permutations.n_rows; r++){
        arma::rowvec perm_i = permutations.row(r); // no .t() copy needed

        // thresholds stats
        arma::uword l = 0;
        for(arma::uword p = 0; p < P; p++){
            for(arma::uword h = 1; h < n_categories(p); h++){
                if(perm_i(which_stats(l)) == category_stats(l)){
                    X(r, l) = 1.0;
                }
                l++;
            }
        }

        // interaction stats
        for(arma::uword k = 0; k < n_interactions; k++){
            arma::uword i = matrix_indices_sigma(0, k);
            arma::uword j = matrix_indices_sigma(1, k);
            X(r, n_thresholds + k) = perm_i(i) * perm_i(j);
        }
    }

    return X;
}

