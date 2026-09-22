#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

#ifndef PAPER_UTILS_H
#define PAPER_UTILS_H


arma::uword sample_index(arma::uword N);

Rcpp::List cpp_gibbs_sampler_omrf(
                        const arma::mat& mu,  // matrix of thresholds (size P x (m-1))
                        const arma::mat& sigma,  // matrix of interactions (size P x P)
                        const arma::uvec& n_categories, // vector of number of categories for each variable (size P)
                        const arma::uword& N, // number of samples to generate
                        const arma::uword& P, // number of variables
                        const arma::uword& iter, // number of iterations for the Gibbs sampler
                        Rcpp::Nullable<arma::imat> X_start  = R_NilValue, // must be a matrix of size N x P
                        bool save_iter = false); // if true, save the state of X at each iteration

arma::rowvec gibbs_single_draw(const arma::mat& mu,
                                      const arma::mat& sigma,
                                      const arma::uvec& n_categories,
                                      arma::rowvec x_start,  
                                      arma::uword P,
                                      arma::uword iter); 
                    
void thresholds_vec_to_mat(arma::mat &thresholds_mat,  
                                    const arma::vec &pars, 
                                    const arma::uvec &n_categories, 
                                    arma::uword P, 
                                    arma::uword H);

void interactions_vec_to_mat(arma::mat &interactions_mat,       
                                    const arma::vec &pars,
                                    const arma::uvec &lower_indices,   
                                    arma::uword P,
                                    arma::uword n_thresholds,
                                    arma::uword n_pars);

void sigma_lower_tri_indices(arma::umat& matrix_indices_sigma, arma::uword P);

void sigma_lower_tri_indices_with_diagonal(arma::umat& matrix_indices_sigma, arma::uword P);

void threshold_stat_index_map(arma::uvec& which_stats,
                                     arma::vec& category_stats,
                                     const arma::uvec &n_categories,
                                     arma::uword n_thresholds,
                                     arma::uword P);

void category_offset_map(arma::uvec& category_offsets,
                        const arma::uvec& n_categories,
                        arma::uword P);                                     

void sufficient_statistics_omrf(arma::vec& obs_stats,
                                        const arma::mat& data,
                                        const arma::uvec &n_categories,
                                        arma::uword P,
                                        arma::uword n_thresholds,
                                        arma::uword n_pars,
                                        double interaction_scale);                                 

arma::vec fishermala_propose(const arma::vec& current_pars,
                            const arma::vec& approximate_gradient_current,
                            const arma::mat& R_n,
                            double sigma2_R);

void fishermala_update_preconditioner(arma::mat& R_n,
                                            double log_a,
                                            const arma::vec& approximate_gradient_current,
                                            const arma::vec& approximate_gradient_proposed,
                                            const arma::mat& I_d,
                                            double lambda,
                                            arma::uword s,
                                            arma::uword adaptive_stage_n_iter);                                           
void fishermala_update_step_size(double& sigma2,
                                        double& sigma2_R,
                                        const arma::mat& R_n,
                                        arma::uword n_pars,
                                        double log_a,
                                        double target_ar,
                                        double learning_rate);                            

// Compute pseudo-posterior based gradient and logZ_ratio
struct PseudoGradient {
    arma::vec   gradient;         // gradient evaluated at 'pars' (length is n_pars)
    double      logZ_ratio;       // logarithm of the ratio between normalizing constants (logZ_current - logZ_proposed)
};

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
                                    double interactions_scale);

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
                                        double interactions_scale);

// Compute approximate gradient and logZ_ratio using Monte Carlo approximation with L samples from the state space (single-data sample approximation)
struct ApproximateGradient {
    arma::vec   gradient;         // gradient evaluated at 'pars' (length is n_pars)
    double      logZ_ratio;       // logarithm of the ratio between normalizing constants (logZ_current - logZ_proposed)
};

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
                                    double interactions_scale);

arma::mat cpp_compute_mc_hessian(const arma::mat &data,
                                const arma::vec &pars,
                                const arma::uvec &n_categories,
                                arma::uword L,
                                arma::uword sampler_n_iter,
                                double thresholds_alpha, 
                                double thresholds_beta, 
                                double interactions_location, 
                                double interactions_scale);

Rcpp::List cpp_compute_robbins_monro(const arma::mat &data,
                                    const arma::vec &pars,
                                    const arma::uvec &n_categories,
                                    double thresholds_alpha, 
                                    double thresholds_beta, 
                                    double interactions_location, 
                                    double interactions_scale,
                                    double rm_step_thresholds = 0.001, // default value is set to 0.001 according to Bouranis et al.
                                    double rm_step_interactions = 0.001, // default value is set to 0.001 according to Bouranis et al.
                                    arma::uword L = 1000, // number of simulated networks at each iteration of the RM-algorithm (used to approximate gradient and hessian from the correct model) --> matches DMH's value
                                    arma::uword sampler_n_iter = 1, // number of iterations of the binary MRF Gibbs sampler (suggested number of iteration is (#nodes)^2)  
                                    arma::uword rm_max_iter = 200, // max number of iterations
                                    double tolerance = 0.0001);                            

arma::mat cpp_build_permutations_stats(const arma::mat &permutations,
                                    arma::uword n_pars,
                                    arma::uword n_thresholds,
                                    arma::uvec n_categories);                                   
#endif
