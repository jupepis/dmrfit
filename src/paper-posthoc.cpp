#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"
#include "paper-utils.h"
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// Post-hoc correction methods

// Post-hoc GHW (Godambe-Huber-White) correction
// @param pars           PMLEs, length n_pars, ordered [thresholds, interactions]
// @param draws          pseudo-posterior draws, n_draws x n_pars
// @param chol_hessian   upper-triangular Cholesky of the (prior-adjusted) negative Hessian -- Lt
// @param GHW            Godambe-Huber-White sandwich covariance (prior-adjusted)
// [[Rcpp::export]]
Rcpp::List cpp_ghw_correction(const arma::vec& pars,
                          const arma::mat& draws,
                          const arma::mat& chol_hessian,
                          const arma::mat& GHW) {
    arma::wall_clock timer;
    timer.tic();

    arma::mat Gamma = arma::chol(GHW, "lower");
    arma::uword n_draws = draws.n_rows;

    arma::mat adjusted_draws(n_draws, pars.n_elem);
    for (arma::uword s = 0; s < n_draws; ++s) {
        arma::vec v = Gamma * (chol_hessian * (draws.row(s).t() - pars));
        adjusted_draws.row(s) = (v + pars).t();
    }

    double elapsed = timer.toc();

    return Rcpp::List::create(
        Rcpp::Named("draws")           = adjusted_draws,
        Rcpp::Named("seconds_elapsed") = elapsed
    );
}


// Post-hoc Monte Carlo Hessian correction
// @param pars           PMLEs, length n_pars, ordered [thresholds, interactions]
// @param draws          pseudo-posterior draws, n_draws x n_pars
// @param data           observed data (needed by the MC Hessian approximation)
// @param chol_hessian   upper-triangular Cholesky of the (prior-adjusted) negative Hessian -- Lt
// @param n_categories   categories per node
// @param L              number of MC samples for the Hessian approximation
// @param sampler_n_iter number of iterations for the Gibbs sampler (default is 1)
// @param thresholds_alpha hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param thresholds_beta hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param interactions_location location parameter for the cauchy prior on the interactions (default is 0.0)
// @param interactions_scale scale parameter for the cauchy prior on the interactions (default is 2.5)
// @return a list with the adjusted draws
// [[Rcpp::export]]
Rcpp::List cpp_mc_hessian_correction(const arma::vec& pars,
                                 const arma::mat& draws,
                                 const arma::mat& data,
                                 const arma::mat& chol_hessian,
                                 const arma::uvec& n_categories,
                                 arma::uword L,
                                 arma::uword sampler_n_iter,
                                 double thresholds_alpha, 
                                 double thresholds_beta, 
                                 double interactions_location, 
                                 double interactions_scale) {
    arma::wall_clock timer;
    timer.tic();

    // Monte Carlo approximation of the Hessian at `pars`
    arma::mat hessian_mc = cpp_compute_mc_hessian(data, pars, n_categories, L, sampler_n_iter, 
                                              thresholds_alpha, thresholds_beta, 
                                              interactions_location, interactions_scale); // default sampler_n_iter = 1
    arma::mat hessian_mc_sym = (hessian_mc + hessian_mc.t()) * 0.5; // symmetrize the Hessian to avoid numerical issues (because MCMC-based)

    // Compute the Cholesky of the inverse of the negative Hessian
    arma::mat Gamma_mc = arma::chol(arma::inv_sympd(-hessian_mc_sym), "lower");

    arma::uword n_draws = draws.n_rows;
    arma::mat adjusted_draws(n_draws, pars.n_elem);
    for (arma::uword s = 0; s < n_draws; ++s) {
        arma::vec v = Gamma_mc * (chol_hessian * (draws.row(s).t() - pars));
        adjusted_draws.row(s) = (v + pars).t();
    }

    double elapsed = timer.toc();

    return Rcpp::List::create(
        Rcpp::Named("draws")           = adjusted_draws,
        Rcpp::Named("seconds_elapsed") = elapsed
    );
}

// Post-hoc Robbins-Monro correction (Bouranis et al.)
// @param pars           PMLEs, length n_pars, ordered [thresholds, interactions]
// @param draws          pseudo-posterior draws, n_draws x n_pars
// @param data           observed data (needed by the RM algorithm)
// @param n_categories   categories per node
// @param chol_hessian   upper-triangular Cholesky of the (prior-adjusted) negative Hessian -- Lt (= M in Bouranis et al.)
// @param L              number of MC samples used inside the RM algorithm
// @param sampler_n_iter number of iterations for the Gibbs sampler (default is 1)
// @param thresholds_alpha hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param thresholds_beta hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param interactions_location location parameter for the cauchy prior on the interactions (default is 0.0)
// @param interactions_scale scale parameter for the cauchy prior on the interactions (default is 2.5)
// @return a list with the adjusted draws
// [[Rcpp::export]]
Rcpp::List cpp_rm_correction(const arma::vec& pars,
                         const arma::mat& draws,
                         const arma::mat& data,
                         const arma::uvec& n_categories,
                         const arma::mat& chol_hessian,
                         arma::uword L,
                         arma::uword sampler_n_iter,
                         double thresholds_alpha, 
                         double thresholds_beta, 
                         double interactions_location, 
                         double interactions_scale,
                         double rm_step_thresholds = 0.01,
                         double rm_step_interactions = 0.001,
                         arma::uword rm_max_iter = 200,
                         double tolerance = 0.001) {
    arma::wall_clock timer;
    timer.tic();

    // Run Robbins-Monro stochastic approximation
    Rcpp::List rm = cpp_compute_robbins_monro(data, pars, n_categories, thresholds_alpha, thresholds_beta, 
                                            interactions_location, interactions_scale, rm_step_thresholds, 
                                            rm_step_interactions, L, sampler_n_iter, rm_max_iter, tolerance);
    arma::vec pars_star  = Rcpp::as<arma::vec>(rm["pars"]);
    arma::mat hessian_rm = Rcpp::as<arma::mat>(rm["hessian"]);
    arma::mat hessian_rm_sym = 0.5 * (hessian_rm + hessian_rm.t()); // symmetrize the Hessian to avoid numerical issues (because MCMC-based)

    // Compute the Cholesky of the inverse of the negative Hessian from RM
    arma::mat Gamma_rm = arma::chol(arma::inv_sympd(-hessian_rm_sym), "lower");

    arma::uword n_draws = draws.n_rows;
    arma::mat adjusted_draws(n_draws, pars.n_elem);
    for (arma::uword s = 0; s < n_draws; ++s) {
        arma::vec v = Gamma_rm * (chol_hessian * (draws.row(s).t() - pars));
        adjusted_draws.row(s) = (v + pars_star).t();   // re-center on pars_star, not pars
    }

    double elapsed = timer.toc();

    return Rcpp::List::create(
        Rcpp::Named("draws")           = adjusted_draws,
        Rcpp::Named("seconds_elapsed") = elapsed
    );
}