#include <string>
#include <RcppArmadillo.h>
#include "paper-utils.h"

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// CoRe edge-selection sampler: spike-and-slab prior on the transformed parameters beta, the pseudo-likelihood
// evaluated at eta(beta) = eta_star + Ainv (beta - eta_star), the curvature-adjusted pseudo-likelihood (Chandler and
// Bate, 2007; Ribatet et al., 2012) with the CoRe anchor and rescaling matrix (Arena and Marsman, 2026). The anchor
// and the map come from the full model and are fixed, so one adjusted likelihood serves every graph. Moves and phases
// as in edge_selection_sampler() (paper-utils.cpp).
// @param X unique response patterns [N x P], categories coded 0, ..., C_p - 1
// @param freq frequency of each pattern [N]
// @param n_categories number of categories per node [P]
// @param pars initial parameters beta [n_pars], thresholds then interactions (all edges start present)
// @param eta_star anchor of the map: pseudo-posterior mode of the full model [n_pars]
// @param Ainv inverse rescaling matrix [n_pars x n_pars]. A = Gamma L', with L L' = -H (H the Hessian of the log
//        pseudo-likelihood at eta_star) and Gamma Gamma' = V the target covariance (likelihood only)
// @param nsim number of kept iterations (after thinning)
// @param burnin number of burn-in iterations (phases 1-3). Phase 3 (edge selection) starts at burnin / 2
// @param adaptive_stage_n_iter iterations of phase 1 (plain MALA, R_n = I). Must be smaller than burnin / 2
// @param thin thinning of the kept iterations
// @param sigma2 initial step size of the FisherMALA block
// @param edge_prior "beta-bernoulli" (inclusion probability ~ Beta(beta_alpha, beta_beta), integrated out) or
//        "bernoulli" (inclusion_prob)
// @param slab "cauchy" or "normal", with location interactions_location and scale interactions_scale
// @param threshold_prior "betaprime" (thresholds_alpha, thresholds_beta, on exp(beta_mu)) or "normal"
//        (threshold_mean [n_thresholds], threshold_sd)
// @param use_likelihood false switches the likelihood off (for checks: the sampler then targets the prior)
// @param progress whether to show a progress bar
// @return a list with the draws of beta [nsim x n_pars] (absent interactions are exactly zero), the
//         Rao-Blackwellized inclusion probabilities, the move counts and acceptance rates, and the run time
// [[Rcpp::export]]
Rcpp::List cpp_core_edge_sampler(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                                 const arma::vec& pars, const arma::vec& eta_star, const arma::mat& Ainv,
                                 arma::uword nsim, arma::uword burnin,
                                 arma::uword adaptive_stage_n_iter = 500, arma::uword thin = 1, double sigma2 = 0.1,
                                 std::string edge_prior = "beta-bernoulli", double inclusion_prob = 0.5,
                                 double beta_alpha = 1.0, double beta_beta = 1.0,
                                 std::string slab = "cauchy", double interactions_location = 0.0,
                                 double interactions_scale = 2.5,
                                 std::string threshold_prior = "betaprime", double thresholds_alpha = 0.5,
                                 double thresholds_beta = 0.5,
                                 Rcpp::Nullable<Rcpp::NumericVector> threshold_mean = R_NilValue,
                                 double threshold_sd = 1.0,
                                 bool use_likelihood = true, bool progress = true) {

    arma::vec thr_mean;
    if(threshold_mean.isNotNull()){
        thr_mean = Rcpp::as<arma::vec>(threshold_mean.get());
    }

    // --- Model with the CoRe map ---
    EdgeModel M = edge_make_model(X, freq, n_categories, true, eta_star, Ainv,
                                  thresholds_alpha, thresholds_beta, interactions_location, interactions_scale,
                                  slab, threshold_prior, thr_mean, threshold_sd,
                                  edge_prior, inclusion_prob, beta_alpha, beta_beta, use_likelihood);

    return edge_selection_sampler(M, pars, nsim, burnin, adaptive_stage_n_iter, thin, sigma2, progress);
}
