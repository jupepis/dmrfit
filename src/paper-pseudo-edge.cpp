#include <string>
#include <RcppArmadillo.h>
#include "paper-utils.h"

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// Pseudo-likelihood edge-selection sampler: spike-and-slab prior on the interactions, the pseudo-likelihood evaluated
// at the sampled parameters eta (the pseudo-posterior of bgms; Marsman et al., 2025), sampled with the moves and
// phases of edge_selection_sampler() (paper-utils.cpp).
// @param X unique response patterns [N x P], categories coded 0, ..., C_p - 1
// @param freq frequency of each pattern [N]
// @param n_categories number of categories per node [P]
// @param pars initial parameters [n_pars], thresholds then interactions (all edges start present)
// @param nsim number of kept iterations (after thinning)
// @param burnin number of burn-in iterations (phases 1-3). Phase 3 (edge selection) starts at burnin / 2
// @param adaptive_stage_n_iter iterations of phase 1 (plain MALA, R_n = I). Must be smaller than burnin / 2
// @param thin thinning of the kept iterations
// @param sigma2 initial step size of the FisherMALA block
// @param edge_prior "beta-bernoulli" (inclusion probability ~ Beta(beta_alpha, beta_beta), integrated out) or
//        "bernoulli" (inclusion_prob)
// @param slab "cauchy" or "normal", with location interactions_location and scale interactions_scale
// @param threshold_prior "betaprime" (thresholds_alpha, thresholds_beta, on exp(mu)) or "normal"
//        (threshold_mean [n_thresholds], threshold_sd)
// @param use_likelihood false switches the likelihood off (for checks: the sampler then targets the prior)
// @param progress whether to show a progress bar
// @return a list with the draws [nsim x n_pars] (absent interactions are exactly zero), the Rao-Blackwellized
//         inclusion probabilities, the move counts and acceptance rates, and the run time
// [[Rcpp::export]]
Rcpp::List cpp_pseudo_edge_sampler(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                                   const arma::vec& pars, arma::uword nsim, arma::uword burnin,
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

    // --- Model without the CoRe map: eta(psi) = psi ---
    EdgeModel M = edge_make_model(X, freq, n_categories, false, arma::vec(), arma::mat(),
                                  thresholds_alpha, thresholds_beta, interactions_location, interactions_scale,
                                  slab, threshold_prior, thr_mean, threshold_sd,
                                  edge_prior, inclusion_prob, beta_alpha, beta_beta, use_likelihood);

    return edge_selection_sampler(M, pars, nsim, burnin, adaptive_stage_n_iter, thin, sigma2, progress);
}
