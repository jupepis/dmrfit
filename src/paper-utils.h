#include <cmath>
#include <string>
#include <vector>
#include <RcppArmadillo.h>

#ifndef PAPER_EDGE_UTILS_H
#define PAPER_EDGE_UTILS_H

// Shared machinery of the spike-and-slab edge-selection samplers (pseudo-likelihood and CoRe).
// Parameter order follows dmrfit: thresholds (node by node), then interactions in trimatl_ind order (column j outer,
// row i > j inner). The sampled vector is psi (eta for the pseudo-likelihood sampler, beta for CoRe). The
// pseudo-likelihood is evaluated at eta(psi), with eta(beta) = eta_star + Ainv (beta - eta_star) for CoRe.

// --- Model: data, layout, map and priors ---
struct EdgeModel {
    arma::mat X;              // N x P unique response patterns (categories 0, ..., C_p - 1)
    arma::vec f;              // frequency of each pattern
    arma::uvec ncat;          // number of categories per node
    arma::uvec thr_start;     // position of node p's first threshold in the parameter vector
    arma::uvec edge_i, edge_j;
    arma::uword P, N, H, n_thr, n_int, n_pars;
    bool core;                // true: CoRe map, eta(beta) = eta_star + Ainv (beta - eta_star)
    arma::vec eta_star;       // CoRe anchor
    arma::mat Ainv;           // inverse CoRe rescaling matrix
    bool slab_normal;         // slab: Cauchy(int_loc, int_scale) or Normal(int_loc, int_scale^2)
    bool thr_normal;          // thresholds: beta-prime(thr_alpha, thr_beta) on exp(mu) or Normal(thr_mean, thr_sd^2)
    double thr_alpha, thr_beta, int_loc, int_scale, thr_sd;
    arma::vec thr_mean;
    bool beta_bernoulli;      // edge prior: beta-Bernoulli(bb_a, bb_b) or Bernoulli(pi)
    double log_prior_odds;    // Bernoulli: log(pi / (1 - pi))
    double bb_a, bb_b;
    bool use_lik;             // false: the likelihood is switched off (prior-only check)
};

// --- State of the chain ---
struct EdgeState {
    arma::vec psi;      // sampled parameters
    arma::vec eta;      // point at which the pseudo-likelihood is evaluated
    arma::mat mu;       // P x H thresholds at eta
    arma::mat R;        // N x P rest scores, R(n, p) = sum_j x_nj eta_pj
    arma::uvec gamma;   // edge indicators
    arma::uword n_present;
};

// --- Direction of a one-dimensional move in psi (change of eta, thresholds and rest scores per unit step) ---
struct EdgeDirection {
    std::vector<arma::uword> nodes; // nodes whose full conditional changes
    arma::mat dmu;                  // P x H
    arma::mat dR;                   // N x P
    arma::vec deta;                 // n_pars
};

// --- Pseudo-likelihood along a direction: value (over the affected nodes) and first two derivatives ---
struct EdgeLineEval {
    double value;
    double d1;
    double d2;
};

// --- Log target and gradient of the FisherMALA block (gamma fixed) ---
struct BlockEval {
    double value;    // log pseudo-likelihood + log prior of the block's parameters
    arma::vec grad;  // gradient in the block's parameters
    arma::mat mu, R; // conditionals at the evaluated point
};

// --- FisherMALA (Titsias, 2024) with the preconditioner learned in phase 2 and frozen afterwards ---
struct BlockFisherMALA {
    arma::mat R_n, I_d;
    arma::mat M_full;   // R_n R_n', set when frozen
    bool frozen = false;
    double sigma2, sigma2_R, lambda = 10.0, target = 0.574, learning_rate = 0.015;
    arma::uword d;
};

// --- Preconditioner of a block: S (proposal covariance shape) and L with L L' = S ---
struct BlockPrecond {
    arma::mat S;
    arma::mat L;
};

// --- Adaptation of the one-dimensional within moves ---
struct EdgeAdaptation {
    arma::vec log_eps2;    // per-edge log step size (initialized from the local information)
    arma::uvec n_adapt;    // number of adaptation steps per edge
};

// --- Counters of one run (after burn-in) ---
struct EdgeCounters {
    arma::umat edge_counts;  // per edge: add try/acc, delete try/acc, within try/acc
    arma::vec rb_sum;        // Rao-Blackwellized inclusion: sum of gamma + (1 - 2 gamma) alpha
    arma::uword block_tried = 0, block_accepted = 0;
};

// model
EdgeModel edge_make_model(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                          bool core, const arma::vec& eta_star, const arma::mat& Ainv,
                          double thresholds_alpha, double thresholds_beta,
                          double interactions_location, double interactions_scale,
                          const std::string& slab, const std::string& threshold_prior,
                          const arma::vec& threshold_mean, double threshold_sd,
                          const std::string& edge_prior, double inclusion_prob, double beta_alpha, double beta_beta,
                          bool use_likelihood);

// priors
double edge_log_dnorm(double x, double mean, double var);
double edge_slab_log(const EdgeModel& M, double x);
double edge_slab_d1(const EdgeModel& M, double x);
double edge_slab_d2(const EdgeModel& M, double x);
double edge_slab_info(const EdgeModel& M);
double edge_threshold_log(const EdgeModel& M, arma::uword l, double x);
double edge_threshold_d1(const EdgeModel& M, arma::uword l, double x);
double edge_log_odds_add(const EdgeModel& M, arma::uword s);
double edge_log_odds_delete(const EdgeModel& M, arma::uword s);

// state and evaluations
void edge_refresh_state(const EdgeModel& M, EdgeState& S);
EdgeDirection edge_direction(const EdgeModel& M, arma::uword k);
EdgeLineEval edge_line_eval(const EdgeModel& M, const EdgeState& S, const EdgeDirection& D, double t);
void edge_apply_move(const EdgeModel& M, EdgeState& S, const EdgeDirection& D, double t);
BlockEval block_eval(const EdgeModel& M, const EdgeState& S, const arma::uvec& A, const arma::vec& delta);

// log acceptance ratios (one per move)
double edge_within_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& current, const EdgeLineEval& proposed,
                                        double psi_current, double psi_proposed,
                                        double mean_forward, double mean_backward, double var);
double edge_add_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& at_zero, const EdgeLineEval& proposed,
                                     double psi_proposed, double q_mean, double q_var, arma::uword n_present);
double edge_delete_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& current, const EdgeLineEval& at_zero,
                                        double psi_current, double q_mean, double q_var, arma::uword n_present);
double block_log_acceptance_ratio(const BlockFisherMALA& F, const BlockPrecond& P, const BlockEval& current,
                                  const BlockEval& proposed, const arma::vec& x, const arma::vec& y);

// moves
void edge_within_move(const EdgeModel& M, const EdgeState& S, const EdgeDirection& D, arma::uword k,
                      double psi_start, double& t, EdgeLineEval& current,
                      EdgeAdaptation& adapt, EdgeCounters& counts, bool adapting, bool keep_stats);
void edge_between_move(const EdgeModel& M, EdgeState& S, const EdgeDirection& D, arma::uword k,
                       double psi_start, double& t, EdgeLineEval& current,
                       EdgeCounters& counts, bool keep_stats);
void fishermala_block_move(const EdgeModel& M, EdgeState& S, BlockFisherMALA& F, const arma::uvec& A,
                           arma::uword s, arma::uword adaptive_stage_n_iter,
                           EdgeCounters& counts, bool adapting, bool keep_stats);

// FisherMALA block utilities (also used by the exact sampler)
void fishermala_block_init(BlockFisherMALA& F, arma::uword d, double sigma2);
void fishermala_block_freeze(BlockFisherMALA& F);
BlockPrecond fishermala_block_precond(BlockFisherMALA& F, const arma::uvec& A);
arma::vec fishermala_block_propose(const BlockFisherMALA& F, const BlockPrecond& P, const arma::vec& x,
                                   const arma::vec& g);
void fishermala_block_adapt(BlockFisherMALA& F, double log_a, const arma::vec& gx, const arma::vec& gy,
                            arma::uword s, arma::uword adaptive_stage_n_iter);
void edge_shuffle(arma::uvec& order);

// sampler (phases and iterations), shared by the pseudo-likelihood and CoRe samplers
Rcpp::List edge_selection_sampler(const EdgeModel& M, const arma::vec& pars,
                                  arma::uword nsim, arma::uword burnin, arma::uword adaptive_stage_n_iter,
                                  arma::uword thin, double sigma2, bool progress);

#endif
