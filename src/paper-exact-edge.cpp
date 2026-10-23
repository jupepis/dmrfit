#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "paper-utils.h"
#include "paper-exact-edge.h"

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// Exact-likelihood edge-selection sampler: the full likelihood of the OMRF, with the normalizing constant computed
// by enumerating all response patterns, and the same spike-and-slab model, moves and phases as the pseudo-likelihood
// and CoRe samplers (paper-utils.cpp). Reference for small networks (number of patterns prod(n_categories)).
//
//   log L(eta) = S' eta - n log Z(eta),   log Z(eta) = log sum_y exp(T(y)' eta),
// with T(y) the sufficient statistics of pattern y (category indicators of the thresholds, products y_i y_j of the
// interactions) and S = sum_n T(x_n) those of the data. The acceptance ratios are the shared functions of
// paper-utils.cpp with log L in place of the log pseudo-likelihood. The evaluations below differ.
//
// The state keeps the vector of pattern log-weights lw = T eta (one entry per pattern), updated after each move and
// recomputed at the start of each iteration.

// The data and state types are in paper-exact-edge.h.

// pattern log-weights T eta (category indicators for the thresholds, y_i y_j for the interactions)
arma::vec exact_log_weights(const EdgeModel& M, const ExactData& E, const arma::vec& eta) {
    arma::uword n_patterns = E.Y.n_rows;
    arma::vec lw(n_patterns, arma::fill::zeros);
    for(arma::uword y = 0; y < n_patterns; y++){
        double v = 0.0;
        for(arma::uword p = 0; p < M.P; p++){
            int h = E.Y(y, p);
            if(h > 0){
                v += eta(M.thr_start(p) + h - 1);
            }
        }
        for(arma::uword k = 0; k < M.n_int; k++){
            double eta_k = eta(M.n_thr + k);
            if(eta_k != 0.0){
                v += eta_k * E.Y(y, M.edge_i(k)) * E.Y(y, M.edge_j(k));
            }
        }
        lw(y) = v;
    }
    return lw;
}

// interaction statistic y_i y_j of edge k for every pattern: the direction of a move in psi_k
arma::vec exact_edge_direction(const EdgeModel& M, const ExactData& E, arma::uword k) {
    arma::vec d = arma::conv_to<arma::vec>::from(E.Y.col(M.edge_i(k)) % E.Y.col(M.edge_j(k)));
    return d;
}

// log L at psi + t e_k relative to psi (S_k t - n log Z), with its first two derivatives in t. As in edge_line_eval(),
// t is the displacement of psi_k from its value in the state St (psi_k = psi_start + t), so the derivatives in t
// are the derivatives in psi_k.
EdgeLineEval exact_line_eval(const EdgeModel& M, const ExactData& E, const ExactState& St,
                                    const arma::vec& dlw, arma::uword k, double t) {
    EdgeLineEval out = {0.0, 0.0, 0.0};
    if(!M.use_lik){
        return out;
    }
    arma::vec lwt = St.lw + t * dlw;
    double mx = lwt.max();
    arma::vec w = arma::exp(lwt - mx);
    double Z = arma::accu(w);
    double e1 = arma::dot(w, dlw) / Z;
    double e2 = arma::dot(w, dlw % dlw) / Z;
    double s_k = E.S(M.n_thr + k);
    out.value = s_k * t - E.n * (mx + std::log(Z));
    out.d1 = s_k - E.n * e1;
    out.d2 = -E.n * (e2 - e1 * e1);
    return out;
}

// log target and gradient of the block A at psi_A + delta (gamma fixed). Also returns the new log-weights
ExactBlockEval exact_block_eval(const EdgeModel& M, const ExactData& E, const ExactState& St,
                                       const arma::uvec& A, const arma::vec& delta) {
    ExactBlockEval out;
    arma::vec psi_A = St.psi(A) + delta;
    if(arma::any(delta != 0.0)){
        arma::vec deta(M.n_pars, arma::fill::zeros);
        deta(A) = delta;
        out.lw = St.lw + exact_log_weights(M, E, deta);
    }
    else{
        out.lw = St.lw;
    }
    out.value = 0.0;
    out.grad.zeros(A.n_elem);
    if(M.use_lik){
        // expected sufficient statistics of the block under the model at psi_A + delta
        double mx = out.lw.max();
        arma::vec w = arma::exp(out.lw - mx);
        double Z = arma::accu(w);
        w /= Z;
        arma::vec expected(M.n_pars, arma::fill::zeros);
        for(arma::uword y = 0; y < E.Y.n_rows; y++){
            for(arma::uword p = 0; p < M.P; p++){
                int h = E.Y(y, p);
                if(h > 0){
                    expected(M.thr_start(p) + h - 1) += w(y);
                }
            }
        }
        for(arma::uword l = 0; l < A.n_elem; l++){
            if(A(l) >= M.n_thr){
                arma::uword k = A(l) - M.n_thr;
                expected(A(l)) = arma::dot(w, arma::conv_to<arma::vec>::from(E.Y.col(M.edge_i(k)) % E.Y.col(M.edge_j(k))));
            }
        }
        out.value = arma::dot(E.S(A), psi_A) - E.n * (mx + std::log(Z));
        out.grad = E.S(A) - E.n * expected(A);
    }
    for(arma::uword l = 0; l < A.n_elem; l++){
        if(A(l) < M.n_thr){
            out.value += edge_threshold_log(M, A(l), psi_A(l));
            out.grad(l) += edge_threshold_d1(M, A(l), psi_A(l));
        }
        else{
            out.value += edge_slab_log(M, psi_A(l));
            out.grad(l) += edge_slab_d1(M, psi_A(l));
        }
    }
    return out;
}

// within move on present edge k, as edge_within_move() with the exact likelihood
void exact_edge_within_move(const EdgeModel& M, const ExactData& E, const ExactState& St, const arma::vec& dlw,
                                   arma::uword k, double psi_start, double& t, EdgeLineEval& current,
                                   EdgeAdaptation& adapt, EdgeCounters& counts, bool adapting, bool keep_stats) {
    const double target_ar = 0.574, rm_exponent = 0.75;
    double psi_current = psi_start + t;
    double grad_current = current.d1 + edge_slab_d1(M, psi_current);
    if(std::isnan(adapt.log_eps2(k))){
        double info = -current.d2 + edge_slab_info(M);
        adapt.log_eps2(k) = std::log(1.0 / info);
    }
    double var = std::exp(adapt.log_eps2(k));
    double mean_forward = psi_current + 0.5 * var * grad_current;
    double psi_proposed = mean_forward + std::sqrt(var) * R::norm_rand();
    EdgeLineEval proposed = exact_line_eval(M, E, St, dlw, k, psi_proposed - psi_start);
    double grad_proposed = proposed.d1 + edge_slab_d1(M, psi_proposed);
    double mean_backward = psi_proposed + 0.5 * var * grad_proposed;
    double log_a = edge_within_log_acceptance_ratio(M, current, proposed, psi_current, psi_proposed,
                                                    mean_forward, mean_backward, var);
    double acc = std::exp(std::min(0.0, log_a));
    if(keep_stats){
        counts.edge_counts(k, 4)++;
    }
    if(R::unif_rand() < acc){
        t = psi_proposed - psi_start;
        current = proposed;
        if(keep_stats){
            counts.edge_counts(k, 5)++;
        }
    }
    if(adapting){
        adapt.log_eps2(k) += std::pow(static_cast<double>(++adapt.n_adapt(k)), -rm_exponent) * (acc - target_ar);
    }
}

// add/delete move of edge k, as edge_between_move() with the exact likelihood
void exact_edge_between_move(const EdgeModel& M, const ExactData& E, ExactState& St, const arma::vec& dlw,
                                    arma::uword k, double psi_start, double& t, EdgeLineEval& current,
                                    EdgeCounters& counts, bool keep_stats) {
    EdgeLineEval at_zero = (St.gamma(k) == 0) ? current : exact_line_eval(M, E, St, dlw, k, -psi_start);
    double g0 = at_zero.d1 + edge_slab_d1(M, 0.0);
    // curvature at psi_k = 0: the likelihood part -d2 (>= 0, the log likelihood is concave) plus the slab curvature,
    // floored at the slab information so that I0 > 0. With location 0 the floor is never active (Cauchy: 2 / s^2
    // against 1 / (2 s^2), normal: equal). It acts only for a Cauchy slab centered away from zero
    // (|location| > 0.68 s), whose curvature at 0 is small or negative.
    double I0 = -at_zero.d2 + std::max(-edge_slab_d2(M, 0.0), edge_slab_info(M));
    double q_mean = g0 / I0;
    double q_var = 1.0 / I0;

    if(St.gamma(k) == 0){
        // edge absent, add it: (0, 0) -> (1, psi_proposed). The forward move draws psi_proposed from q_k and the
        // reverse move (delete) is deterministic, so the ratio has -log q_k(psi_proposed)
        double psi_proposed = q_mean + std::sqrt(q_var) * R::norm_rand();
        EdgeLineEval proposed = exact_line_eval(M, E, St, dlw, k, psi_proposed);
        double log_a = edge_add_log_acceptance_ratio(M, at_zero, proposed, psi_proposed, q_mean, q_var, St.n_present);
        double acc = std::isnan(log_a) ? 0.0 : std::exp(std::min(0.0, log_a));
        if(keep_stats){
            counts.edge_counts(k, 0)++;
            counts.rb_sum(k) += acc;
        }
        if(R::unif_rand() < acc){
            St.gamma(k) = 1;
            St.n_present++;
            t = psi_proposed;
            current = proposed;
            if(keep_stats){
                counts.edge_counts(k, 1)++;
            }
        }
    }
    else{
        // edge present, delete it: (1, psi_current) -> (0, 0). The forward move is deterministic and the reverse move
        // (add) would draw psi_current from q_k, so the ratio has +log q_k(psi_current). q_k is built at psi_k = 0
        // with the other parameters unchanged, so it is the same density the reverse add would use
        double psi_current = psi_start + t;
        double log_a = edge_delete_log_acceptance_ratio(M, current, at_zero, psi_current, q_mean, q_var,
                                                        St.n_present);
        double acc = std::isnan(log_a) ? 0.0 : std::exp(std::min(0.0, log_a));
        if(keep_stats){
            counts.edge_counts(k, 2)++;
            counts.rb_sum(k) += 1.0 - acc;
        }
        if(R::unif_rand() < acc){
            St.gamma(k) = 0;
            St.n_present--;
            t = -psi_start;
            current = at_zero;
            if(keep_stats){
                counts.edge_counts(k, 3)++;
            }
        }
    }
}

// FisherMALA block move with the exact likelihood
void exact_fishermala_block_move(const EdgeModel& M, const ExactData& E, ExactState& St, BlockFisherMALA& F,
                                        const arma::uvec& A, arma::uword s, arma::uword adaptive_stage_n_iter,
                                        EdgeCounters& counts, bool adapting, bool keep_stats) {
    BlockPrecond P = fishermala_block_precond(F, A);
    ExactBlockEval current = exact_block_eval(M, E, St, A, arma::vec(A.n_elem, arma::fill::zeros));
    arma::vec x = St.psi(A);
    arma::vec y = fishermala_block_propose(F, P, x, current.grad);
    ExactBlockEval proposed = exact_block_eval(M, E, St, A, y - x);
    BlockEval cur_b, prop_b; // the shared acceptance ratio uses the value and gradient only
    cur_b.value = current.value;
    cur_b.grad = current.grad;
    prop_b.value = proposed.value;
    prop_b.grad = proposed.grad;
    double log_a = block_log_acceptance_ratio(F, P, cur_b, prop_b, x, y);
    if(keep_stats){
        counts.block_tried++;
    }
    if(R::unif_rand() < std::exp(log_a)){
        St.psi(A) = y;
        St.lw = proposed.lw;
        if(keep_stats){
            counts.block_accepted++;
        }
    }
    if(adapting){
        fishermala_block_adapt(F, log_a, current.grad, proposed.grad, s + 1, adaptive_stage_n_iter);
    }
}

// Exact-likelihood edge-selection sampler
// @param Y all response patterns [n_patterns x P], categories 0, ..., C_p - 1
// @param S sufficient statistics of the data [n_pars]: category counts of the thresholds, sums of x_i x_j of the
//        interactions
// @param n number of observations
// @param X unique observed response patterns (used for the parameter layout only)
// @param n_categories number of categories per node [P]
// @param pars initial parameters [n_pars], thresholds then interactions (all edges start present)
// @param nsim, burnin, adaptive_stage_n_iter, thin, sigma2 as in cpp_pseudo_edge_sampler()
// @param edge_prior, inclusion_prob, beta_alpha, beta_beta, slab, interactions_location, interactions_scale,
//        threshold_prior, thresholds_alpha, thresholds_beta, threshold_mean, threshold_sd priors, as in
//        cpp_pseudo_edge_sampler()
// @param use_likelihood false switches the likelihood off (for checks: the sampler then targets the prior)
// @param progress whether to show a progress bar
// @return a list with the draws [nsim x n_pars], the Rao-Blackwellized inclusion probabilities, the move counts and
//         acceptance rates, and the run time
// [[Rcpp::export]]
Rcpp::List cpp_exact_edge_sampler(const arma::imat& Y, const arma::vec& S, double n,
                                  const arma::mat& X, const arma::uvec& n_categories,
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

    // --- Model (layout and priors) and data of the exact likelihood ---
    EdgeModel M = edge_make_model(X, arma::ones<arma::vec>(X.n_rows), n_categories, false, arma::vec(), arma::mat(),
                                  thresholds_alpha, thresholds_beta, interactions_location, interactions_scale,
                                  slab, threshold_prior, thr_mean, threshold_sd,
                                  edge_prior, inclusion_prob, beta_alpha, beta_beta, use_likelihood);
    ExactData E;
    E.Y = Y;
    E.S = S;
    E.n = n;
    if(Y.n_cols != M.P){
        Rcpp::stop("Y must have one column per variable.");
    }
    if(S.n_elem != M.n_pars || pars.n_elem != M.n_pars){
        Rcpp::stop("S and pars must have length n_pars.");
    }

    // --- Phases ---
    arma::uword selection_start = burnin / 2;
    arma::uword n_iter = burnin + nsim * thin;
    if(adaptive_stage_n_iter < 1 || thin < 1){
        Rcpp::stop("adaptive_stage_n_iter and thin must be at least 1.");
    }
    if(adaptive_stage_n_iter >= selection_start){
        Rcpp::stop("adaptive_stage_n_iter must be smaller than burnin / 2 (phase 2 would be empty).");
    }

    // --- Initial state: full network ---
    ExactState St;
    St.psi = pars;
    St.gamma.ones(M.n_int);
    St.n_present = M.n_int;
    const arma::uvec all_pars = arma::regspace<arma::uvec>(0, M.n_pars - 1);
    const arma::uvec thresholds = arma::regspace<arma::uvec>(0, M.n_thr - 1);

    // --- Output and adaptation ---
    arma::mat draws(M.n_pars, nsim, arma::fill::zeros);
    arma::vec sigma2_trace(n_iter, arma::fill::zeros);
    EdgeCounters counts;
    counts.edge_counts.zeros(M.n_int, 6);
    counts.rb_sum.zeros(M.n_int);
    EdgeAdaptation adapt;
    adapt.log_eps2.set_size(M.n_int);
    adapt.log_eps2.fill(arma::datum::nan);
    adapt.n_adapt.zeros(M.n_int);
    BlockFisherMALA F;
    fishermala_block_init(F, M.n_pars, sigma2);
    arma::uvec edge_order = arma::regspace<arma::uvec>(0, M.n_int - 1);

    Progress bar(n_iter, progress);
    arma::uword print_every = std::max<arma::uword>(1, n_iter / 100);
    arma::wall_clock timer;
    timer.tic();

    for(arma::uword s = 0; s < n_iter; s++){
        Rcpp::checkUserInterrupt();
        bool adapting = s < burnin;
        bool keep_stats = s >= burnin;
        bool selecting = s >= selection_start;
        if(s == selection_start){
            fishermala_block_freeze(F);
        }
        St.lw = exact_log_weights(M, E, St.psi); // recompute from scratch (removes drift)

        // --- (1) edge pass (phase 3 and kept iterations): within move, then add/delete, edges in a random order ---
        if(selecting){
            edge_shuffle(edge_order);
            for(arma::uword idx = 0; idx < M.n_int; idx++){
                arma::uword k = edge_order(idx);
                arma::vec dlw = exact_edge_direction(M, E, k);
                double psi_start = St.psi(M.n_thr + k);
                double t = 0.0;
                EdgeLineEval current = exact_line_eval(M, E, St, dlw, k, 0.0);
                if(St.gamma(k) == 1){
                    exact_edge_within_move(M, E, St, dlw, k, psi_start, t, current, adapt, counts, adapting, keep_stats);
                }
                exact_edge_between_move(M, E, St, dlw, k, psi_start, t, current, counts, keep_stats);
                if(t != 0.0){
                    St.lw += t * dlw;
                }
                St.psi(M.n_thr + k) = (St.gamma(k) == 1) ? psi_start + t : 0.0;
            }
        }

        // --- (2) FisherMALA block: all parameters (phases 1-2), or thresholds and present interactions ---
        arma::uvec A = selecting ? arma::uvec(arma::join_cols(thresholds, arma::find(St.gamma == 1) + M.n_thr)) : all_pars;
        exact_fishermala_block_move(M, E, St, F, A, s, adaptive_stage_n_iter, counts, adapting, keep_stats);
        sigma2_trace(s) = F.sigma2;

        if(s >= burnin && (s - burnin) % thin == 0){
            draws.col((s - burnin) / thin) = St.psi;
        }
        if(progress && (s + 1) % print_every == 0){
            if(Progress::check_abort()){
                Rcpp::stop("Sampler aborted by user.");
            }
            bar.increment(print_every);
        }
    }

    double seconds = timer.toc();

    return Rcpp::List::create(
        Rcpp::Named("draws") = draws.t(),
        Rcpp::Named("edge_counts") = counts.edge_counts,
        Rcpp::Named("inclusion_rb") = counts.rb_sum / static_cast<double>(std::max<arma::uword>(1, nsim * thin)),
        Rcpp::Named("block_acceptance") = static_cast<double>(counts.block_accepted) /
                                          std::max<double>(1.0, counts.block_tried),
        Rcpp::Named("sigma2") = sigma2_trace,
        Rcpp::Named("edge_step") = arma::exp(adapt.log_eps2),
        Rcpp::Named("seconds_elapsed") = seconds);
}
