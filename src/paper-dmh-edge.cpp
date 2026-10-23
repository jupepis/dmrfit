#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "paper-utils.h"
#include "paper-exact-edge.h"
#include "samplers-utils.h"

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// DMH-MoMS edge-selection sampler: the full likelihood of the OMRF, with the ratio of normalizing constants estimated
// from auxiliary draws as in dmrfit's DMH sampler (cpp_dmh_sampler), and the same spike-and-slab model, moves and
// phases as the other edge-selection samplers.
//
//   log L(eta) = S' eta - n log Z(eta).
//
// L auxiliary draws y_l are made at a parameter eta_0, each from a randomly chosen observation followed by a few Gibbs
// iterations, as in dmrfit's DMH. With the log-weights lw_l = T(y_l)' (eta - eta_0), the ratio
// Z(eta') / Z(eta) is estimated by sum_l exp(lw_l(eta')) / sum_l exp(lw_l(eta)), the auxiliary ratio of DMH averaged
// over the draws (the Monte Carlo Metropolis-Hastings estimator of Liang and Jin, 2013). The evaluations are those of
// the exact sampler with the draws in place of all response patterns.
//
// One iteration:
//   (1) edge pass with the current draws. Each move reweights them, and new draws are made at the current state when
//       the effective sample size of the weights falls below L / 2.
//   (2) FisherMALA block move as in dmrfit's DMH: new draws at the proposal give the gradient there and the ratio of
//       normalizing constants. If the move is accepted, they become the current draws.
//
// Per iteration the sampler records the smallest effective sample size of the weights in the edge pass and the number
// of new sets of draws made there.

// L draws at eta, each from an observed row followed by inner_sampler_n_iter Gibbs iterations (dmrfit's
// gibbs_single_draw)
arma::imat dmh_draws(const EdgeModel& M, const arma::mat& data, const arma::vec& eta,
                     arma::uword L, arma::uword inner_sampler_n_iter) {
    // thresholds [P x H] and symmetric interactions [P x P] at eta
    arma::mat mu(M.P, M.H, arma::fill::zeros);
    for(arma::uword p = 0; p < M.P; p++){
        for(arma::uword h = 1; h < M.ncat(p); h++){
            mu(p, h - 1) = eta(M.thr_start(p) + h - 1);
        }
    }
    arma::mat sigma(M.P, M.P, arma::fill::zeros);
    for(arma::uword k = 0; k < M.n_int; k++){
        sigma(M.edge_i(k), M.edge_j(k)) = eta(M.n_thr + k);
        sigma(M.edge_j(k), M.edge_i(k)) = eta(M.n_thr + k);
    }
    arma::imat Y(L, M.P);
    for(arma::uword l = 0; l < L; l++){
        // starting point: a randomly chosen observation (dmrfit's sample_index, R's random number generator)
        arma::uword v = sample_index(data.n_rows);
        arma::rowvec x = gibbs_single_draw(mu, sigma, M.ncat, data.row(v), M.P, inner_sampler_n_iter);
        for(arma::uword p = 0; p < M.P; p++){
            Y(l, p) = static_cast<int>(x(p));
        }
    }
    return Y;
}

// new draws at the current state, after which the log-weights are zero (eta_0 = psi)
void dmh_refresh_draws(const EdgeModel& M, ExactData& E, ExactState& St, arma::vec& eta_0, const arma::mat& data,
                       arma::uword L, arma::uword inner_sampler_n_iter) {
    eta_0 = St.psi;
    E.Y = dmh_draws(M, data, eta_0, L, inner_sampler_n_iter);
    St.lw.zeros(L);
}

// In a move of edge k, a draw enters only through v = y_i y_j, an integer between 0 and (C_i - 1)(C_j - 1). The
// weights are summed by v once per edge, so that every evaluation along the line costs one term per value of v.
struct DMHLine {
    arma::vec W;        // sum of the weights exp(lw - mx) of the draws with each value of v
    arma::vec W2;       // sum of their squares (effective sample size)
    double mx;          // largest log-weight (shift)
    double s_k;         // sufficient statistic of the data for the edge
    double n;
    bool use_lik;
};

DMHLine dmh_line_prepare(const EdgeModel& M, const ExactData& E, const ExactState& St, arma::uword k,
                         arma::uword n_values) {
    DMHLine Lk;
    Lk.W.zeros(n_values);
    Lk.W2.zeros(n_values);
    Lk.mx = St.lw.max();
    Lk.s_k = E.S(M.n_thr + k);
    Lk.n = E.n;
    Lk.use_lik = M.use_lik;
    arma::uword i = M.edge_i(k), j = M.edge_j(k);
    for(arma::uword l = 0; l < E.Y.n_rows; l++){
        arma::uword v = static_cast<arma::uword>(E.Y(l, i) * E.Y(l, j));
        double w = std::exp(St.lw(l) - Lk.mx);
        Lk.W(v) += w;
        Lk.W2(v) += w * w;
    }
    return Lk;
}

// log L at psi + t e_k relative to psi, with its first two derivatives in t (as exact_line_eval(), t is the
// displacement of psi_k from its value in the state, so the derivatives in t are the derivatives in psi_k)
EdgeLineEval dmh_line_eval(const DMHLine& Lk, double t) {
    EdgeLineEval out = {0.0, 0.0, 0.0};
    if(!Lk.use_lik){
        return out;
    }
    // log-sum-exp over the values of v
    double am = -arma::datum::inf;
    for(arma::uword v = 0; v < Lk.W.n_elem; v++){
        if(Lk.W(v) > 0.0){
            am = std::max(am, std::log(Lk.W(v)) + t * static_cast<double>(v));
        }
    }
    double Z = 0.0, e1 = 0.0, e2 = 0.0;
    for(arma::uword v = 0; v < Lk.W.n_elem; v++){
        if(Lk.W(v) > 0.0){
            double dv = static_cast<double>(v);
            double w = std::exp(std::log(Lk.W(v)) + t * dv - am);
            Z += w;
            e1 += w * dv;
            e2 += w * dv * dv;
        }
    }
    e1 /= Z;
    e2 /= Z;
    out.value = Lk.s_k * t - Lk.n * (Lk.mx + am + std::log(Z));
    out.d1 = Lk.s_k - Lk.n * e1;
    out.d2 = -Lk.n * (e2 - e1 * e1);
    return out;
}

// effective sample size of the weights after a move of t, (sum w)^2 / sum w^2
double dmh_line_ess(const DMHLine& Lk, double t) {
    double sum_w = 0.0, sum_w2 = 0.0;
    for(arma::uword v = 0; v < Lk.W.n_elem; v++){
        double dv = static_cast<double>(v);
        sum_w += Lk.W(v) * std::exp(t * dv);
        sum_w2 += Lk.W2(v) * std::exp(2.0 * t * dv);
    }
    return sum_w * sum_w / sum_w2;
}

// log-weights after a move of t of edge k: only the statistic y_i y_j of the edge enters
void dmh_edge_update_log_weights(const EdgeModel& M, const ExactData& E, ExactState& St, arma::uword k, double t) {
    if(t == 0.0){
        return;
    }
    arma::uword i = M.edge_i(k), j = M.edge_j(k);
    for(arma::uword l = 0; l < E.Y.n_rows; l++){
        St.lw(l) += t * static_cast<double>(E.Y(l, i) * E.Y(l, j));
    }
}

// within move on present edge k, as exact_edge_within_move() with the estimated likelihood
void dmh_edge_within_move(const EdgeModel& M, const DMHLine& Lk, arma::uword k, double psi_start, double& t,
                          EdgeLineEval& current, EdgeAdaptation& adapt, EdgeCounters& counts, bool adapting,
                          bool keep_stats) {
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
    EdgeLineEval proposed = dmh_line_eval(Lk, psi_proposed - psi_start);
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

// add/delete move of edge k, as exact_edge_between_move() with the estimated likelihood
void dmh_edge_between_move(const EdgeModel& M, ExactState& St, const DMHLine& Lk, arma::uword k, double psi_start,
                           double& t, EdgeLineEval& current, EdgeCounters& counts, bool keep_stats) {
    EdgeLineEval at_zero = (St.gamma(k) == 0) ? current : dmh_line_eval(Lk, -psi_start);
    double g0 = at_zero.d1 + edge_slab_d1(M, 0.0);
    // curvature at psi_k = 0, floored as in exact_edge_between_move()
    double I0 = -at_zero.d2 + std::max(-edge_slab_d2(M, 0.0), edge_slab_info(M));
    double q_mean = g0 / I0;
    double q_var = 1.0 / I0;

    if(St.gamma(k) == 0){
        // edge absent, add it: (0, 0) -> (1, psi_proposed). The forward move draws psi_proposed from q_k and the
        // reverse move (delete) is deterministic, so the ratio has -log q_k(psi_proposed)
        double psi_proposed = q_mean + std::sqrt(q_var) * R::norm_rand();
        EdgeLineEval proposed = dmh_line_eval(Lk, psi_proposed);
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

// FisherMALA block move as in dmrfit's DMH: the proposal uses the gradient at the current state from the current
// draws, and new draws at the proposal give the gradient there (for the reverse proposal) and the ratio of
// normalizing constants. If the move is accepted, the new draws become the current draws (eta_0 = proposal).
void dmh_fishermala_block_move(const EdgeModel& M, ExactData& E, ExactState& St, arma::vec& eta_0,
                               const arma::mat& data, arma::uword L, arma::uword inner_sampler_n_iter,
                               BlockFisherMALA& F, const arma::uvec& A, arma::uword s,
                               arma::uword adaptive_stage_n_iter, EdgeCounters& counts, bool adapting,
                               bool keep_stats) {
    BlockPrecond P = fishermala_block_precond(F, A);
    const arma::vec no_change(A.n_elem, arma::fill::zeros);
    ExactBlockEval current = exact_block_eval(M, E, St, A, no_change);
    arma::vec x = St.psi(A);
    arma::vec y = fishermala_block_propose(F, P, x, current.grad);
    arma::vec eta_proposed = St.psi;
    eta_proposed(A) = y;

    // new draws at the proposal, and the log-weights of the current state relative to it
    ExactData E_new;
    E_new.Y = dmh_draws(M, data, eta_proposed, L, inner_sampler_n_iter);
    E_new.S = E.S;
    E_new.n = E.n;
    ExactState St_new = St;
    St_new.lw = exact_log_weights(M, E_new, St.psi - eta_proposed);
    ExactBlockEval current_new = exact_block_eval(M, E_new, St_new, A, no_change);
    ExactBlockEval proposed = exact_block_eval(M, E_new, St_new, A, y - x);

    BlockEval cur_b, prop_b; // the shared acceptance ratio uses the value and gradient only
    cur_b.value = current_new.value;
    cur_b.grad = current.grad;
    prop_b.value = proposed.value;
    prop_b.grad = proposed.grad;
    double log_a = block_log_acceptance_ratio(F, P, cur_b, prop_b, x, y);
    if(keep_stats){
        counts.block_tried++;
    }
    if(R::unif_rand() < std::exp(log_a)){
        St.psi(A) = y;
        E.Y = E_new.Y;
        eta_0 = eta_proposed;
        St.lw.zeros(L);
        if(keep_stats){
            counts.block_accepted++;
        }
    }
    if(adapting){
        fishermala_block_adapt(F, log_a, current.grad, proposed.grad, s + 1, adaptive_stage_n_iter);
    }
}

// DMH-MoMS edge-selection sampler
// @param S sufficient statistics of the data [n_pars]: category counts of the thresholds, sums of x_i x_j of the
//        interactions
// @param n number of observations
// @param X unique observed response patterns (used for the parameter layout only)
// @param freq frequency of each pattern
// @param data all observations [n x P], categories 0, ..., C_p - 1: the starting points of the auxiliary draws
// @param n_categories number of categories per node [P]
// @param pars initial parameters [n_pars], thresholds then interactions (all edges start present)
// @param nsim, burnin, adaptive_stage_n_iter, thin, sigma2 as in cpp_pseudo_edge_sampler()
// @param edge_prior, inclusion_prob, beta_alpha, beta_beta, slab, interactions_location, interactions_scale,
//        threshold_prior, thresholds_alpha, thresholds_beta, threshold_mean, threshold_sd priors, as in
//        cpp_pseudo_edge_sampler()
// @param L number of auxiliary draws
// @param inner_sampler_n_iter number of Gibbs iterations per auxiliary draw
// @param use_likelihood false switches the likelihood off (for checks: the sampler then targets the prior)
// @param progress whether to show a progress bar
// @return a list with the draws [nsim x n_pars], the Rao-Blackwellized inclusion probabilities, the move counts and
//         acceptance rates, the smallest effective sample size of the weights and the new sets of draws per
//         iteration, and the run time
// [[Rcpp::export]]
Rcpp::List cpp_dmh_edge_sampler(const arma::vec& S, double n, const arma::mat& X, const arma::vec& freq,
                                const arma::mat& data, const arma::uvec& n_categories, const arma::vec& pars,
                                arma::uword nsim, arma::uword burnin, arma::uword adaptive_stage_n_iter = 500,
                                arma::uword thin = 1,
                                double sigma2 = 0.1,
                                std::string edge_prior = "beta-bernoulli", double inclusion_prob = 0.5,
                                double beta_alpha = 1.0, double beta_beta = 1.0,
                                std::string slab = "cauchy", double interactions_location = 0.0,
                                double interactions_scale = 2.5,
                                std::string threshold_prior = "betaprime", double thresholds_alpha = 0.5,
                                double thresholds_beta = 0.5,
                                Rcpp::Nullable<Rcpp::NumericVector> threshold_mean = R_NilValue,
                                double threshold_sd = 1.0, arma::uword L = 25000, arma::uword inner_sampler_n_iter = 5,
                                bool use_likelihood = true, bool progress = true) {

    arma::vec thr_mean;
    if(threshold_mean.isNotNull()){
        thr_mean = Rcpp::as<arma::vec>(threshold_mean.get());
    }

    // --- Model (layout and priors) and data ---
    EdgeModel M = edge_make_model(X, freq, n_categories, false, arma::vec(), arma::mat(),
                                  thresholds_alpha, thresholds_beta, interactions_location, interactions_scale,
                                  slab, threshold_prior, thr_mean, threshold_sd,
                                  edge_prior, inclusion_prob, beta_alpha, beta_beta, use_likelihood);
    if(S.n_elem != M.n_pars || pars.n_elem != M.n_pars){
        Rcpp::stop("S and pars must have length n_pars.");
    }
    if(L < 2 || inner_sampler_n_iter < 1){
        Rcpp::stop("L must be at least 2 and inner_sampler_n_iter at least 1.");
    }
    if(data.n_cols != M.P || data.n_rows < 1){
        Rcpp::stop("data must have one column per variable and at least one row.");
    }
    const arma::uword max_cat = n_categories.max() - 1;
    const arma::uword n_values = max_cat * max_cat + 1; // values of y_i y_j

    // --- Phases ---
    arma::uword selection_start = burnin / 2;
    arma::uword n_iter = burnin + nsim * thin;
    if(adaptive_stage_n_iter < 1 || thin < 1){
        Rcpp::stop("adaptive_stage_n_iter and thin must be at least 1.");
    }
    if(adaptive_stage_n_iter >= selection_start){
        Rcpp::stop("adaptive_stage_n_iter must be smaller than burnin / 2 (phase 2 would be empty).");
    }

    // --- Initial state: full network, draws at the initial parameters ---
    ExactState St;
    St.psi = pars;
    St.gamma.ones(M.n_int);
    St.n_present = M.n_int;
    arma::vec eta_0 = pars;
    ExactData E;
    E.S = S;
    E.n = n;
    E.Y = dmh_draws(M, data, eta_0, L, inner_sampler_n_iter);
    const arma::uvec all_pars = arma::regspace<arma::uvec>(0, M.n_pars - 1);
    const arma::uvec thresholds = arma::regspace<arma::uvec>(0, M.n_thr - 1);

    // --- Output and adaptation ---
    arma::mat draws(M.n_pars, nsim, arma::fill::zeros);
    arma::vec sigma2_trace(n_iter, arma::fill::zeros);
    arma::vec ess_min(n_iter, arma::fill::ones);   // smallest effective sample size of the weights, relative to L
    arma::uvec refreshes(n_iter, arma::fill::zeros); // new sets of draws made in the edge pass
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
        St.lw = exact_log_weights(M, E, St.psi - eta_0); // recompute from scratch (removes drift)

        // --- (1) edge pass (phase 3 and kept iterations): within move, then add/delete, edges in a random order ---
        if(selecting){
            edge_shuffle(edge_order);
            for(arma::uword idx = 0; idx < M.n_int; idx++){
                arma::uword k = edge_order(idx);
                double psi_start = St.psi(M.n_thr + k);
                double t = 0.0;
                DMHLine Lk = dmh_line_prepare(M, E, St, k, n_values);
                EdgeLineEval current = dmh_line_eval(Lk, 0.0);
                if(St.gamma(k) == 1){
                    dmh_edge_within_move(M, Lk, k, psi_start, t, current, adapt, counts, adapting, keep_stats);
                }
                dmh_edge_between_move(M, St, Lk, k, psi_start, t, current, counts, keep_stats);
                dmh_edge_update_log_weights(M, E, St, k, t);
                St.psi(M.n_thr + k) = (St.gamma(k) == 1) ? psi_start + t : 0.0;

                // new draws at the current state when the weights have degraded
                double ess = dmh_line_ess(Lk, t) / static_cast<double>(L);
                ess_min(s) = std::min(ess_min(s), ess);
                if(ess < 0.5){
                    dmh_refresh_draws(M, E, St, eta_0, data, L, inner_sampler_n_iter);
                    refreshes(s)++;
                }
            }
        }

        // --- (2) FisherMALA block: all parameters (phases 1-2), or thresholds and present interactions ---
        arma::uvec A = selecting ? arma::uvec(arma::join_cols(thresholds, arma::find(St.gamma == 1) + M.n_thr))
                                 : all_pars;
        dmh_fishermala_block_move(M, E, St, eta_0, data, L, inner_sampler_n_iter, F, A, s, adaptive_stage_n_iter,
                                  counts, adapting, keep_stats);
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
        Rcpp::Named("dmh") = Rcpp::List::create(
            Rcpp::Named("L") = L,
            Rcpp::Named("inner_sampler_n_iter") = inner_sampler_n_iter,
            Rcpp::Named("ess_min") = ess_min,
            Rcpp::Named("refreshes") = refreshes),
        Rcpp::Named("seconds_elapsed") = seconds);
}
