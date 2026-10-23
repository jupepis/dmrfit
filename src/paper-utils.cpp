#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "paper-utils.h"

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// Spike-and-slab edge selection for ordinal MRFs: machinery shared by the pseudo-likelihood sampler
// (paper-pseudo-edge.cpp) and the CoRe sampler (paper-core-edge.cpp). For these two samplers the acceptance ratios
// are the same functions of the pseudo-likelihood values. The CoRe map only changes the point eta(beta) at which the
// pseudo-likelihood is evaluated, which is handled by edge_direction() and block_eval().
//
// Burn-in in three phases. Phase 1: full network, plain MALA on all parameters (R_n = I), sigma2 adapted. Phase 2:
// full network, FisherMALA on all parameters, learning the preconditioner M = R_n R_n'. Phase 3: edge selection
// with M frozen. Kept iterations follow. One edge-selection iteration:
//   (1) interactions in a new random order. For each edge k, a within move (1-D MALA) if the edge is present, then
//       an add/delete move (MoMS, Gottardo and Raftery, 2008, and Marsman et al., 2025) with a Newton proposal built at
//       psi_k = 0.
//   (2) one FisherMALA block move on the thresholds and the present interactions, gamma fixed, with the Schur
//       (conditional) block of the frozen M as preconditioner.

// Priors

double edge_softplus(double x) {
    return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

double edge_sigmoid(double x) {
    return (x >= 0.0) ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x));
}

// log density of mu when exp(mu) follows a beta-prime(a, b) distribution, and its derivative
double log_beta_prime_scalar(double x, double a, double b) {
    return a * x - (a + b) * edge_softplus(x) - R::lbeta(a, b);
}
double log_beta_prime_scalar_d1(double x, double a, double b) {
    return a - (a + b) * edge_sigmoid(x);
}

// Cauchy(l, g) log density (normalized: it enters the add/delete acceptance ratios) and its derivatives
double log_cauchy_scalar(double x, double l, double g) {
    double z = (x - l) / g;
    return -std::log(M_PI * g * (1.0 + z * z));
}
double log_cauchy_scalar_d1(double x, double l, double g) {
    double d = x - l;
    return -2.0 * d / (g * g + d * d);
}
double log_cauchy_scalar_d2(double x, double l, double g) {
    double d2 = (x - l) * (x - l), g2 = g * g;
    return 2.0 * (d2 - g2) / ((g2 + d2) * (g2 + d2));
}

double edge_log_dnorm(double x, double mean, double var) {
    return -0.5 * std::log(2.0 * M_PI * var) - 0.5 * (x - mean) * (x - mean) / var;
}

double edge_slab_log(const EdgeModel& M, double x) {
    return M.slab_normal ? R::dnorm(x, M.int_loc, M.int_scale, true) : log_cauchy_scalar(x, M.int_loc, M.int_scale);
}
double edge_slab_d1(const EdgeModel& M, double x) {
    return M.slab_normal ? -(x - M.int_loc) / (M.int_scale * M.int_scale) : log_cauchy_scalar_d1(x, M.int_loc, M.int_scale);
}
double edge_slab_d2(const EdgeModel& M, double x) {
    return M.slab_normal ? -1.0 / (M.int_scale * M.int_scale) : log_cauchy_scalar_d2(x, M.int_loc, M.int_scale);
}
// Fisher information of the slab location: the prior part of the within-move metric (always positive)
double edge_slab_info(const EdgeModel& M) {
    return M.slab_normal ? 1.0 / (M.int_scale * M.int_scale) : 1.0 / (2.0 * M.int_scale * M.int_scale);
}
double edge_threshold_log(const EdgeModel& M, arma::uword l, double x) {
    return M.thr_normal ? R::dnorm(x, M.thr_mean(l), M.thr_sd, true) : log_beta_prime_scalar(x, M.thr_alpha, M.thr_beta);
}
double edge_threshold_d1(const EdgeModel& M, arma::uword l, double x) {
    return M.thr_normal ? -(x - M.thr_mean(l)) / (M.thr_sd * M.thr_sd) : log_beta_prime_scalar_d1(x, M.thr_alpha, M.thr_beta);
}

// Prior odds p(gamma') / p(gamma) on the log scale. Add: s edges present before the move (k excluded). Delete: s
// edges present, k included. Beta-Bernoulli(a, b): p(gamma) = B(a + s, b + K - s) / B(a, b).
double edge_log_odds_add(const EdgeModel& M, arma::uword s) {
    if(!M.beta_bernoulli){
        return M.log_prior_odds;
    }
    return std::log(M.bb_a + s) - std::log(M.bb_b + M.n_int - 1 - s);
}
double edge_log_odds_delete(const EdgeModel& M, arma::uword s) {
    if(!M.beta_bernoulli){
        return -M.log_prior_odds;
    }
    return std::log(M.bb_b + M.n_int - s) - std::log(M.bb_a + s - 1);
}

// Model

EdgeModel edge_make_model(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                          bool core, const arma::vec& eta_star, const arma::mat& Ainv,
                          double thresholds_alpha, double thresholds_beta,
                          double interactions_location, double interactions_scale,
                          const std::string& slab, const std::string& threshold_prior,
                          const arma::vec& threshold_mean, double threshold_sd,
                          const std::string& edge_prior, double inclusion_prob, double beta_alpha, double beta_beta,
                          bool use_likelihood) {
    if(slab != "cauchy" && slab != "normal"){
        Rcpp::stop("slab must be 'cauchy' or 'normal'.");
    }
    if(threshold_prior != "betaprime" && threshold_prior != "normal"){
        Rcpp::stop("threshold_prior must be 'betaprime' or 'normal'.");
    }
    if(edge_prior != "bernoulli" && edge_prior != "beta-bernoulli"){
        Rcpp::stop("edge_prior must be 'bernoulli' or 'beta-bernoulli'.");
    }

    EdgeModel M;
    M.X = X;
    M.f = freq;
    M.ncat = n_categories;
    M.P = X.n_cols;
    M.N = X.n_rows;
    M.H = arma::max(n_categories) - 1;
    M.n_thr = arma::accu(n_categories - 1);
    M.n_int = M.P * (M.P - 1) / 2;
    M.n_pars = M.n_thr + M.n_int;

    // --- Layout: threshold offsets and edge endpoints ---
    M.thr_start.set_size(M.P);
    arma::uword pos = 0;
    for(arma::uword p = 0; p < M.P; p++){
        M.thr_start(p) = pos;
        pos += n_categories(p) - 1;
    }
    M.edge_i.set_size(M.n_int);
    M.edge_j.set_size(M.n_int);
    arma::uword k = 0;
    for(arma::uword j = 0; j < M.P - 1; j++){
        for(arma::uword i = j + 1; i < M.P; i++){
            M.edge_i(k) = i;
            M.edge_j(k) = j;
            k++;
        }
    }

    // --- CoRe map ---
    M.core = core;
    if(core){
        if(eta_star.n_elem != M.n_pars || Ainv.n_rows != M.n_pars || Ainv.n_cols != M.n_pars){
            Rcpp::stop("eta_star and Ainv must have length / dimension n_pars.");
        }
        M.eta_star = eta_star;
        M.Ainv = Ainv;
    }

    // --- Priors ---
    M.thr_alpha = thresholds_alpha;
    M.thr_beta = thresholds_beta;
    M.int_loc = interactions_location;
    M.int_scale = interactions_scale;
    M.slab_normal = (slab == "normal");
    M.thr_normal = (threshold_prior == "normal");
    M.thr_sd = threshold_sd;
    if(M.thr_normal){
        if(threshold_mean.n_elem != M.n_thr){
            Rcpp::stop("threshold_mean must have length n_thresholds.");
        }
        M.thr_mean = threshold_mean;
    }
    M.beta_bernoulli = (edge_prior == "beta-bernoulli");
    M.log_prior_odds = std::log(inclusion_prob) - std::log1p(-inclusion_prob);
    M.bb_a = beta_alpha;
    M.bb_b = beta_beta;
    M.use_lik = use_likelihood;
    return M;
}

// State and evaluations

void edge_eta_to_mu(const EdgeModel& M, const arma::vec& eta, arma::mat& mu) {
    mu.zeros(M.P, M.H);
    for(arma::uword p = 0; p < M.P; p++){
        for(arma::uword h = 1; h < M.ncat(p); h++){
            mu(p, h - 1) = eta(M.thr_start(p) + h - 1);
        }
    }
}

arma::mat edge_eta_to_sigma(const EdgeModel& M, const arma::vec& eta) {
    arma::mat S(M.P, M.P, arma::fill::zeros);
    for(arma::uword k = 0; k < M.n_int; k++){
        S(M.edge_i(k), M.edge_j(k)) = eta(M.n_thr + k);
        S(M.edge_j(k), M.edge_i(k)) = eta(M.n_thr + k);
    }
    return S;
}

arma::vec edge_psi_to_eta(const EdgeModel& M, const arma::vec& psi) {
    return M.core ? arma::vec(M.eta_star + M.Ainv * (psi - M.eta_star)) : psi;
}

// recompute eta, thresholds and rest scores from psi (also removes drift of the incremental updates)
void edge_refresh_state(const EdgeModel& M, EdgeState& S) {
    S.eta = edge_psi_to_eta(M, S.psi);
    edge_eta_to_mu(M, S.eta, S.mu);
    S.R = M.X * edge_eta_to_sigma(M, S.eta);
}

// direction of interaction k (0-based). Pseudo-likelihood: only the conditionals of i and j change. CoRe: a unit
// step in beta_k moves eta along column k of Ainv
EdgeDirection edge_direction(const EdgeModel& M, arma::uword k) {
    EdgeDirection D;
    if(M.core){
        D.deta = M.Ainv.col(M.n_thr + k);
        edge_eta_to_mu(M, D.deta, D.dmu);
        arma::mat S = edge_eta_to_sigma(M, D.deta);
        D.dR = M.X * S;
        for(arma::uword p = 0; p < M.P; p++){
            if(arma::any(D.dmu.row(p) != 0.0) || arma::any(S.col(p) != 0.0)){
                D.nodes.push_back(p);
            }
        }
        return D;
    }
    arma::uword i = M.edge_i(k), j = M.edge_j(k);
    D.deta.zeros(M.n_pars);
    D.deta(M.n_thr + k) = 1.0;
    D.dmu.zeros(M.P, M.H);
    D.dR.zeros(M.N, M.P);
    D.dR.col(i) = M.X.col(j);
    D.dR.col(j) = M.X.col(i);
    D.nodes = {i, j};
    return D;
}

// log pseudo-likelihood at eta + t * deta, summed over the affected nodes (differences are exact), with its first
// and second derivatives in t (the second is <= 0: the pseudo-likelihood is concave in eta). t is the displacement of
// psi_k from its value in the state S, which is not updated during the edge update: the point evaluated has
// psi_k = psi_start + t, so t = -psi_start gives psi_k = 0. psi_start is fixed, so the derivatives in t are the
// derivatives in psi_k (in beta_k for CoRe).
EdgeLineEval edge_line_eval(const EdgeModel& M, const EdgeState& S, const EdgeDirection& D, double t) {
    EdgeLineEval out = {0.0, 0.0, 0.0};
    if(!M.use_lik){
        return out;
    }
    std::vector<double> lp(M.H + 1), dlp(M.H + 1);
    for(arma::uword p : D.nodes){
        arma::uword C = M.ncat(p);
        for(arma::uword n = 0; n < M.N; n++){
            double dr = D.dR(n, p);
            double r = S.R(n, p) + t * dr;
            lp[0] = 0.0;
            dlp[0] = 0.0;
            double mx = 0.0;
            for(arma::uword h = 1; h < C; h++){
                lp[h] = S.mu(p, h - 1) + t * D.dmu(p, h - 1) + static_cast<double>(h) * r;
                dlp[h] = D.dmu(p, h - 1) + static_cast<double>(h) * dr;
                if(lp[h] > mx){
                    mx = lp[h];
                }
            }
            double z = 0.0;
            for(arma::uword h = 0; h < C; h++){
                z += std::exp(lp[h] - mx);
            }
            double logZ = mx + std::log(z);
            double e1 = 0.0, e2 = 0.0;
            for(arma::uword h = 0; h < C; h++){
                double pr = std::exp(lp[h] - logZ);
                e1 += pr * dlp[h];
                e2 += pr * dlp[h] * dlp[h];
            }
            arma::uword x = static_cast<arma::uword>(M.X(n, p));
            out.value += M.f(n) * (lp[x] - logZ);
            out.d1 += M.f(n) * (dlp[x] - e1);
            out.d2 -= M.f(n) * (e2 - e1 * e1);
        }
    }
    return out;
}

// apply a step of size t along D to eta, the thresholds and the rest scores. psi_k itself is updated in
// edge_selection_sampler() right after this call (to psi_start + t, or exactly zero if the edge was deleted)
void edge_apply_move(const EdgeModel& M, EdgeState& S, const EdgeDirection& D, double t) {
    if(t == 0.0){
        return;
    }
    S.eta += t * D.deta;
    for(arma::uword p : D.nodes){
        S.mu.row(p) += t * D.dmu.row(p);
        S.R.col(p) += t * D.dR.col(p);
    }
}

// log pseudo-likelihood and its gradient in eta. The interaction part of the gradient is computed only if
// with_interactions is true (otherwise those entries stay zero)
double edge_full_eval(const EdgeModel& M, const arma::mat& mu, const arma::mat& R, bool with_interactions,
                      arma::vec& grad) {
    grad.zeros(M.n_pars);
    if(!M.use_lik){
        return 0.0;
    }
    double value = 0.0;
    arma::mat Ex(M.N, M.P, arma::fill::zeros);
    std::vector<double> lp(M.H + 1);
    for(arma::uword p = 0; p < M.P; p++){
        arma::uword C = M.ncat(p);
        for(arma::uword n = 0; n < M.N; n++){
            lp[0] = 0.0;
            double mx = 0.0;
            for(arma::uword h = 1; h < C; h++){
                lp[h] = mu(p, h - 1) + static_cast<double>(h) * R(n, p);
                if(lp[h] > mx){
                    mx = lp[h];
                }
            }
            double z = 0.0;
            for(arma::uword h = 0; h < C; h++){
                z += std::exp(lp[h] - mx);
            }
            double logZ = mx + std::log(z);
            arma::uword x = static_cast<arma::uword>(M.X(n, p));
            value += M.f(n) * (lp[x] - logZ);
            for(arma::uword h = 1; h < C; h++){
                double pr = std::exp(lp[h] - logZ);
                grad(M.thr_start(p) + h - 1) += M.f(n) * (static_cast<double>(x == h) - pr);
                Ex(n, p) += static_cast<double>(h) * pr;
            }
        }
    }
    // speed shortcut: the interaction part costs an extra N x P^2 matrix product and is skipped when the caller needs
    // the threshold gradient only. In the edge-selection samplers block_eval() skips it only for the pseudo-likelihood
    // sampler when no edge is present (block = thresholds), which in practice happens only in prior-only runs or with
    // very weak data. CoRe always needs it (A^{-1} maps every block coordinate to all of eta).
    if(with_interactions){
        // d l / d eta_ij = sum_n f_n [x_nj (x_ni - E_ni) + x_ni (x_nj - E_nj)]
        arma::mat res = M.X - Ex;
        res.each_col() %= M.f;
        arma::mat G = M.X.t() * res; // G(a, b) = sum_n x_na f_n res_nb
        for(arma::uword k = 0; k < M.n_int; k++){
            arma::uword i = M.edge_i(k), j = M.edge_j(k);
            grad(M.n_thr + k) = G(j, i) + G(i, j);
        }
    }
    return value;
}

// log target and gradient of the block A (thresholds first, then present interactions) at psi_A + delta
BlockEval block_eval(const EdgeModel& M, const EdgeState& S, const arma::uvec& A, const arma::vec& delta) {
    BlockEval out;
    out.mu = S.mu;
    out.R = S.R;
    arma::vec psi_A = S.psi(A) + delta;
    bool has_int = A.n_elem > M.n_thr;
    if(arma::any(delta != 0.0)){
        arma::vec deta(M.n_pars, arma::fill::zeros);
        if(M.core){
            deta = M.Ainv.cols(A) * delta;
        }
        else{
            deta(A) = delta;
        }
        arma::mat dmu;
        edge_eta_to_mu(M, deta, dmu);
        out.mu += dmu;
        if(M.core || has_int){
            out.R += M.X * edge_eta_to_sigma(M, deta);
        }
    }
    arma::vec grad_eta;
    out.value = edge_full_eval(M, out.mu, out.R, M.core || has_int, grad_eta);
    out.grad = M.core ? arma::vec(M.Ainv.cols(A).t() * grad_eta) : arma::vec(grad_eta(A));
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

// FisherMALA block (Titsias, 2024)

void fishermala_block_init(BlockFisherMALA& F, arma::uword d, double sigma2) {
    F.d = d;
    F.R_n.eye(d, d);
    F.I_d.eye(d, d);
    F.sigma2 = sigma2;
    F.sigma2_R = sigma2;
}

// end of phase 2: freeze the learned preconditioner M = R_n R_n'
void fishermala_block_freeze(BlockFisherMALA& F) {
    F.M_full = F.R_n * F.R_n.t();
    F.frozen = true;
}

// preconditioner of the block A and step size normalized on it, sigma2 / (tr(S_A) / |A|). Before freezing (phases
// 1-2, A = all parameters): S = R_n R_n'. After: the Schur block S_A = M_AA - M_AB M_BB^{-1} M_BA (B = outside A),
// the conditional covariance of the block with the absent edges at exactly zero.
BlockPrecond fishermala_block_precond(BlockFisherMALA& F, const arma::uvec& A) {
    BlockPrecond P;
    if(!F.frozen){
        P.L = F.R_n.rows(A);
        P.S = P.L * P.L.t();
    }
    else{
        P.S = F.M_full.submat(A, A);
        if(A.n_elem < F.d){
            arma::uvec in_A(F.d, arma::fill::zeros);
            in_A(A).ones();
            arma::uvec B = arma::find(in_A == 0);
            arma::mat M_AB = F.M_full.submat(A, B);
            P.S -= M_AB * arma::solve(arma::symmatu(F.M_full.submat(B, B)), M_AB.t(), arma::solve_opts::likely_sympd);
        }
        P.S = arma::symmatu(P.S);
        if(!arma::chol(P.L, P.S, "lower")){
            // numerical fallback: the marginal block
            P.S = arma::symmatu(F.M_full.submat(A, A));
            P.L = arma::chol(P.S, "lower");
        }
    }
    F.sigma2_R = F.sigma2 / (arma::trace(P.S) / static_cast<double>(A.n_elem));
    return P;
}

arma::vec fishermala_block_propose(const BlockFisherMALA& F, const BlockPrecond& P, const arma::vec& x,
                                          const arma::vec& g) {
    arma::vec z(P.L.n_cols);
    for(arma::uword i = 0; i < P.L.n_cols; i++){
        z(i) = R::norm_rand();
    }
    return x + (F.sigma2_R / 2.0) * (P.S * g) + std::sqrt(F.sigma2_R) * (P.L * z);
}

// R_n is adapted only before it is frozen, sigma2 whenever adapting
void fishermala_block_adapt(BlockFisherMALA& F, double log_a, const arma::vec& gx, const arma::vec& gy,
                                   arma::uword s, arma::uword adaptive_stage_n_iter) {
    if(!F.frozen && gy.is_finite()){
        if(s == adaptive_stage_n_iter){
            double gg = arma::dot(gy, gy);
            double r = 1.0 / (1.0 + std::sqrt(F.lambda / (F.lambda + gg)));
            F.R_n = (1.0 / std::sqrt(F.lambda)) * (F.I_d - r * (gy * gy.t()) / (F.lambda + gg));
        }
        else if(s > adaptive_stage_n_iter){
            arma::vec signal = std::sqrt(std::exp(log_a)) * (gy - gx);
            arma::vec phi = F.R_n.t() * signal;
            double pp = arma::dot(phi, phi);
            double r = 1.0 / (1.0 + std::sqrt(1.0 / (1.0 + pp)));
            F.R_n -= r * ((F.R_n * phi) * phi.t()) / (1.0 + pp);
        }
    }
    F.sigma2 *= 1.0 + F.learning_rate * (std::exp(log_a) - F.target);
}

// Log acceptance ratios

// within move on a present edge: psi_current -> psi_proposed, Langevin proposal with means mean_forward (from
// psi_current) and mean_backward (from psi_proposed) and variance var
double edge_within_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& current, const EdgeLineEval& proposed,
                                        double psi_current, double psi_proposed,
                                        double mean_forward, double mean_backward, double var) {
    double log_a = proposed.value - current.value
                   + edge_slab_log(M, psi_proposed)
                   - edge_slab_log(M, psi_current)
                   + edge_log_dnorm(psi_current, mean_backward, var) - edge_log_dnorm(psi_proposed, mean_forward, var);
    if(std::isnan(log_a)){
        log_a = -arma::datum::inf;
    }
    return log_a;
}

// add an absent edge: 0 -> psi_proposed ~ q_k = N(q_mean, q_var). The current state is at_zero (psi_k = 0), with
// n_present edges present
double edge_add_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& at_zero, const EdgeLineEval& proposed,
                                     double psi_proposed, double q_mean, double q_var, arma::uword n_present) {
    return proposed.value - at_zero.value + edge_slab_log(M, psi_proposed)
           - edge_log_dnorm(psi_proposed, q_mean, q_var) + edge_log_odds_add(M, n_present);
}

// delete a present edge: psi_current -> 0, with n_present edges present (edge k included). The reverse (add) move
// would propose psi_current from q_k
double edge_delete_log_acceptance_ratio(const EdgeModel& M, const EdgeLineEval& current, const EdgeLineEval& at_zero,
                                        double psi_current, double q_mean, double q_var, arma::uword n_present) {
    return at_zero.value - current.value + edge_log_dnorm(psi_current, q_mean, q_var)
           - edge_slab_log(M, psi_current) + edge_log_odds_delete(M, n_present);
}

// FisherMALA block x -> y: target difference plus the Langevin proposal ratio, log q(x | y) - log q(y | x),
// in the form of Titsias (2024, Proposition 1), which needs only products with S
double block_log_acceptance_ratio(const BlockFisherMALA& F, const BlockPrecond& P, const BlockEval& current,
                                  const BlockEval& proposed, const arma::vec& x, const arma::vec& y) {
    arma::vec step_x = x - y - (F.sigma2_R / 4.0) * (P.S * proposed.grad);
    arma::vec step_y = y - x - (F.sigma2_R / 4.0) * (P.S * current.grad);
    double log_q_ratio = 0.5 * (arma::dot(step_x, proposed.grad) - arma::dot(step_y, current.grad));
    double log_a = proposed.value - current.value + log_q_ratio;
    if(std::isnan(log_a)){
        log_a = -arma::datum::inf;
    }
    if(log_a > 0.0){
        log_a = 0.0;
    }
    return log_a;
}

// Moves

// within move on present edge k: one-dimensional MALA with a per-edge step size. The edge sits at psi_start + t (t is
// its displacement in the current edge update). The step starts at the inverse local information and, during burn-in,
// is tuned towards acceptance 0.574 by Robbins-Monro, log eps2 += n^(-0.75) (alpha - 0.574) with n the number of
// updates of this edge (Marsman et al., 2025). The decreasing gain corrects the start quickly, even for edges that
// are rarely present, and then settles. The step is fixed after burn-in.
void edge_within_move(const EdgeModel& M, const EdgeState& S, const EdgeDirection& D, arma::uword k,
                      double psi_start, double& t, EdgeLineEval& current,
                      EdgeAdaptation& adapt, EdgeCounters& counts, bool adapting, bool keep_stats) {
    const double target_ar = 0.574, rm_exponent = 0.75;
    double psi_current = psi_start + t;
    double grad_current = current.d1 + edge_slab_d1(M, psi_current);
    if(std::isnan(adapt.log_eps2(k))){
        // first visit: step from the local information
        double info = -current.d2 + edge_slab_info(M);
        adapt.log_eps2(k) = std::log(1.0 / info);
    }
    double var = std::exp(adapt.log_eps2(k));
    double mean_forward = psi_current + 0.5 * var * grad_current;
    double psi_proposed = mean_forward + std::sqrt(var) * R::norm_rand();
    EdgeLineEval proposed = edge_line_eval(M, S, D, psi_proposed - psi_start);
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

// add/delete move of edge k (MoMS). The Newton proposal q_k = N(g0 / I0, 1 / I0) is built at psi_k = 0 with the
// rest fixed, so the add move and its reverse delete move use the same density. Updates gamma, n_present, t and
// current, and accumulates the Rao-Blackwellized inclusion gamma + (1 - 2 gamma) alpha.
void edge_between_move(const EdgeModel& M, EdgeState& S, const EdgeDirection& D, arma::uword k,
                       double psi_start, double& t, EdgeLineEval& current,
                       EdgeCounters& counts, bool keep_stats) {
    EdgeLineEval at_zero = (S.gamma(k) == 0) ? current : edge_line_eval(M, S, D, -psi_start);
    double g0 = at_zero.d1 + edge_slab_d1(M, 0.0);
    // curvature at psi_k = 0: the likelihood part -d2 (>= 0, the log likelihood is concave) plus the slab curvature,
    // floored at the slab information so that I0 > 0. With location 0 the floor is never active (Cauchy: 2 / s^2
    // against 1 / (2 s^2), normal: equal). It acts only for a Cauchy slab centered away from zero
    // (|location| > 0.68 s), whose curvature at 0 is small or negative.
    double I0 = -at_zero.d2 + std::max(-edge_slab_d2(M, 0.0), edge_slab_info(M));
    double q_mean = g0 / I0;
    double q_var = 1.0 / I0;

    if(S.gamma(k) == 0){
        // edge absent, add it: (0, 0) -> (1, psi_proposed). The forward move draws psi_proposed from q_k and the
        // reverse move (delete) is deterministic, so the ratio has -log q_k(psi_proposed)
        double psi_proposed = q_mean + std::sqrt(q_var) * R::norm_rand();
        EdgeLineEval proposed = edge_line_eval(M, S, D, psi_proposed);
        double log_a = edge_add_log_acceptance_ratio(M, at_zero, proposed, psi_proposed, q_mean, q_var, S.n_present);
        double acc = std::isnan(log_a) ? 0.0 : std::exp(std::min(0.0, log_a));
        if(keep_stats){
            counts.edge_counts(k, 0)++;
            counts.rb_sum(k) += acc;
        }
        if(R::unif_rand() < acc){
            S.gamma(k) = 1;
            S.n_present++;
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
                                                        S.n_present);
        double acc = std::isnan(log_a) ? 0.0 : std::exp(std::min(0.0, log_a));
        if(keep_stats){
            counts.edge_counts(k, 2)++;
            counts.rb_sum(k) += 1.0 - acc;
        }
        if(R::unif_rand() < acc){
            S.gamma(k) = 0;
            S.n_present--;
            t = -psi_start;
            current = at_zero;
            if(keep_stats){
                counts.edge_counts(k, 3)++;
            }
        }
    }
}

// FisherMALA block move on the parameters A, gamma fixed
void fishermala_block_move(const EdgeModel& M, EdgeState& S, BlockFisherMALA& F, const arma::uvec& A,
                           arma::uword s, arma::uword adaptive_stage_n_iter,
                           EdgeCounters& counts, bool adapting, bool keep_stats) {
    BlockPrecond P = fishermala_block_precond(F, A);
    BlockEval current = block_eval(M, S, A, arma::vec(A.n_elem, arma::fill::zeros));
    arma::vec x = S.psi(A);
    arma::vec y = fishermala_block_propose(F, P, x, current.grad);
    BlockEval proposed = block_eval(M, S, A, y - x);
    double log_a = block_log_acceptance_ratio(F, P, current, proposed, x, y);
    if(keep_stats){
        counts.block_tried++;
    }
    if(R::unif_rand() < std::exp(log_a)){
        S.psi(A) = y;
        S.mu = proposed.mu;
        S.R = proposed.R;
        S.eta = edge_psi_to_eta(M, S.psi);
        if(keep_stats){
            counts.block_accepted++;
        }
    }
    if(adapting){
        fishermala_block_adapt(F, log_a, current.grad, proposed.grad, s + 1, adaptive_stage_n_iter);
    }
}

// Sampler: phases and iterations

// shuffle the edge order in place (using the Fisher-Yates shuffle with R's random number generator)
void edge_shuffle(arma::uvec& order) {
    for(arma::uword a = order.n_elem - 1; a > 0; a--){
        arma::uword b = static_cast<arma::uword>(std::floor(R::unif_rand() * static_cast<double>(a + 1)));
        if(b > a){
            b = a;
        }
        std::swap(order(a), order(b));
    }
}

Rcpp::List edge_selection_sampler(const EdgeModel& M, const arma::vec& pars,
                                  arma::uword nsim, arma::uword burnin, arma::uword adaptive_stage_n_iter,
                                  arma::uword thin, double sigma2, bool progress) {
    if(pars.n_elem != M.n_pars){
        Rcpp::stop("pars must have length n_pars.");
    }

    // --- Phases: (1) iterations 0 .. adaptive_stage_n_iter - 1, (2) up to burnin / 2, (3) up to burnin, then kept ---
    arma::uword selection_start = burnin / 2;
    arma::uword n_iter = burnin + nsim * thin;
    if(adaptive_stage_n_iter < 1 || thin < 1){
        Rcpp::stop("adaptive_stage_n_iter and thin must be at least 1.");
    }
    if(adaptive_stage_n_iter >= selection_start){
        Rcpp::stop("adaptive_stage_n_iter must be smaller than burnin / 2 (phase 2 would be empty).");
    }

    // --- Initial state: full network ---
    EdgeState S;
    S.psi = pars;
    S.gamma.ones(M.n_int);
    S.n_present = M.n_int;
    edge_refresh_state(M, S);
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

    // --- Progress bar ---
    Progress bar(n_iter, progress);
    arma::uword print_every = std::max<arma::uword>(1, n_iter / 100);

    arma::wall_clock timer;
    timer.tic();

    for(arma::uword s = 0; s < n_iter; s++){
        Rcpp::checkUserInterrupt();
        bool adapting = s < burnin;          // all adaptation stops after burn-in
        bool keep_stats = s >= burnin;
        bool selecting = s >= selection_start;
        if(s == selection_start){
            fishermala_block_freeze(F);
        }
        // once per iteration: recompute eta, the thresholds and the rest scores from psi, so that the rounding errors
        // of the incremental updates (edge_apply_move(), accepted block moves) do not accumulate
        edge_refresh_state(M, S);

        // --- (1) edge pass (phase 3 and kept iterations): within move, then add/delete, edges in a random order ---
        if(selecting){
            edge_shuffle(edge_order);
            for(arma::uword idx = 0; idx < M.n_int; idx++){
                arma::uword k = edge_order(idx);
                EdgeDirection D = edge_direction(M, k);
                double psi_start = S.psi(M.n_thr + k);
                double t = 0.0;
                EdgeLineEval current = edge_line_eval(M, S, D, 0.0);
                if(S.gamma(k) == 1){
                    edge_within_move(M, S, D, k, psi_start, t, current, adapt, counts, adapting, keep_stats);
                }
                edge_between_move(M, S, D, k, psi_start, t, current, counts, keep_stats);
                edge_apply_move(M, S, D, t);
                S.psi(M.n_thr + k) = (S.gamma(k) == 1) ? psi_start + t : 0.0;
            }
        }

        // --- (2) FisherMALA block: all parameters (phases 1-2), or thresholds and present interactions ---
        arma::uvec A = selecting ? arma::uvec(arma::join_cols(thresholds, arma::find(S.gamma == 1) + M.n_thr)) : all_pars;
        fishermala_block_move(M, S, F, A, s, adaptive_stage_n_iter, counts, adapting, keep_stats);
        sigma2_trace(s) = F.sigma2;

        // --- store ---
        if(s >= burnin && (s - burnin) % thin == 0){
            draws.col((s - burnin) / thin) = S.psi;
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

// Helpers exported for testing

// log pseudo-likelihood and its gradient in eta
// @param X unique response patterns [N x P], categories 0, ..., C_p - 1
// @param freq frequency of each pattern
// @param n_categories number of categories per node
// @param eta parameters, thresholds then interactions
// [[Rcpp::export]]
Rcpp::List cpp_paper_edge_logpl(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                                const arma::vec& eta) {
    EdgeModel M = edge_make_model(X, freq, n_categories, false, arma::vec(), arma::mat(), 0.5, 0.5, 0.0, 2.5,
                                  "cauchy", "betaprime", arma::vec(), 1.0, "bernoulli", 0.5, 1.0, 1.0, true);
    arma::mat mu;
    edge_eta_to_mu(M, eta, mu);
    arma::mat R = M.X * edge_eta_to_sigma(M, eta);
    arma::vec grad;
    double value = edge_full_eval(M, mu, R, true, grad);
    return Rcpp::List::create(Rcpp::Named("value") = value, Rcpp::Named("gradient") = grad);
}

// log pseudo-likelihood along interaction k (1-based) of psi at displacement t: value over the affected nodes and
// first two derivatives (core = true: CoRe map with anchor eta_star and inverse rescaling Ainv)
// [[Rcpp::export]]
Rcpp::List cpp_paper_edge_line(const arma::mat& X, const arma::vec& freq, const arma::uvec& n_categories,
                               const arma::vec& psi, bool core, const arma::vec& eta_star, const arma::mat& Ainv,
                               arma::uword k, double t) {
    EdgeModel M = edge_make_model(X, freq, n_categories, core, eta_star, Ainv, 0.5, 0.5, 0.0, 2.5,
                                  "cauchy", "betaprime", arma::vec(), 1.0, "bernoulli", 0.5, 1.0, 1.0, true);
    EdgeState S;
    S.psi = psi;
    edge_refresh_state(M, S);
    EdgeDirection D = edge_direction(M, k - 1);
    EdgeLineEval e = edge_line_eval(M, S, D, t);
    return Rcpp::List::create(Rcpp::Named("value") = e.value, Rcpp::Named("d1") = e.d1, Rcpp::Named("d2") = e.d2,
                              Rcpp::Named("nodes") = D.nodes.size());
}
