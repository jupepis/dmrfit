#include <RcppArmadillo.h>
#include "paper-utils.h"

#ifndef PAPER_EXACT_EDGE_H
#define PAPER_EXACT_EDGE_H

// Types and evaluations of the exact-likelihood edge-selection sampler (paper-exact-edge.cpp), also used by the
// DMH-MoMS sampler (paper-dmh-edge.cpp). There the rows of Y are the auxiliary draws instead of all response patterns,
// and the log-weights are relative to the parameter at which the draws were made.

struct ExactData {
    arma::imat Y;     // all response patterns [n_patterns x P], categories 0, ..., C_p - 1 (DMH: auxiliary draws)
    arma::vec S;      // sufficient statistics of the data [n_pars]
    double n;         // number of observations
};

struct ExactState {
    arma::vec psi;      // parameters (thresholds, interactions)
    arma::uvec gamma;   // edge indicators
    arma::uword n_present;
    arma::vec lw;       // pattern log-weights T psi (DMH: T (psi - eta_0))
};

// log target and gradient of the block A at psi_A + delta (gamma fixed). Also returns the new log-weights
struct ExactBlockEval {
    double value;
    arma::vec grad;
    arma::vec lw;
};

arma::vec exact_log_weights(const EdgeModel& M, const ExactData& E, const arma::vec& eta);
ExactBlockEval exact_block_eval(const EdgeModel& M, const ExactData& E, const ExactState& St,
                                const arma::uvec& A, const arma::vec& delta);

#endif
