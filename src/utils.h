#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

#ifndef UTILS_H
#define UTILS_H

// utility functions for discrete Markov random fields

// function to calculate the negative pseudologlikelihood and its derivatives (gradient and hessian) for a discrete MRF model, used by the optimization algorithm
Rcpp::List dmrf_deriv(const arma::vec &pars,
                      const arma::mat &data, // this is already the matrix of sufficient statistics, by row (person) it looks like: {X_1, ..., X_P,2X_1X_2,...,2X_{P-1}X_P}
                      const arma::uword &P,
                      const arma::uvec &n_categories,
                      bool with_prior = false, // whether to include (TRUE) or not (FALSE) prior information (non-informative priors are used)
                      int ncores = 1);

#endif