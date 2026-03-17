#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

#ifndef UTILS_H
#define UTILS_H

// utility functions for discrete Markov random fields

// Generate a random normal matrix 
arma::mat rnorm_arma(int nrow, int ncol);

// Generate n multivariate normal samples
arma::mat mvnrnd_arma(const arma::vec &mu, const arma::mat &Sigma, int n);

// function to calculate the negative pseudologlikelihood and its derivatives (gradient and hessian) for a discrete MRF model, used by the optimization algorithm
Rcpp::List dmrf_deriv(const arma::vec &pars,
                      const arma::mat &data,
                      const arma::uword &P,
                      const arma::uvec &n_categories,
                      const bool &with_prior = false,
                      const int &ncores = 1,
                      const double &thresholds_alpha = 0.5,
                      const double &thresholds_beta = 0.5,
                      const double &interactions_location = 0.0,
                      const double &interactions_scale = 2.5);

// function to calculate the negative pseudologlikelihood for a discrete MRF model
double npseudologlik(const arma::vec &pars,
                    const arma::mat &data,
                    const arma::uword &P,
                    const arma::uvec &n_categories,
                    const bool &with_prior = false,
                    const int &ncores = 1,
                    const double &thresholds_alpha = 0.5,
                    const double &thresholds_beta = 0.5,
                    const double &interactions_location = 0.0,
                    const double &interactions_scale = 2.5);


#endif