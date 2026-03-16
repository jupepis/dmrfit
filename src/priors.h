#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

#ifndef PRIORS_H
#define PRIORS_H

// Beta function 
double Beta_fun(double x, double y);

// log density of the Cauchy distribution with location l and scale gamma, evaluated at vector x
arma::vec log_dcauchy(arma::vec x, double l = 0.0, double gamma = 2.5);                      

// log density of the Cauchy distribution evaluated at scalar x
double log_pdf_cauchy(double x, double l = 0.0, double gamma = 2.5);

// Cauchy : log first derivative
arma::vec log_dcauchy_first_derivative(arma::vec x, double l = 0.0, double gamma = 2.5);

// Cauchy : log second derivative
arma::vec log_dcauchy_second_derivative(arma::vec x, double l = 0.0, double gamma = 2.5);

// log density of the beta-prime distribution with parameters alpha and beta, evaluated at vector x (note: reparametrization y = exp(x))
arma::vec log_beta_prime(arma::vec x, double alpha = 0.5, double beta = 0.5);

// log density of the beta-prime distribution with parameters alpha and beta, evaluated at scalar x (note: reparametrization y = exp(x))
double log_pdf_beta_prime(double x, double alpha = 0.5, double beta = 0.5);

// log beta-prime : log first derivative
arma::vec log_beta_prime_first_derivative(arma::vec x, double alpha = 0.5, double beta = 0.5);

// log beta-prime : log second derivative
arma::vec log_beta_prime_second_derivative(arma::vec x, double alpha = 0.5, double beta = 0.5);

#endif
