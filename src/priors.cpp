#include <string>
#include <RcppArmadillo.h>
#include "priors.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Beta function and log beta-prime distribution functions
double Beta_fun(double x, double y){
    return std::tgamma(x)*std::tgamma(y)/std::tgamma(x+y);
}

// log_dcauchy : function, first and second derivative

// log density of the Cauchy distribution with location l and scale gamma, evaluated at vector x
arma::vec log_dcauchy(
    arma::vec x, 
    double l, // location parameter
    double gamma // scale parameter
    ){
    x -= l;
    x /= gamma;
    x %= x;
    x += 1.0;
    x *= (gamma*arma::datum::pi);
    return -arma::log(x);
}

// log density of the Cauchy distribution evaluated at scalar x
double log_pdf_cauchy(
    double x, 
    double l, // location parameter
    double gamma // scale parameter
    ){
    x -= l;
    x /= gamma;
    x *= x;
    x += 1.0;
    x *= (gamma*arma::datum::pi);
    return -std::log(x);
}

// Cauchy : log first derivative
arma::vec log_dcauchy_first_derivative(
    arma::vec x, 
    double l, // location parameter
    double gamma // scale parameter
    ){
    x -= l; // (x-l)
    return  ( -2.0 ) * ( x / (gamma * gamma + x % x) );
}

// Cauchy : log second derivative
arma::vec log_dcauchy_second_derivative(
    arma::vec x, 
    double l, // location parameter
    double gamma){ // scale parameter
    x -= l; // (x-l)
    x %= x; // (x-l)^2
    gamma *= gamma; // gamma^2
    arma::vec denom = gamma + x;
    return   2.0 * ( x - gamma ) / ( denom % denom );
}


// log beta-prime : function, first and second derivative

// log density of the beta-prime distribution with parameters alpha and beta, evaluated at scalar x 
double log_pdf_beta_prime(
    double x, // note: reparametrization y = exp(x) 
    double alpha,
    double beta
    ){
    return alpha*x-(alpha+beta)*std::log(1.0+std::exp(x))-std::log(Beta_fun(alpha,beta));
}

// log density of the beta-prime distribution with parameters alpha and beta, evaluated at vector x 
arma::vec log_beta_prime(
    arma::vec x, // note: reparametrization y = exp(x) 
    double alpha,
    double beta
    ){
    return alpha*x-(alpha+beta)*arma::log(1.0+arma::exp(x))-std::log(Beta_fun(alpha,beta));
}

// log beta-prime : log first derivative
arma::vec log_beta_prime_first_derivative(
    arma::vec x, 
    double alpha,
    double beta
    ){
    x = arma::exp(x);
    return alpha - ((alpha+beta) * (x % (1.0/(1.0+x))));
}

// log beta-prime : log second derivative
arma::vec log_beta_prime_second_derivative(
    arma::vec x, 
    double alpha,
    double beta
    ){
    x = arma::exp(x);
    x = x % (1.0/(1.0+x));
    return (-alpha-beta) * (x % (1.0-x));
}

