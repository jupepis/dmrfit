#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Utility functions for optimization algorithm via trust regions
double beta_root_fun(double b, 
                    double r, 
                    double c1, 
                    double c2, 
                    const arma::vec &beta_vec,
                    const arma::vec &gq) {
    double out  = -1.0/r;
    if (b == 0.0) {
        if (c2 > 0.0){
            return out;
        }
        else{
            out += std::sqrt(1.0/c1);
            return out;
        }
    }
    else{
        arma::vec inverse_pb = 1/(beta_vec+b);
        arma::vec lip = gq % inverse_pb;
        out += std::sqrt(1.0/arma::dot(lip,lip));
    return  out;
    }
}

// A function to find the root of `beta_root_fun(b) = 0` using the bisection method
double find_root(double lower, 
                 double upper, 
                 double r, 
                 double c1, 
                 double c2, 
                 const arma::vec &beta_vec, 
                 const arma::vec &gq,
                 double tol = 1e-8) {
  // initial evaluations of `beta_root_fun(b)` at the bounds
  double f_lower = beta_root_fun(lower, r, c1, c2, beta_vec, gq);
  double f_upper = beta_root_fun(upper, r, c1, c2, beta_vec, gq);
  // Make sure that the function values at the bounds have opposite signs
  if (f_lower * f_upper > 0) {
    Rcpp::stop("Function values at the bounds do not have opposite signs.");
  }
  // using bisection method to find the root
  double mid, f_mid;
  while ((upper - lower) / 2.0 > tol) {
    mid = (lower + upper) / 2.0;
    f_mid = beta_root_fun(mid, r, c1, c2, beta_vec, gq);
    
    if (f_mid == 0.0) {
      return mid;  // exact root was found
    } else if (f_lower * f_mid < 0) {
      upper = mid;
      f_upper = f_mid;
    } else {
      lower = mid;
      f_lower = f_mid;
    }
  }
  // return best approximation of the root
  return (lower + upper) / 2.0;
}

// optimization routine for discrete MRFs via trust region algorithm

// [[Rcpp::export]]
Rcpp::List optimize(
    const arma::mat &data,
    const arma::vec &parinit,
    const arma::uvec &n_categories,
    const arma::uword &P,
    const double &f_term,
    const double &m_term,
    const arma::uword &n_iter_max = 100,
    const double &rinit = 1.0,
    const double &rmax =  10.0,
    const bool &with_prior = false,
    const double &epsilon = 1e-06,
    const int &ncores = 1,
    const double &thresholds_alpha = 0.5,
    const double &thresholds_beta = 0.5,
    const double &interactions_location = 0.0,
    const double &interactions_scale = 2.5){

    double r = rinit;
    arma::vec pars = parinit;
    arma::uword n_pars = pars.n_elem;
    Rcpp::List deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha,thresholds_beta, interactions_location,interactions_scale);

    double rho = 0.0;
    double f = 0.0;
    double min_eigval = 0.0;
    bool accept = true;
    bool is_terminate = false;
    bool is_newton = false;
    arma::uword i;

    arma::vec g(n_pars,arma::fill::zeros);
    arma::vec p(n_pars,arma::fill::zeros);
    arma::vec gq(n_pars,arma::fill::zeros);
    arma::mat B(n_pars,n_pars,arma::fill::zeros);

    arma::vec eigval;
    arma::mat eigvec;

    for(i = 0; i < n_iter_max; i++){
            
        if(accept){
            f = deriv["value"];
            g = Rcpp::as<arma::vec>(deriv["gradient"]);
            B = Rcpp::as<arma::mat>(deriv["hessian"]);
            arma::eig_sym(eigval, eigvec, B); 
            gq = eigvec.t() * g;
        }

        // solve trust-region subproblem
        is_newton = false;
        min_eigval = eigval.min();

        // if all eigenvalues are positive... we try to find a Newton solution
        if (min_eigval > 0.0) {
            // Calculate Newton solution 
            p = -(eigvec * (gq / eigval));
            if(arma::norm(p) <= r) { // return solution
                is_newton = true;
            }
        }

        // .. otherwise we try a non-Newton solution
        if(!is_newton){
            arma::vec beta_vec = eigval - min_eigval;
            arma::uvec beta_eq_0 = arma::find(beta_vec == 0.0);
            arma::uvec beta_neq_0 = arma::find(beta_vec != 0.0);


            // Handle the easy and hard cases
            arma::vec beta_inverse = arma::ones(beta_vec.n_elem);
            for (arma::uword j = 0; j < beta_vec.n_elem; j++) {
                if (beta_vec(j) > 1e-8) {
                    beta_inverse(j) = 1.0 / beta_vec(j);
                } 
                else {
                    beta_inverse(j) = 0.0; // avoiding division by zero for small beta
                }
            }

            // Compute constants
            arma::vec gq_beta_neq_0 = gq(beta_neq_0) % beta_inverse(beta_neq_0);
            double c1 = arma::dot(gq_beta_neq_0,gq_beta_neq_0);
            double c2 = arma::dot(gq(beta_eq_0),gq(beta_eq_0));
            double c3 = arma::dot(gq, gq);

            if((c2 > 0.0) || (c1 > (r*r))){ // easy-easy and easy-hard cases
                // handling easy cases
                double beta_down = std::sqrt(c2)/r;
                double beta_up =  std::sqrt(c3)/r;
                double root_fred = 0.0;

                if(beta_root_fun(beta_up,r,c1,c2,beta_vec,gq) <= 0.0){
                    root_fred = beta_up;
                } 
                else if(beta_root_fun(beta_down,r,c1,c2,beta_vec,gq) >= 0.0){
                    root_fred = beta_down;
                } 
                else{
                    root_fred = find_root(beta_down,beta_up,r,c1,c2,beta_vec,gq);
                }
                arma::vec w = gq / (beta_vec + root_fred);
                p = - (eigvec * w);
            }
            else{ // very hard cases
                arma::vec w = gq / beta_vec;
                w(beta_eq_0) *= 0.0;
                p = - (eigvec * w);
                double u = std::sqrt(r*r - arma::dot(p,p));
                if(u > 0.0){
                    arma::mat V = eigvec.cols(beta_eq_0);
                    p += (u * V.col(0));
                }
            }
        }

        double preddiff = arma::dot(p ,(g + (B * p)/2.0));
        arma::vec pars_try = pars + p;


        deriv = dmrf_deriv(pars_try,data,P,n_categories,with_prior,ncores,thresholds_alpha,thresholds_beta, interactions_location,interactions_scale);
        double f_try = deriv["value"];
        rho = (f_try-f)/preddiff;

        if(f_try < arma::datum::inf){
            is_terminate = (std::abs(f_try-f) < f_term) || (abs(preddiff) < m_term);
        } 
        else{
            is_terminate = false;
            rho = -arma::datum::inf;
        }

        if(is_terminate){
            if(f_try <f){
                accept = true;
                pars = pars_try;
            }
        } 
        else{
            if(rho < 0.25){
                accept = false;
                r /= 4.0;
            }
            else{
                accept = true;
                pars = pars_try;
                if(rho > 0.75 && (!is_newton)){
                    r = std::min(2.0*r,rmax);
                }
            }
        }

        if(is_terminate){
            break;
        }

    }

    deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);
    return Rcpp::List::create(Rcpp::Named("argument") = pars, Rcpp::Named("utils") = deriv);
}


// optimization routine for discrete MRFs with fixed structure via trust region algorithm

// [[Rcpp::export]]
Rcpp::List optimize_with_structure(
    const arma::mat &data,
    const arma::vec &parinit,
    const arma::vec &structure,
    const arma::uvec &n_categories,
    const arma::uword &P,
    const double &f_term,
    const double &m_term,
    const arma::uword &n_iter_max = 100,
    const double &rinit = 1.0,
    const double &rmax =  10.0,
    const bool &with_prior = false,
    const double &epsilon = 1e-06,
    const int &ncores = 1,
    const double &thresholds_alpha = 0.5,
    const double &thresholds_beta = 0.5,
    const double &interactions_location = 0.0,
    const double &interactions_scale = 2.5){

    double r = rinit;
    arma::vec pars = parinit % structure;
    arma::uword n_pars = pars.n_elem;
    Rcpp::List deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);

    double rho = 0.0;
    double f = 0.0;
    double min_eigval = 0.0;
    bool accept = true;
    bool is_terminate = false;
    bool is_newton = false;
    arma::uword i;

    arma::vec g(n_pars,arma::fill::zeros);
    arma::vec p(n_pars,arma::fill::zeros);
    arma::vec gq(n_pars,arma::fill::zeros);
    arma::mat B(n_pars,n_pars,arma::fill::zeros);

    arma::vec eigval;
    arma::mat eigvec;

    for(i = 0; i < n_iter_max; i++){
            
        if(accept){
            f = deriv["value"];
            g = Rcpp::as<arma::vec>(deriv["gradient"]) % structure;
            B = Rcpp::as<arma::mat>(deriv["hessian"]) % (structure * structure.t());
            arma::eig_sym(eigval, eigvec, B); 
            gq = eigvec.t() * g;
        }

        // solve trust-region subproblem
        is_newton = false;
        min_eigval = eigval.min();

        // if all eigenvalues are positive... we try to find a Newton solution
        if (min_eigval > 0.0) {
            // Calculate Newton solution 
            p = -(eigvec * (gq / eigval));
            if(arma::norm(p) <= r) { // return solution
                is_newton = true;
            }
        }

        // .. otherwise we try a non-Newton solution
        if(!is_newton){
            arma::vec beta_vec = eigval - min_eigval;
            arma::uvec beta_eq_0 = arma::find(beta_vec == 0.0);
            arma::uvec beta_neq_0 = arma::find(beta_vec != 0.0);


            // Handle the easy and hard cases
            arma::vec beta_inverse = arma::ones(beta_vec.n_elem);
            for (arma::uword j = 0; j < beta_vec.n_elem; j++) {
                if (beta_vec(j) > 1e-8) {
                    beta_inverse(j) = 1.0 / beta_vec(j);
                } 
                else {
                    beta_inverse(j) = 0.0; // avoiding division by zero for small beta
                }
            }

            // Compute constants
            arma::vec gq_beta_neq_0 = gq(beta_neq_0) % beta_inverse(beta_neq_0);
            double c1 = arma::dot(gq_beta_neq_0,gq_beta_neq_0);
            double c2 = arma::dot(gq(beta_eq_0),gq(beta_eq_0));
            double c3 = arma::dot(gq, gq);

            if((c2 > 0.0) || (c1 > (r*r))){ // easy-easy and easy-hard cases
                // handling easy cases
                double beta_down = std::sqrt(c2)/r;
                double beta_up =  std::sqrt(c3)/r;
                double root_fred = 0.0;

                if(beta_root_fun(beta_up,r,c1,c2,beta_vec,gq) <= 0.0){
                    root_fred = beta_up;
                } 
                else if(beta_root_fun(beta_down,r,c1,c2,beta_vec,gq) >= 0.0){
                    root_fred = beta_down;
                } 
                else{
                    root_fred = find_root(beta_down,beta_up,r,c1,c2,beta_vec,gq);
                }
                arma::vec w = gq / (beta_vec + root_fred);
                p = - (eigvec * w);
            }
            else{ // very hard cases
                arma::vec w = gq / beta_vec;
                w(beta_eq_0) *= 0.0;
                p = - (eigvec * w);
                double u = std::sqrt(r*r - arma::dot(p,p));
                if(u > 0.0){
                    arma::mat V = eigvec.cols(beta_eq_0);
                    p += (u * V.col(0));
                }
            }
        }

        double preddiff = arma::dot(p ,(g + (B * p)/2.0));
        arma::vec pars_try = pars + p;
        pars_try %= structure;


        deriv = dmrf_deriv(pars_try,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);
        double f_try = deriv["value"];
        rho = (f_try-f)/preddiff;

        if(f_try < arma::datum::inf){
            is_terminate = (std::abs(f_try-f) < f_term) || (abs(preddiff) < m_term);
        } 
        else{
            is_terminate = false;
            rho = -arma::datum::inf;
        }

        if(is_terminate){
            if(f_try <f){
                accept = true;
                pars = pars_try;
            }
        } 
        else{
            if(rho < 0.25){
                accept = false;
                r /= 4.0;
            }
            else{
                accept = true;
                pars = pars_try;
                if(rho > 0.75 && (!is_newton)){
                    r = std::min(2.0*r,rmax);
                }
            }
        }

        if(is_terminate){
            break;
        }

    }

    deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);
    return Rcpp::List::create(Rcpp::Named("argument") = pars, Rcpp::Named("utils") = deriv);
}



// trust-region algorithm profile negative loglikelihood via parconstr

// [[Rcpp::export]]
Rcpp::List optimize_profile(
    const arma::mat &data,
    const arma::vec &parinit,
    const arma::uvec &which_parconstr,
    const arma::vec &parconstr,
    const arma::uvec &n_categories,
    const arma::uword &P,
    const double &f_term,
    const double &m_term,
    const arma::uword &n_iter_max = 100,
    const double &rinit = 1.0,
    const double &rmax =  10.0,
    const bool &with_prior = false,
    const double &epsilon = 1e-06,
    const int &ncores = 1,
    const double &thresholds_alpha = 0.5,
    const double &thresholds_beta = 0.5,
    const double &interactions_location = 0.0,
    const double &interactions_scale = 2.5){

    double r = rinit;
    arma::vec pars = parinit;
    arma::uword n_pars = pars.n_elem;
    // impose constrained value for certain parameters
    pars(which_parconstr) = parconstr;
    Rcpp::List deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);

    double rho = 0.0;
    double f = 0.0;
    double min_eigval = 0.0;
    bool accept = true;
    bool is_terminate = false;
    bool is_newton = false;
    arma::uword i;

    arma::vec g(n_pars,arma::fill::zeros);
    arma::vec p(n_pars,arma::fill::zeros);
    arma::vec gq(n_pars,arma::fill::zeros);
    arma::mat B(n_pars,n_pars,arma::fill::zeros);

    arma::vec eigval;
    arma::mat eigvec;

    for(i = 0; i < n_iter_max; i++){
            
        if(accept){
            f = deriv["value"];
            g = Rcpp::as<arma::vec>(deriv["gradient"]);
            B = Rcpp::as<arma::mat>(deriv["hessian"]);
            arma::eig_sym(eigval, eigvec, B); 
            gq = eigvec.t() * g;
        }

        // solve trust-region subproblem
        is_newton = false;
        min_eigval = eigval.min();

        // if all eigenvalues are positive... we try to find a Newton solution
        if (min_eigval > 0.0) {
            // Calculate Newton solution 
            p = -(eigvec * (gq / eigval));
            if(arma::norm(p) <= r) { // return solution
                is_newton = true;
            }
        }

        // .. otherwise we try a non-Newton solution
        if(!is_newton){
            arma::vec beta_vec = eigval - min_eigval;
            arma::uvec beta_eq_0 = arma::find(beta_vec == 0.0);
            arma::uvec beta_neq_0 = arma::find(beta_vec != 0.0);


            // Handle the easy and hard cases
            arma::vec beta_inverse = arma::ones(beta_vec.n_elem);
            for (arma::uword j = 0; j < beta_vec.n_elem; j++) {
                if (beta_vec(j) > 1e-8) {
                    beta_inverse(j) = 1.0 / beta_vec(j);
                } 
                else {
                    beta_inverse(j) = 0.0; // avoiding division by zero for small beta
                }
            }

            // Compute constants
            arma::vec gq_beta_neq_0 = gq(beta_neq_0) % beta_inverse(beta_neq_0);
            double c1 = arma::dot(gq_beta_neq_0,gq_beta_neq_0);
            double c2 = arma::dot(gq(beta_eq_0),gq(beta_eq_0));
            double c3 = arma::dot(gq, gq);

            if((c2 > 0.0) || (c1 > (r*r))){ // easy-easy and easy-hard cases
                // handling easy cases
                double beta_down = std::sqrt(c2)/r;
                double beta_up =  std::sqrt(c3)/r;
                double root_fred = 0.0;

                if(beta_root_fun(beta_up,r,c1,c2,beta_vec,gq) <= 0.0){
                    root_fred = beta_up;
                } 
                else if(beta_root_fun(beta_down,r,c1,c2,beta_vec,gq) >= 0.0){
                    root_fred = beta_down;
                } 
                else{
                    root_fred = find_root(beta_down,beta_up,r,c1,c2,beta_vec,gq);
                }
                arma::vec w = gq / (beta_vec + root_fred);
                p = - (eigvec * w);
            }
            else{ // very hard cases
                arma::vec w = gq / beta_vec;
                w(beta_eq_0) *= 0.0;
                p = - (eigvec * w);
                double u = std::sqrt(r*r - arma::dot(p,p));
                if(u > 0.0){
                    arma::mat V = eigvec.cols(beta_eq_0);
                    p += (u * V.col(0));
                }
            }
        }

        double preddiff = arma::dot(p ,(g + (B * p)/2.0));
        arma::vec pars_try = pars + p;
        // impose constrained value for certain parameters
        pars_try(which_parconstr) = parconstr;


        deriv = dmrf_deriv(pars_try,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);
        double f_try = deriv["value"];
        rho = (f_try-f)/preddiff;

        if(f_try < arma::datum::inf){
            is_terminate = (std::abs(f_try-f) < f_term) || (abs(preddiff) < m_term);
        } 
        else{
            is_terminate = false;
            rho = -arma::datum::inf;
        }

        if(is_terminate){
            if(f_try <f){
                accept = true;
                pars = pars_try;
            }
        } 
        else{
            if(rho < 0.25){
                accept = false;
                r /= 4.0;
            }
            else{
                accept = true;
                pars = pars_try;
                if(rho > 0.75 && (!is_newton)){
                    r = std::min(2.0*r,rmax);
                }
            }
        }

        if(is_terminate){
            break;
        }

    }

    // impose constrained value for certain parameters
    pars(which_parconstr) = parconstr;
    deriv = dmrf_deriv(pars,data,P,n_categories,with_prior,ncores,thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);
    return Rcpp::List::create(Rcpp::Named("argument") = pars, Rcpp::Named("utils") = deriv);
}
