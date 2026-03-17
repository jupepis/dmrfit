#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Generate a random normal matrix 
arma::mat rnorm_arma(int nrow, int ncol) {
  Rcpp::NumericVector v = Rcpp::rnorm(nrow * ncol);      
  return arma::mat(v.begin(), nrow, ncol);  // fill arma::mat column-wise
}

// Generate n multivariate normal samples
// [[Rcpp::export]]
arma::mat mvnrnd_arma(const arma::vec &mu, const arma::mat &Sigma, int n) {
  int p = mu.n_elem;
  arma::mat Z = rnorm_arma(p, n);   // each column is a sample
  arma::mat C = arma::chol(Sigma, "lower");
  arma::mat X = C * Z; // from standard to correlated normal
  X.each_col() += mu;
  return X;
}

// Utils Discrete MRF Pseudolikelihood 
//
// pseudolikelihood of a discrete mrf function value and derivatives value at specific parameters (used by the optimization algortihm)
// the function returns a list of named objects
// "loglik" is the negative pseudologlikelihood calculated at the parameters supplied via the input pars
// "gradient" is the negative gradient value at pars
// "hessian" is the negative hessian value (its inverse returns the matrix of variances and covariances of the model parameters)
// "HW" is the Huber-White sandwich variance estimator (it is already a matrix of variances and covariances)
// note: loglik, gradient and hessian are used by the trust algorithm to find the MPLEs
Rcpp::List dmrf_deriv(
    const arma::vec &pars,
    const arma::mat &data, // this is already the matrix of sufficient statistics, by row (person) it looks like: {X_1, ..., X_P,2X_1X_2,...,2X_{P-1}X_P}
    const arma::uword &P,
    const arma::uvec &n_categories,
    const bool &with_prior, // whether to include (TRUE) or not (FALSE) prior information (non-informative priors are used)
    const int &ncores,
    const double &thresholds_alpha,
    const double &thresholds_beta,
    const double &interactions_location,
    const double &interactions_scale)
{
    arma::uword n,p,i,j,h;
    arma::uword n_pars = pars.n_elem; // this must be equal to P+P*(P-1)/2 or to pars.n_elem
    arma::uword N = data.n_rows;
    arma::uword n_thresholds = arma::accu(n_categories-1);
    arma::vec thresholds = pars(arma::span(0,n_thresholds-1));  // vector of thresholds parameters
    arma::vec interactions_vec = pars(arma::span(n_thresholds,n_pars-1));// vector of P*(P-1)/2 interaction parameters

    // building the symmetric matrix of interaction parameters from its lower triangular
    arma::mat lower_matrix_interactions(P,P,arma::fill::zeros);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(lower_matrix_interactions), -1); // element indices of the lower triangular excluding the diagonal elements 
    lower_matrix_interactions(lower_indices) = interactions_vec;
    arma::mat interactions = arma::symmatl(lower_matrix_interactions);

    // check on n_pars
    if((n_thresholds+P*(P-1)/2) != n_pars){ // check that data has the same columns as the number of parameters in the model
        Rcpp::Rcout << "The length of pars must be equal to " << (n_thresholds+P*(P-1)/2) << "\n";
        // stop algorithm (this check can be handled at R-level)
    }

    // utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects 
    arma::umat matrix_indices_sigma(2,P*(P-1)/2);
    arma::uword l = 0;
    for(j = 0; j < (P-1); j++){ // column j
        for(i = (j+1); i < P; i++){ // row i
            matrix_indices_sigma(0,l) = i;
            matrix_indices_sigma(1,l) = j;
            l++;
        }
    }

    // utility vectors indicating which stats and which category (of length n_thresholds - used in the hessian computation)
    arma::uvec which_stats(n_thresholds,arma::fill::zeros);
    arma::vec category_stats(n_thresholds,arma::fill::zeros);
    l = 0; // reset index operator l
    for(p = 0; p < P; p++){
        for(h = 1; h < n_categories(p); h++){
            which_stats(l) = p;
            category_stats(l) = static_cast<double>(h);
            l++;
        }
    }


    // creating empty objects where to save loglik, gradient and hessian computed per each person (statistical unit) , this is useful for the parallelization step
    arma::vec loglik_vec(N,arma::fill::zeros);
    arma::mat gradient_mat(n_pars,N,arma::fill::zeros);
    arma::cube hessian_cube(n_pars,n_pars,N,arma::fill::zeros);
    arma::cube square_score_cube(n_pars,n_pars,N,arma::fill::zeros);

    // # pragma omp parallel for if(ncores>1) private() shared() ... parallelize here (??) the for loop over people
    // loop over people  
    #ifdef _OPENMP
    omp_set_dynamic(0);         
    omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for if(ncores>1) private(n,p,h,i,j) shared(N,P,n_thresholds,n_pars,data,n_categories,thresholds,interactions_vec,interactions,matrix_indices_sigma,category_stats,which_stats,loglik_vec,gradient_mat,hessian_cube)
    #endif
    for(n = 0; n < N; n++){
        // select n-th person statistics
        arma::vec stats_n = data.row(n).t(); // stats for n-th person {X1,X2,...,2XiXj}
        // processing observed category per each X
        arma::uvec stats_X_n = arma::conv_to<arma::uvec>::from(stats_n(arma::span(0,P-1)));
        arma::vec pars_n(P+P*(P-1)/2,arma::fill::zeros);
        for(p = 0; p < P; p++){
            arma::uword which_threshold = stats_X_n(p);
            if(which_threshold > 0){
                arma::uword index_threshold_Xp = (which_threshold-1);
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                pars_n(p) = thresholds(index_threshold_Xp); // which_threshold-1 because in the thresholds vector we omit the baseline
            }
        }
        // filling in the interaction parameters
        pars_n(arma::span(P,pars_n.n_elem-1)) = interactions_vec;
        // calculating the numerator of the pseudo-loglikelihood for participant n
        double loglik_n = arma::accu(pars_n(arma::span(0,P-1))); // sum_p{sum_h{threshold_h*I(x == h)}}
        loglik_n += arma::accu((pars_n(arma::span(P,pars_n.n_elem-1)).t() * stats_n(arma::span(P,pars_n.n_elem-1)))); // sum_p{sum{2x_p*x_j*sigma_pj}

        //for each p we have to compute the support (normalizing constant, denominator) for each item, that is ln[1+sum_h{exp(mu_h+h*sum_{j!=p}{x_jsigma_pj})}]
        double support_n = 0.0;
        arma::vec probs_n(n_thresholds,arma::fill::zeros);
        arma::vec var_n(n_thresholds,arma::fill::zeros);
        for(p = 0; p < P; p++){
            double denom_p = 1.0;
            arma::vec stats_excl_p = stats_n(arma::span(0,P-1));
            stats_excl_p(p) = 0.0; // this could be avoided because the interaction matrix has 0.0 in the diagonal
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                double success_event = arma::accu(thresholds(index_threshold_Xp) + static_cast<double>(h)*(stats_excl_p.t() * interactions.col(p)));                 
                // calculating denom_p, and probs_n numerator
                denom_p += std::exp(success_event);
                probs_n(index_threshold_Xp) = std::exp(success_event); // success probability for p-th variable
            }
            support_n += std::log(denom_p);
        
            // calculate pr(Xp = h) and variances 
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                probs_n(index_threshold_Xp) /= denom_p;
                var_n(index_threshold_Xp) = probs_n(index_threshold_Xp)*(1.0-probs_n(index_threshold_Xp)); // variance (for a bernoulli) p-th variable for level h
            }
        }
        
        // updating loglikelihood
        loglik_n -= support_n;
        loglik_vec(n) -= loglik_n; // (-=) because negative loglikelihood       

        // calculate (negative) gradient
        arma::vec gradient_n(n_pars,arma::fill::zeros);

        // gradient for thresholds
        for(p = 0; p < P; p++){
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                gradient_n(index_threshold_Xp) += static_cast<double>((stats_n(p) == h)) - probs_n(index_threshold_Xp);
            }
        }

        // gradient for interactions
        gradient_n(arma::span(n_thresholds,n_pars-1)) += stats_n(arma::span(P,P*(P-1)/2+P-1)); // summing first part of the gradient, that is the sufficient statistic 2XiXj

        arma::uword index_ij = n_thresholds;
        for(j = 0; j < (P-1); j++){
            // calculate expected Xj
            double expected_Xj = 0.0;
            for(h = 1; h < n_categories(j); h++){
                arma::uword index_threshold_Xj = h-1;
                if(j>0){
                    index_threshold_Xj += arma::accu(n_categories(arma::span(0,j-1))-1);
                }
                expected_Xj += static_cast<double>(h)*probs_n(index_threshold_Xj);
            }
            for(i = (j+1); i < P; i++){
                // calculate expected Xi
                double expected_Xi = 0.0;
                for(h = 1; h < n_categories(i); h++){
                    arma::uword index_threshold_Xi = h-1;
                    if(i>0){
                        index_threshold_Xi += arma::accu(n_categories(arma::span(0,i-1))-1);
                    }
                    expected_Xi += static_cast<double>(h)*probs_n(index_threshold_Xi);
                }
                gradient_n(index_ij) -= (stats_n(j)*expected_Xi + stats_n(i)*expected_Xj); // subtracting the second part that is -Xj*\sum_h{h*P(Xi=h)}-Xi*\sum_h{h(P(Xj=h)}
                index_ij++;
            }
        }

        gradient_mat.col(n) -= gradient_n; // (-=) because negative gradient

        // covariance matrix of score n-th person
        square_score_cube.slice(n) = gradient_mat.col(n) * gradient_mat.col(n).t();
 
        // calculate (negative) Hessian
        arma::mat hessian_n(n_pars,n_pars,arma::fill::zeros);

        // hessian
        for(j = 0; j < n_pars; j++){ // column j
            for(i = j; i < n_pars; i++){ // row i
                if((j < n_thresholds) && (i < n_thresholds)){ // hessian for (thresholds) - only when i == j (diagonal elements), off diagonal elements remain 0.0
                    if(i == j){ // for (threshold_h,threshold_h) of the same X
                        hessian_n(i,j) -= var_n(j);
                        //if(n==0){
                        //    Rcpp::Rcout << "hessian value for mu(" << which_stats(j) << "," << category_stats(j) << ") = " << var_n(j) << "\n";
                        //}
                    }
                    else if(which_stats(i) == which_stats(j)){ // for(threshold_h,threshold_k) of the same node X
                        hessian_n(i,j) += probs_n(i)*probs_n(j);
                        hessian_n(j,i) = hessian_n(i,j);
                    }
                }
                else if((j < n_thresholds) && (i >= n_thresholds) && arma::any(matrix_indices_sigma.col(i-n_thresholds) == which_stats(j))){ // hessian for (thresholds,interactions) - only for (mu_k,sigma_{kl}) or (mu_l,sigma_{kl}), derivatives for (mu_k,sigma_{rl}) remain 0.0
                    arma::uword s = (!(matrix_indices_sigma(1,i-n_thresholds) == which_stats(j)))*1; // selecting position index of the element different from which_stats(j)
                    arma::uword which_index = matrix_indices_sigma(s,i-n_thresholds);
                    // calculate variance of variable X_which_stats(j)
                    double expected_X = 0.0;
                    for(h = 1; h < n_categories(which_stats(j)); h++){
                        arma::uword index_threshold_h = h-1;
                        if(which_stats(j)>0){
                            index_threshold_h += arma::accu(n_categories(arma::span(0,which_stats(j)-1))-1);
                        }
                        expected_X += static_cast<double>(h)*probs_n(index_threshold_h);
                    }
                    // updating hessian
                    hessian_n(i,j) -= stats_n(which_index)*probs_n(j)*(category_stats(j)-expected_X); 
                    hessian_n(j,i) = hessian_n(i,j);
                } 
                else if((i >= n_thresholds) && (j >= n_thresholds)){ // hessian for (interactions) 
                    arma::uvec indices_i = matrix_indices_sigma.col(i-n_thresholds);
                    arma::uvec indices_j = matrix_indices_sigma.col(j-n_thresholds);
                    if(i == j){ // derivative for (sigma_{ij},sigma_{ij}) - here using either indices_j or indices_i is the same as they refer to the same sigma_{kl}
                        // calculate variance of variable X indices_j(0)
                        double expected_X_j0 = 0.0;
                        double expected_X_square_j0 = 0.0;
                        for(h = 1; h < n_categories(indices_j(0)); h++){
                            arma::uword index_threshold_h = h-1;
                            if(indices_j(0)>0){
                                index_threshold_h += arma::accu(n_categories(arma::span(0,indices_j(0)-1))-1);
                            }
                            expected_X_j0 += static_cast<double>(h)*probs_n(index_threshold_h);
                            expected_X_square_j0 += std::pow(static_cast<double>(h),2)*probs_n(index_threshold_h);
                        }
                        double var_j0 = expected_X_square_j0 - std::pow(expected_X_j0,2);

                        // calculate variance of variable X indices_j(1)
                        double expected_X_j1 = 0.0;
                        double expected_X_square_j1 = 0.0;
                        for(h = 1; h < n_categories(indices_j(1)); h++){
                            arma::uword index_threshold_h = h-1;
                            if(indices_j(1)>0){
                                index_threshold_h += arma::accu(n_categories(arma::span(0,indices_j(1)-1))-1);
                            }
                            expected_X_j1 += static_cast<double>(h)*probs_n(index_threshold_h);
                            expected_X_square_j1 += std::pow(static_cast<double>(h),2)*probs_n(index_threshold_h);
                        }
                        double var_j1 = expected_X_square_j1 - std::pow(expected_X_j1,2);

                        hessian_n(i,j) -= (std::pow(stats_n(indices_j(1)),2)*var_j0+std::pow(stats_n(indices_j(0)),2)*var_j1);
                    }
                    else if(arma::any(indices_i == indices_j(0)) || arma::any(indices_i == indices_j(1))){ // derivative for (sigma_{kl},sigma_{lr}), at this stage only one of the two conditions on indices_i can be true
                        arma::uvec m = arma::join_cols(arma::find(indices_i == indices_j(0)),arma::find(indices_i == indices_j(1)));
                        arma::uword index_in_common =  indices_i(m(0));
                        arma::uvec which_indices_j = arma::find(indices_j != index_in_common);
                        arma::uword index_j = indices_j(which_indices_j(0));
                        arma::uvec which_indices_i = arma::find(indices_i != index_in_common);
                        arma::uword index_i = indices_i(which_indices_i(0));

                        // calculate variance of variable or the X in common (if sigma_ij and sigma_il, is the variance of Xi)
                        double expected_X = 0.0;
                        double expected_X_square = 0.0;
                        for(h = 1; h < n_categories(index_in_common); h++){
                            arma::uword index_threshold_h = h-1;
                            if(index_in_common>0){
                                index_threshold_h += arma::accu(n_categories(arma::span(0,index_in_common-1))-1);
                            }
                            expected_X += static_cast<double>(h)*probs_n(index_threshold_h);
                            expected_X_square += std::pow(static_cast<double>(h),2)*probs_n(index_threshold_h);
                        }
                        double var_X = expected_X_square - std::pow(expected_X,2);

                        hessian_n(i,j) -= stats_n(index_i)*stats_n(index_j)*var_X;
                        hessian_n(j,i) = hessian_n(i,j);
                    }
                }
            }
        }
        hessian_cube.slice(n) -= hessian_n; // (-=) because negative hessian
        
    }
    
    // sum over people
    double loglik = arma::accu(loglik_vec); 
    arma::vec gradient = arma::sum(gradient_mat,1); 
    arma::mat hessian = arma::sum(hessian_cube,2); 

    // calculate Huber-White sandwich variance estimator

    arma::mat inverse_negative_hessian = arma::inv_sympd(hessian);
    arma::mat square_score = arma::sum(square_score_cube,2);  //gradient * gradient.t(); 
    arma::mat HW = inverse_negative_hessian * square_score;
    HW *= inverse_negative_hessian;

    // Adding prior information on loglik, gradient and hessian 
    if(with_prior){
        // prior on pseudologlikelihood
        loglik -= (arma::accu(log_beta_prime(thresholds, thresholds_alpha, thresholds_beta)) + arma::accu(log_dcauchy(interactions_vec, interactions_location, interactions_scale)));
        // prior on gradient
        arma::vec log_prior_gradient_thresholds = log_beta_prime_first_derivative(thresholds, thresholds_alpha, thresholds_beta);
        arma::vec log_prior_gradient_interactions = log_dcauchy_first_derivative(interactions_vec, interactions_location, interactions_scale);
        arma::vec gradient_prior = arma::join_cols(log_prior_gradient_thresholds,log_prior_gradient_interactions);
        gradient -= gradient_prior;
        // prior on the hessian
        arma::vec log_prior_hessian_thresholds = log_beta_prime_second_derivative(thresholds, thresholds_alpha, thresholds_beta);
        arma::vec log_prior_hessian_interactions = log_dcauchy_second_derivative(interactions_vec, interactions_location, interactions_scale);
        arma::vec log_prior_hessian = arma::join_cols(log_prior_hessian_thresholds,log_prior_hessian_interactions);
        arma::mat log_prior_hessian_mat = arma::diagmat(log_prior_hessian);
        hessian -= log_prior_hessian_mat;

        // adding prior information to HW
        arma::mat HW_inv = arma::inv_sympd(HW);
        HW_inv -= log_prior_hessian_mat;
        HW = arma::inv_sympd(HW_inv);

        // prior on the score
        //square_score += (gradient_prior * gradient_prior.t()); // we do not use this definition anywhere [[TO REMOVE]]
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("value") = loglik,
        Rcpp::Named("gradient") = gradient,
        Rcpp::Named("hessian") = hessian,
        Rcpp::Named("HW") = HW,
        Rcpp::Named("FisherInfo") = square_score);
  
   return out;
}



// function to calculate the negative pseudologlikelihood for a discrete MRF model
// [[Rcpp::export]]
double npseudologlik(
    const arma::vec &pars,
    const arma::mat &data, // this is already the matrix of sufficient statistics, by row (person) it looks like: {X_1, ..., X_P,2X_1X_2,...,2X_{P-1}X_P}
    const arma::uword &P,
    const arma::uvec &n_categories,
    const bool &with_prior, // whether to include (TRUE) or not (FALSE) prior information (non-informative priors are used)
    const int &ncores,
    const double &thresholds_alpha,
    const double &thresholds_beta,
    const double &interactions_location,
    const double &interactions_scale)
{
    arma::uword n,p,h;
    arma::uword n_pars = pars.n_elem; // this must be equal to P+P*(P-1)/2 or to pars.n_elem
    arma::uword N = data.n_rows;
    arma::uword n_thresholds = arma::accu(n_categories-1);
    arma::vec thresholds = pars(arma::span(0,n_thresholds-1));  // vector of thresholds parameters
    arma::vec interactions_vec = pars(arma::span(n_thresholds,n_pars-1));// vector of P*(P-1)/2 interaction parameters

    // building the symmetric matrix of interaction parameters from its lower triangular
    arma::mat lower_matrix_interactions(P,P,arma::fill::zeros);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(lower_matrix_interactions), -1); // element indices of the lower triangular excluding the diagonal elements 
    lower_matrix_interactions(lower_indices) = interactions_vec;
    arma::mat interactions = arma::symmatl(lower_matrix_interactions);

    // check on n_pars
    if((n_thresholds+P*(P-1)/2) != n_pars){ // check that data has the same columns as the number of parameters in the model
        Rcpp::Rcout << "The length of pars must be equal to " << (n_thresholds+P*(P-1)/2) << "\n";
        // stop algorithm (this check can be handled at R-level)
    }

    // creating empty objects where to save loglik, gradient and hessian computed per each person (statistical unit) , this is useful for the parallelization step
    arma::vec loglik_vec(N,arma::fill::zeros);

    // # pragma omp parallel for if(ncores>1) private() shared() ... parallelize here (??) the for loop over people
    // loop over people  
    #ifdef _OPENMP
    omp_set_dynamic(0);         
    omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for if(ncores>1) private(n,p,h,i,j) shared(N,P,n_thresholds,n_pars,data,n_categories,thresholds,interactions_vec,interactions,loglik_vec)
    #endif
    for(n = 0; n < N; n++){
        // select n-th person statistics
        arma::vec stats_n = data.row(n).t(); // stats for n-th person {X1,X2,...,2XiXj}
        // processing observed category per each X
        arma::uvec stats_X_n = arma::conv_to<arma::uvec>::from(stats_n(arma::span(0,P-1)));
        arma::vec pars_n(P+P*(P-1)/2,arma::fill::zeros);
        for(p = 0; p < P; p++){
            arma::uword which_threshold = stats_X_n(p);
            if(which_threshold > 0){
                arma::uword index_threshold_Xp = (which_threshold-1);
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                pars_n(p) = thresholds(index_threshold_Xp); // which_threshold-1 because in the thresholds vector we omit the baseline
            }
        }
        // filling in the interaction parameters
        pars_n(arma::span(P,pars_n.n_elem-1)) = interactions_vec;
        // calculating the numerator of the pseudo-loglikelihood for participant n
        double loglik_n = arma::accu(pars_n(arma::span(0,P-1))); // sum_p{sum_h{threshold_h*I(x == h)}}
        loglik_n += arma::accu((pars_n(arma::span(P,pars_n.n_elem-1)).t() * stats_n(arma::span(P,pars_n.n_elem-1)))); // sum_p{sum{2x_p*x_j*sigma_pj}

        //for each p we have to compute the support (normalizing constant, denominator) for each item, that is ln[1+sum_h{exp(mu_h+h*sum_{j!=p}{x_jsigma_pj})}]
        double support_n = 0.0;
        for(p = 0; p < P; p++){
            double denom_p = 1.0;
            arma::vec stats_excl_p = stats_n(arma::span(0,P-1));
            stats_excl_p(p) = 0.0; // this could be avoided because the interaction matrix has 0.0 in the diagonal
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                double success_event = arma::accu(thresholds(index_threshold_Xp) + static_cast<double>(h)*(stats_excl_p.t() * interactions.col(p)));                 
                // calculating denom_p
                denom_p += std::exp(success_event);
            }
            support_n += std::log(denom_p);
        }
        // updating loglikelihood
        loglik_n -= support_n;
        loglik_vec(n) -= loglik_n; // (-=) because negative loglikelihood          
    }
    
    // sum over people
    double loglik = arma::accu(loglik_vec); 

    // Adding prior information on loglik, gradient and hessian 
    if(with_prior){
        // prior on pseudologlikelihood
        loglik -= (arma::accu(log_beta_prime(thresholds, thresholds_alpha, thresholds_beta)) + arma::accu(log_dcauchy(interactions_vec, interactions_location, interactions_scale)));
    }
   return loglik;
}