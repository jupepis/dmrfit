#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects 
// [[Rcpp::export]]
arma::umat get_indices_sigma(const arma::uword &P){
    arma::umat matrix_indices_sigma(2,P*(P-1)/2);
    arma::uword l = 0;
    for(arma::uword j = 0; j < (P-1); j++){ // column j
        for(arma::uword i = (j+1); i < P; i++){ // row i
            matrix_indices_sigma(0,l) = i;
            matrix_indices_sigma(1,l) = j;
            l++;
        }
    }
    return matrix_indices_sigma;
}

// utility vector of which stats are involved in the thresholds (of length n_thresholds - used in the hessian computation)
// [[Rcpp::export]]
arma::uvec get_which_stats(const arma::uvec &n_categories, const arma::uword &n_thresholds){
    arma::uword P = n_categories.n_elem;
    arma::uvec which_stats(n_thresholds,arma::fill::zeros);
    arma::uword l = 0; 
    for(arma::uword p = 0; p < P; p++){     
        for(arma::uword h = 1; h < n_categories(p); h++){
            which_stats(l) = p;
            l++;
        }
    }
    return which_stats;
}

// utility vector of category statistics (of length n_thresholds - used in the hessian computation)
// [[Rcpp::export]]
arma::vec get_category_stats(const arma::uvec &n_categories, const arma::uword &n_thresholds){
    arma::uword P = n_categories.n_elem;
    arma::vec category_stats(n_thresholds,arma::fill::zeros);
    arma::uword l = 0; 
    for(arma::uword p = 0; p < P; p++){
        for(arma::uword h = 1; h < n_categories(p); h++){
            category_stats(l) = static_cast<double>(h);
            l++;
        }
    }
    return category_stats;
}

// sufficient statistics for omrf
// [[Rcpp::export]]
arma::vec get_sufficient_stats_omrf(
    const arma::mat &data, // P x N matrix of data (each column is a person, each row is a variable)
    const arma::uvec &which_stats,
    const arma::vec &category_stats,
    const arma::umat &matrix_indices_sigma,
    const arma::uword &n_thresholds,
    const arma::uword &n_pars
){
    arma::vec obs_stats(n_pars,arma::fill::zeros);
    for(arma::uword i = 0; i < n_pars; i++){
        if(i < n_thresholds){
            // calculate threshold statistic for the auxiliary data
            arma::uword which_stats_idx = which_stats(i);
            double category_stats_i = category_stats(i);
            obs_stats(i) = arma::accu(data.row(which_stats_idx) == category_stats_i);
        }
        else{
            // calculate interaction statistic for the auxiliary data
            arma::uvec ij = matrix_indices_sigma.col(i-n_thresholds);
            obs_stats(i) = 2.0 * arma::accu(data.row(ij(0)) % data.row(ij(1))); // * 2.0 because pseudo-likelihood parametrization
        }
    }
    return obs_stats;
}

// calculate log proposal ratio
// [[Rcpp::export]]
double get_log_proposal_ratio(
    const arma::vec &current_pars,
    const arma::vec &proposed_pars,
    const arma::vec &gradient_pseudo_current,
    const arma::vec &gradient_pseudo_proposed,
    const double &sigma2_R,
    const arma::mat &R_n) {

    arma::vec step_current = current_pars - proposed_pars - ((sigma2_R / 4.0) * (R_n * (R_n.t() * gradient_pseudo_proposed)));
    double h_current_pars = arma::accu(step_current.t() * gradient_pseudo_proposed);
    arma::vec step_proposed = proposed_pars - current_pars - ((sigma2_R / 4.0) * (R_n * (R_n.t() * gradient_pseudo_current)));
    double h_proposed_pars = arma::accu(step_proposed.t() * gradient_pseudo_current);
    return 0.5 * (h_current_pars - h_proposed_pars);
}

// caculate log-likelihood ratio (fast way to calculate without inverting vcov matrix - based on Proposition 1 from Michalis Titsias)
// [[Rcpp::export]]
double get_log_prior_ratio(
    const arma::vec &theta_current, 
    const arma::vec &theta_proposed,
    const arma::uword &n_thresholds, 
    const arma::uword &n_pars,
    const double &thresholds_alpha,
    const double &thresholds_beta,
    const double &interactions_location,
    const double &interactions_scale
){
    return arma::accu(
            log_beta_prime(theta_proposed(arma::span(0,n_thresholds-1)), thresholds_alpha, thresholds_beta) - log_beta_prime(theta_current(arma::span(0,n_thresholds-1)), thresholds_alpha, thresholds_beta)
        ) + 
        arma::accu(
            log_dcauchy(theta_proposed(arma::span(n_thresholds,n_pars-1)), interactions_location, interactions_scale) - log_dcauchy(theta_current(arma::span(n_thresholds,n_pars-1)), interactions_location, interactions_scale)
        ); 
}

// update preconditioning matrix R_n based on the current acceptance rate
void update_optimal_preconditioner(
    const arma::uword &s,
    const arma::uword &adaptive_stage_n_iter,
    const double &lambda,
    const arma::uword &n_pars,
    const arma::vec &gradient_pseudo_proposed,
    const arma::vec &gradient_pseudo_current,
    const double &log_a,
    const double &target_ar,
    arma::mat &R_n
){
    if(s - adaptive_stage_n_iter == 0){
        double r_proposed = 1.0 / (1.0 + std::sqrt(lambda / (lambda + arma::accu(gradient_pseudo_proposed.t() * gradient_pseudo_proposed))));
        arma::mat I_d = arma::eye<arma::mat>(n_pars,n_pars);
        R_n = (1.0 / std::sqrt(lambda)) * (I_d - r_proposed * (gradient_pseudo_proposed * gradient_pseudo_proposed.t()) / (lambda + (arma::accu(gradient_pseudo_proposed.t() * gradient_pseudo_proposed)))); // Proposition 4 Formula (12)
    }
    else if(s > adaptive_stage_n_iter){
        // Compute adaptation signal
        arma::vec adapt_signal = std::sqrt(std::exp(log_a)) * (gradient_pseudo_proposed - gradient_pseudo_current);
        arma::vec phi_proposed = R_n.t() * adapt_signal;
        double r_proposed =  1.0 / (1.0 + std::sqrt(1.0 / (1.0 + arma::accu(phi_proposed.t() * phi_proposed))));
        R_n -= r_proposed * ( ((R_n * phi_proposed) * phi_proposed.t())/(1 + (arma::accu(phi_proposed.t() * phi_proposed))) ); // Proposition 4 Formula (13)
    }
}

// update global step size sigma2 and normalize it (sigma2_R)
void update_global_step_size(
    double &sigma2,
    double &sigma2_R,
    const arma::mat &R_n,
    const double &log_a,
    const double &target_ar,
    const double &learning_rate,
    const arma::uword &n_pars){
    // Adapt
    sigma2 = sigma2 * (1.0 + learning_rate * (std::exp(log_a) - target_ar));
    // Normalize
    sigma2_R = sigma2 / ((1 / static_cast<double>(n_pars)) * arma::trace(R_n * R_n.t())); // check if arma::accu(R % R) is faster than arma::trace(R * R.t())
}

// gradient vector and log_Z_ratio Ordinal MRF Pseudolikelihood
// [[Rcpp::export]]
arma::field<arma::vec> omrf_pl_grad_logZ(
    const arma::vec &pars,
    const arma::vec &current_pars,
    const arma::mat &data, // P x N matrix of data (each column is a person, each row is a variable)
    const arma::vec &obs_stats,
    const arma::uword &P,
    const arma::umat &matrix_indices_sigma,
    const arma::uvec &which_stats,
    const arma::vec &category_stats,
    const arma::uvec &n_categories,
    const arma::uvec &lower_indices,
    const arma::uword &n_thresholds,
    const arma::uword &n_pars,
    const int &ncores,
    const bool &with_prior,
    const double &thresholds_alpha,
    const double &thresholds_beta,
    const double &interactions_location,
    const double &interactions_scale
)
{
    arma::vec thresholds = pars(arma::span(0,n_thresholds-1)); 
    arma::vec thresholds_current = current_pars(arma::span(0,n_thresholds-1)); 

    // interactions (proposed): building the symmetric matrix of interaction parameters from its lower triangular
    arma::mat lower_matrix_interactions(P,P,arma::fill::zeros);
    lower_matrix_interactions(lower_indices) = pars(arma::span(n_thresholds,n_pars-1));
    arma::mat interactions = arma::symmatl(lower_matrix_interactions);

    // interactions_current: building the symmetric matrix of interaction parameters from its lower triangular
    arma::mat lower_matrix_interactions_current(P,P,arma::fill::zeros);
    lower_matrix_interactions_current(lower_indices) = current_pars(arma::span(n_thresholds,n_pars-1));
    arma::mat interactions_current = arma::symmatl(lower_matrix_interactions_current);

    // creating empty objects where to save log_Z_ratio and gradient
    arma::vec log_Z_ratio(1,arma::fill::zeros); 
    arma::vec gradient(n_pars,arma::fill::zeros);
    arma::vec gradient_n(n_pars,arma::fill::zeros); // this vector is initialized here to avoid re-declaring it in the loop over n, it is reset to 0 at each iteration of the loop
    arma::vec probs_n(n_thresholds,arma::fill::zeros); // vector of probabilities  P(Xp = h)


    // for loop over observations (rowise access to data) 
    for(arma::uword n = 0; n < data.n_cols; n++){
        // select the statistics for the n-th observation
        arma::vec stats_n = data.col(n); // stats for n-th person {X1,X2,...,XP}

        //for each p we have to compute the support (normalizing constant, denominator) for each item, that is ln[1+sum_h{exp(mu_h+h*sum_{j!=p}{x_jsigma_pj})}]
        probs_n.zeros(); // reset probs_n for each observation
        for(arma::uword p = 0; p < P; p++){
            double denom_p = 1.0;
            double denom_p_current_pars = 1.0;
            arma::vec stats_excl_p = stats_n;
            stats_excl_p(p) = 0.0; // [verify that the comment on the right can be done and apply] this could be avoided because the interaction matrix has 0.0 in the diagonal
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }

                // at the proposed pars: pars
                double success_event = arma::accu(thresholds(index_threshold_Xp) + static_cast<double>(h)*(stats_excl_p.t() * interactions.col(p)));                 
                // calculating denom_p, and probs_n numerator for proposed pars
                denom_p += std::exp(success_event);
                probs_n(index_threshold_Xp) = std::exp(success_event); // success probability for p-th variable

                // at the current pars: current_pars
                double success_event_current_pars = arma::accu(thresholds_current(index_threshold_Xp) + static_cast<double>(h)*(stats_excl_p.t() * interactions_current.col(p))); 
                // calculating denom_p for current pars
                denom_p_current_pars += std::exp(success_event_current_pars);

            }

            // calculate pr(Xp = h)  
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                probs_n(index_threshold_Xp) /= denom_p; // divide probs_n for the calculated denom_p
            }

            log_Z_ratio(0) += std::log(denom_p_current_pars) - std::log(denom_p);
        }

        // calculate (negative) gradient
        gradient_n.zeros();

        // gradient for thresholds
        for(arma::uword p = 0; p < P; p++){
            for(arma::uword h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h-1;
                if(p>0){
                    index_threshold_Xp += arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                gradient_n(index_threshold_Xp) -= probs_n(index_threshold_Xp);
            }
        }

        // gradient for interactions
        arma::uword index_ij = n_thresholds;
        for(arma::uword j = 0; j < (P-1); j++){
            // calculate expected Xj
            double expected_Xj = 0.0;
            for(arma::uword h = 1; h < n_categories(j); h++){
                arma::uword index_threshold_Xj = h-1;
                if(j>0){
                    index_threshold_Xj += arma::accu(n_categories(arma::span(0,j-1))-1);
                }
                expected_Xj += static_cast<double>(h)*probs_n(index_threshold_Xj);
            }
            for(arma::uword i = (j+1); i < P; i++){
                // calculate expected Xi
                double expected_Xi = 0.0;
                for(arma::uword h = 1; h < n_categories(i); h++){
                    arma::uword index_threshold_Xi = h-1;
                    if(i>0){
                        index_threshold_Xi += arma::accu(n_categories(arma::span(0,i-1))-1);
                    }
                    expected_Xi += static_cast<double>(h)*probs_n(index_threshold_Xi);
                }
                gradient_n(index_ij) -= (stats_n(j)*expected_Xi + stats_n(i)*expected_Xj); // subtracting the second part that is: - Xj * \sum_h{h * P(Xi=h)} - Xi * \sum_h{h(P(Xj=h)}
                index_ij++;
            }
        }

        gradient += gradient_n; 
        
    }

    
    // add sufficient statistics S(X) to the gradient (gradient = S(X) - E_{theta}[S(X)])
    gradient += obs_stats;

    if(with_prior){
        // log_prior 
        arma::vec log_prior(n_pars,arma::fill::zeros);
        log_prior(arma::span(0,n_thresholds-1)) = log_beta_prime_first_derivative(pars(arma::span(0,n_thresholds-1)), thresholds_alpha, thresholds_beta);
        log_prior(arma::span(n_thresholds,n_pars-1)) = log_dcauchy_first_derivative(pars(arma::span(n_thresholds,n_pars-1)), interactions_location, interactions_scale);
        gradient += log_prior;
    }

    arma::field<arma::vec> out_field(2);
    out_field(0) = gradient;
    out_field(1) = log_Z_ratio;

    return out_field;
}

// [[Rcpp::export]]
Rcpp::List omrf_core_sampler(
  const arma::mat &data, // matrix of P rows (variables: X_1, ..., X_P) and N columns (each column is a person)
  const arma::vec &pars, // (thresholds, interactions)
  const arma::uvec &n_categories, // number of categories per each node
  const arma::uword &P, // number of nodes (items)
  arma::uword &nsim, //number of iterations after burnin
  const arma::uword &burnin, // burnin iterations
  const arma::vec &pmles, // pmles or robbins-monro adjusted estimates
  const arma::mat &current_scale, //  cholesky of variance and covariance matrix of pseudo-likelihood calculate ad the 'pmles' 
  const arma::mat &new_scale, // transposed cholesky of inverse variance and covariance matrix - MC-hessian, Huber-White or Robbins-Monro's
  const arma::uword &adaptive_stage_n_iter = 500, // iterations for the initial adaptive stage (simple MALA)
  double sigma2 = 1.0, // for ordinal MRF, set this value to [[0.01]] or higher 0.1, by default this is set to 1.0 but with the ordinal MRF the algorithms runs fine with 0.1 or 0.01 (this sigma2 is adaptive and a good starting value is needed, however, it is allowed to vary over iterations)
  const double &thresholds_alpha = 0.5, 
  const double &thresholds_beta = 0.5, 
  const double &interactions_location = 0.0, 
  const double &interactions_scale = 2.5
){
    arma::uword s = 1,
                step_s = 0,
                rejections = 0,
                accepted = 0;

    // for the function omrf_pl_grad_logZ
    bool with_prior = true; 
    int ncores = 1; 

    // processing input arguments
    arma::uword n_thresholds = arma::accu(n_categories-1);
    arma::uword n_pars = pars.n_elem;
    nsim += burnin + adaptive_stage_n_iter;
    arma::mat draws(n_pars,nsim,arma::fill::zeros);

    // precomputing matrices for the coordinate rescaling step
    arma::mat A = new_scale * current_scale;
    arma::mat Ainv = arma::inv(A);
    arma::mat AinvT = Ainv.t();

    // lower triangular indices
    arma::mat lower_matrix_interactions(P,P,arma::fill::zeros);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(lower_matrix_interactions), -1); // lower triangular matrix indices excluding diagonal elements [this can be passed by reference to the gradient_and_log_Z_ratio_... function]
    
    // utility matrix indicating which variables are involved in the interaction effects (of size 2 x P*(P-1)/2)
    arma::umat matrix_indices_sigma = get_indices_sigma(P); 

    // utility vectors indicating which stats and which category (of length n_thresholds - used in the hessian computation)
    arma::uvec which_stats = get_which_stats(n_categories, n_thresholds); 
    arma::vec category_stats = get_category_stats(n_categories, n_thresholds);

    // sufficient statistics for omrf
    arma::vec obs_stats = get_sufficient_stats_omrf(data, which_stats, category_stats, matrix_indices_sigma, n_thresholds, n_pars); 

    // setting initial values for the adaptive MALA algorithm
    arma::vec mu(n_pars,arma::fill::zeros); // vector of mean zero for the random normal perturbation in the proposal distribution
    arma::mat vcov_mat = arma::eye<arma::mat>(n_pars,n_pars); // variance-covariance matrix for the random normal perturbation in the proposal distribution
    double target_ar = 0.574; // target acceptance rate for MALA with Fisher's optimal preconditioning matrix
    double lambda = 10.0; // default value taken by Algorithm 1 in Fisher adaptive MALA
    double learning_rate = 0.015; // default value for cost learning rate
    arma::mat R_n = arma::eye<arma::mat>(n_pars,n_pars); // square root of the preconditioning matrix, initialized to the identity matrix and updated over the iterations
    double sigma2_R = sigma2;  // normalized global step size parameter, initialize sigma2_R with the value of sigma2 (default is 1.0), inside the sampler sigma2_R = sigma2 / ((1/static_cast<double>(n_pars)) * arma::trace(R_n * R_n.t()))

    // starting the chain from the pmles or from the robbins-monro adjusted estimates
    arma::vec current_pars = (A * (pars - pmles)) + pmles;
    draws.col(s-1) = current_pars;
    arma::vec theta_current = pars;
    arma::field<arma::vec> gradient_pseudo_current_theta = omrf_pl_grad_logZ(theta_current,theta_current,data,obs_stats,P,matrix_indices_sigma,which_stats,category_stats,n_categories,lower_indices,n_thresholds,n_pars,ncores,with_prior,thresholds_alpha,thresholds_beta,interactions_location,interactions_scale);
    arma::mat gradient_pseudo_current = AinvT * gradient_pseudo_current_theta(0);

    // while s <  nsim
    bool burnin_started = false;
    bool sampling_started = false;
    Progress pb(nsim, true);
    REprintf("  [Adaptive stage] (%d iterations)\n", adaptive_stage_n_iter);
    arma::wall_clock timer_chain;
    timer_chain.tic();
    while(s < nsim){
        if(Progress::check_abort()) return Rcpp::List();

        // print phase messages before incrementing
        if(!burnin_started && s == adaptive_stage_n_iter) {
            burnin_started = true;
            REprintf("  [Burn-in] (%d iterations)\n", burnin);
        }
        if(!sampling_started && s == adaptive_stage_n_iter + burnin) {
            sampling_started = true;
            REprintf("  [Sampling] (%d iterations)\n", nsim);
        }

        pb.increment();

        // (1) Propose new parameters
        arma::vec eta_n = Rcpp::rnorm(n_pars); //arma::mvnrnd(mu,vcov_mat);
        arma::vec proposed_pars = current_pars + (sigma2_R / 2.0) * (R_n * (R_n.t() * gradient_pseudo_current)) + std::sqrt(sigma2_R) * (R_n * eta_n);

        // (2) Exact gradient at new parameters
        arma::vec theta_proposed = (Ainv * (proposed_pars - pmles)) + pmles;
        arma::field<arma::vec> gradient_pseudo_proposed_theta = omrf_pl_grad_logZ(theta_proposed,theta_current,data,obs_stats,P,matrix_indices_sigma,which_stats,category_stats,n_categories,lower_indices,n_thresholds,n_pars,ncores,with_prior,thresholds_alpha,thresholds_beta,interactions_location,interactions_scale);
        arma::mat gradient_pseudo_proposed = AinvT * gradient_pseudo_proposed_theta(0);
        
        // (3) Compute acceptance ratio

        // (3.1) compute log_proposal_ratio 
        double log_proposal_ratio = get_log_proposal_ratio(current_pars, proposed_pars, gradient_pseudo_current, gradient_pseudo_proposed, sigma2_R, R_n); 
    
        // (3.2) compute loglik_data_ratio
        double loglik_data_ratio = arma::accu((theta_proposed - theta_current) % obs_stats) + arma::as_scalar(gradient_pseudo_proposed_theta(1)); // this is equivalent to the log-likelihood ratio because the log_Z_ratio is already included in the gradient_pseudo_proposed_theta(1) (see the output of the function omrf_pl_grad_logZ)

        // (3.3)  compute log_prior_ratio
        double log_prior_ratio = get_log_prior_ratio(theta_proposed, theta_current, n_thresholds, n_pars, thresholds_alpha, thresholds_beta, interactions_location, interactions_scale);

        // log acceptance ratio 
        double log_a = log_prior_ratio + loglik_data_ratio + log_proposal_ratio;
        if(log_a > 0.0){ 
            log_a = 0.0; // equivalent to log min{1,acceptance_ratio}
        }

        // (4) update preconditioning matrix R_n
        update_optimal_preconditioner(s, adaptive_stage_n_iter, lambda, n_pars, gradient_pseudo_proposed, gradient_pseudo_current, log_a, target_ar, R_n);

        // (5) Adapt and normalize step size sigma2
        update_global_step_size(sigma2, sigma2_R, R_n, log_a, target_ar, learning_rate, n_pars);

        // (6) draw random U(0,1) and decide whether to accept or reject current draw 
        double u = Rcpp::as<double>(Rcpp::runif(1));

        // if(u < std::exp(log_a)) then set current_pars = proposed_pars update counter s++
        if(u < std::exp(log_a)){
            current_pars = proposed_pars; 
            theta_current = theta_proposed;
            gradient_pseudo_current = gradient_pseudo_proposed;
            draws.col(s) = proposed_pars;
            if(s >= (burnin + adaptive_stage_n_iter)){
                accepted++;
            }
        }
        else{
            draws.col(s) = current_pars; 
            if(s > burnin + adaptive_stage_n_iter){
                rejections++;
            }
        }

        if(s >= burnin + adaptive_stage_n_iter){ 
            step_s++;
        }

        //Rcpp::Rcout << "("<< accepted <<" of " << step_s <<" - acceptance = " << std::round((static_cast<double>(accepted)/static_cast<double>(step_s)) * 1000.0) / 1000.0 << ") || iteration " << s <<  " || a = " << std::exp(log_a) <<" || sigma2 (" << sigma2 << ") with (de)increment = " << learning_rate*(std::exp(log_a)-target_ar) <<"\n";

        s++;
    }
    double seconds_elapsed_chain = timer_chain.toc();

    // Acceptance probability s/(s+rejections)
    double acceptance_probability = static_cast<double>(accepted)/(static_cast<double>(accepted)+static_cast<double>(rejections)); //static_cast<double>(n_iter)/(static_cast<double>(n_iter)+static_cast<double>(rejections));

    // Output
    arma::mat output_draws = draws.cols(arma::span(burnin+adaptive_stage_n_iter,nsim-1));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("draws") = output_draws,
        Rcpp::Named("acceptance") = acceptance_probability,
        Rcpp::Named("seconds_elapsed") = seconds_elapsed_chain
    );
    return out;
}