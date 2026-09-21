#include <string>
#include <RcppArmadillo.h>
#include "utils.h"
#include "priors.h"
#include "paper-utils.h"
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

// [[Rcpp::export]]
void update_target_rescaling(
    const arma::mat &I,
    const arma::vec &pars,
    const arma::mat &data, // this is already the matrix of sufficient statistics, by row (person) it looks like: {X_1, ..., X_P,2X_1X_2,...,2X_{P-1}X_P}
    const arma::vec &frequency, // frequency of data by row
    const arma::uword &P,
    const arma::uword &N,
    const arma::uword &n_pars,
    const arma::uword &n_thresholds,
    const arma::uvec &n_categories,
    const arma::uvec &lower_indices,
    arma::mat &interactions,
    const arma::umat &matrix_indices_sigma,
    const arma::uvec &which_stats,
    const arma::vec &category_stats,
    const arma::uvec &category_offsets,
    arma::mat &Score,
    arma::mat &hessian,
    arma::mat &invHW,
    arma::mat &HW,
    arma::mat &Gamma,
    arma::mat &invGamma,
    arma::mat &Lt,
    arma::mat &invLt,
    arma::vec &log_prior_curvature,
    arma::mat &log_prior_curvature_mat
)
{
    arma::uword n,p,i,j,h;
    arma::vec thresholds = pars(arma::span(0,n_thresholds-1));  // vector of thresholds parameters
    arma::vec interactions_vec = pars(arma::span(n_thresholds,n_pars-1));// vector of P*(P-1)/2 interaction parameters

    // --- Convert thresholds interactions from vector to matrix form ---
    interactions_vec_to_mat(interactions, pars, lower_indices, P, n_thresholds, n_pars); // This overwrites matrix 'interactions'

    // Creating empty objects where to save gradient and hessian computed per each person (statistical unit)
    arma::vec probs_n(n_thresholds,arma::fill::zeros);
    arma::vec var_n(n_thresholds,arma::fill::zeros);

    arma::vec gradient_n(n_pars,arma::fill::zeros);
    arma::mat hessian_n(n_pars,n_pars,arma::fill::zeros);

    hessian.zeros();
    Score.zeros();

    // Loop over people - only computing Hessian and Score 
    for(n = 0; n < N; n++){
        gradient_n.zeros(); // resetting gradient vector for n-th sample
        hessian_n.zeros(); // resetting hessian matrix for n-th sample
        probs_n.zeros(); // resetting vector of probabilities for n-th sample
        var_n.zeros(); // resetting vector of variances for n-th sample

        // Select n-th person statistics
        arma::vec stats_n = data.row(n).t(); // stats for n-th unique state {X1,X2,...,2XiXj}
 
        // For each p we have to compute the support (normalizing constant, denominator) for each item, that is ln[1+sum_h{exp(mu_h+h*sum_{j!=p}{x_jsigma_pj})}] --> used in 'probs' and 'var'
        for(p = 0; p < P; p++){
            double denom_p = 1.0;
            arma::vec stats_excl_p = stats_n(arma::span(0,P - 1));
            stats_excl_p(p) = 0.0; // this could be avoided because the interaction matrix has 0.0 in the diagonal
            double xixj_sigma = arma::dot(stats_excl_p,interactions.col(p));
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h - 1;
                if(p>0){
                    index_threshold_Xp += category_offsets(p - 1); //arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                double success_event = thresholds(index_threshold_Xp) + static_cast<double>(h) * xixj_sigma;                 
                // Calculating denom_p, and probs_n numerator
                denom_p += std::exp(success_event);
                probs_n(index_threshold_Xp) = std::exp(success_event); // success probability for p-th variable
            }
        
            // Calculate pr(Xp = h) and variances 
            for(h = 1; h < n_categories(p); h++){
                arma::uword index_threshold_Xp = h - 1;
                if(p>0){
                    index_threshold_Xp += category_offsets(p - 1); //arma::accu(n_categories(arma::span(0,p-1))-1);
                }
                probs_n(index_threshold_Xp) /= denom_p;
                gradient_n(index_threshold_Xp) += static_cast<double>((stats_n(p) == h)) - probs_n(index_threshold_Xp); // gradient for thresholds = X - E(X)
                var_n(index_threshold_Xp) = probs_n(index_threshold_Xp) * (1.0 - probs_n(index_threshold_Xp)); // variance (for a bernoulli) p-th variable for level h
            }
        }
        

        // Precompute quantities before gradient and hessian loops: expected_X and expected_X_square for every node once before the Hessian double loop
        arma::vec expected_X(P, arma::fill::zeros), 
                  expected_X_sq(P, arma::fill::zeros);
        for(p = 0; p < P; p++){
            for(h = 1; h < n_categories(p); h++){
                arma::uword idx = h - 1;
                if(p > 0) {
                    idx += category_offsets(p-1);
                }
                expected_X(p)    += static_cast<double>(h) * probs_n(idx);
                expected_X_sq(p) += std::pow(static_cast<double>(h), 2) * probs_n(idx);
            }
        }
        arma::vec var_X = expected_X_sq - arma::square(expected_X);   // Var(X_p) for every node


        // Calculate gradient for interactions
        gradient_n(arma::span(n_thresholds,n_pars-1)) += stats_n(arma::span(P,P * (P - 1) / 2 + P - 1)); // summing first part of the gradient, that is the sufficient statistic 2XiXj

        // Now the (j,i) double loop just reads from expected_X: second part of the gradient for interactions, subtracting expected values
        arma::uword index_ij = n_thresholds;
        for(j = 0; j < (P-1); j++){
            for(i = (j+1); i < P; i++){
                gradient_n(index_ij) -= (stats_n(j) * expected_X(i) + stats_n(i) * expected_X(j));
                index_ij++;
            }
        }

        // Covariance matrix of score n-th person (Score contribute)
        Score += frequency(n) * (gradient_n * gradient_n.t()); // gr * gr.t() 
 
        // Calculate Hessian (cumulate the negative hessian_n, because we need it for the GHW)

        // hessian
        for(j = 0; j < n_pars; j++){ // column j
            for(i = j; i < n_pars; i++){ // row i
                if((j < n_thresholds) && (i < n_thresholds)){ // hessian for (thresholds) - only when i == j (diagonal elements), off diagonal elements remain 0.0
                    if(i == j){ // for (threshold_h,threshold_h) of the same X
                        hessian_n(i,j) -= var_n(j);
                    }
                    else if(which_stats(i) == which_stats(j)){ // for(threshold_h,threshold_k) of the same node X
                        hessian_n(i,j) += probs_n(i) * probs_n(j);
                        //hessian_n(j,i) = hessian_n(i,j);
                    }
                }
                else if((j < n_thresholds) && (i >= n_thresholds) && arma::any(matrix_indices_sigma.col(i-n_thresholds) == which_stats(j))){ // hessian for (thresholds,interactions) - only for (mu_k,sigma_{kl}) or (mu_l,sigma_{kl}), derivatives for (mu_k,sigma_{rl}) remain 0.0
                    arma::uword s = (!(matrix_indices_sigma(1,i - n_thresholds) == which_stats(j))) * 1; // selecting position index of the element different from which_stats(j)
                    arma::uword which_index = matrix_indices_sigma(s,i - n_thresholds);
                    hessian_n(i,j) -= stats_n(which_index) * probs_n(j) * (category_stats(j) - expected_X(which_stats(j))); 
                    //hessian_n(j,i) = hessian_n(i,j);
                } 
                else if((i >= n_thresholds) && (j >= n_thresholds)){ // hessian for (interactions) 
                    arma::uvec indices_i = matrix_indices_sigma.col(i-n_thresholds);
                    arma::uvec indices_j = matrix_indices_sigma.col(j-n_thresholds);
                    if(i == j){ // derivative for (sigma_{ij},sigma_{ij}) - here using either indices_j or indices_i is the same as they refer to the same sigma_{kl}
                        hessian_n(i,j) -= (std::pow(stats_n(indices_j(1)),2) * var_X(indices_j(0)) + std::pow(stats_n(indices_j(0)),2) * var_X(indices_j(1)));
                    }
                    else if((indices_i(0) == indices_j(0)|| indices_i(1) == indices_j(0)) || (indices_i(0) == indices_j(1)|| indices_i(1) == indices_j(1))){ // arma::any(indices_i == indices_j(0)) || arma::any(indices_i == indices_j(1))  // derivative for (sigma_{kl},sigma_{lr}), at this stage only one of the two conditions on indices_i can be true

                        /* [[OLD CODE for finding index_in_common, index_i and index_j]]
                        arma::uvec m = arma::join_cols(arma::find(indices_i == indices_j(0)),arma::find(indices_i == indices_j(1)));
                        arma::uword index_in_common =  indices_i(m(0));
                        arma::uvec which_indices_j = arma::find(indices_j != index_in_common);
                        arma::uword index_j = indices_j(which_indices_j(0));
                        arma::uvec which_indices_i = arma::find(indices_i != index_in_common);
                        arma::uword index_i = indices_i(which_indices_i(0));
                        */

                        // find the common element
                        arma::uword index_in_common = (indices_i(0) == indices_j(0) || indices_i(0) == indices_j(1)) 
                                                        ? indices_i(0) 
                                                        : indices_i(1);

                        // find the "other" element in indices_j
                        arma::uword index_j = (indices_j(0) == index_in_common)
                                                ? indices_j(1)
                                                : indices_j(0);

                        // find the "other" element in indices_i
                        arma::uword index_i = (indices_i(0) == index_in_common)
                                                ? indices_i(1)
                                                : indices_i(0);

                        hessian_n(i,j) -= stats_n(index_i) * stats_n(index_j) * var_X(index_in_common); 
                        //hessian_n(j,i) = hessian_n(i,j);
                    }
                }
            }
        }
        hessian_n = arma::symmatl(hessian_n); // make matrix symmetrical by reflecting lower-triangle to upper-triangle
        hessian -= frequency(n) * hessian_n; // (-=) because negative hessian
    } 

    // Adding prior information on inverse HW hessian 

    // prior curvature
    log_prior_curvature = arma::join_cols(log_beta_prime_second_derivative(thresholds),log_dcauchy_second_derivative(interactions_vec));
    log_prior_curvature_mat = arma::diagmat(log_prior_curvature);

    // adding prior information to HW and find Gamma
    //invScore = arma::inv_sympd(Score);
    invHW = (hessian * arma::solve(Score, hessian)) - log_prior_curvature_mat; // (-H) * S^{-1} * (-H) = (HW)^{-1} --> so we can include the prior curvature
    HW = arma::inv_sympd(invHW);
    Gamma = arma::chol(HW,"lower");

    // add prior information to hessian and find L^T
    hessian -= log_prior_curvature_mat;
    Lt = arma::chol(hessian);

    // update matrices to compute A^{-1}
    invGamma = arma::solve(arma::trimatl(Gamma), I);
    invLt = arma::solve(arma::trimatu(Lt), I); 

}


// AdaCoRe posterior sampler
// @param data matrix of size [N x (P + P(P-1)/2)] with the observed data (P columns) and the P(P-1)/2 double-crossproducts (lower triangular order)
// @param pars vector of size [n_pars] with the initial values of the parameters, in order [thresholds,interactions]
// @param n_categories vector of size [P] with the number of categories for each node
// @param pmles Pseudo-maximum likelihood estimates (pmles) 
// @param nsim number of iterations after burnin
// @param burnin number of burnin iterations
// @param adaptive_stage_n_iter number of iterations for the initial adaptive stage (simple MALA)
// @param sigma2 initial value for the step size (this is adaptive and a good starting value is needed (defailt is 1.0), however, it is allowed to vary over iterations) --> for OMRF is set to 0.01 (still needed?)
// @param thresholds_alpha hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param thresholds_beta hyperparameter for the beta prior on the thresholds (default is 0.5)
// @param interactions_location location parameter for the cauchy prior on the interactions (default is 0.0)
// @param interactions_scale scale parameter for the cauchy prior on the interactions (default is 2.5)
// @param verbose whether to print out the progress of each iteration (default is false). Also, verbose = true only if progress = false
// @param progress whether to show a progress bar (default is true)
// @return a list with the draws
// [[Rcpp::export]]
Rcpp::List cpp_adacore_sampler(const arma::mat &data,
                                const arma::vec &pars, 
                                const arma::uvec &n_categories,
                                const arma::vec &pmles, 
                                arma::uword nsim, 
                                arma::uword burnin, 
                                arma::uword adaptive_stage_n_iter = 500, 
                                double sigma2 = 1.0,
                                double thresholds_alpha = 0.5, 
                                double thresholds_beta = 0.5, 
                                double interactions_location = 0.0, 
                                double interactions_scale = 2.5,
                                bool verbose = false,
                                bool progress = true) {

    // --- Number of observations, variables, categories, and parameters ---
    arma::uword P = n_categories.n_elem; // number of variables (formerly data.n_cols but in this case data has P + P(P-1)/2 columns)
    arma::uword H = arma::max(n_categories - 1);
    arma::uword n_thresholds = arma::accu(n_categories - 1);
    arma::uword n_pars = pars.n_elem;

    // --- Number of total iterations (including burnin and adaptive stage) ---
    nsim += burnin + adaptive_stage_n_iter;
    arma::uword past_burnin = burnin + adaptive_stage_n_iter;
    arma::uword n_keep = nsim - past_burnin;

    // --- Find unique rows of the data and their frequencies ---
    arma::mat unique_data; // used for the update of the rescaling matrix (contains also the 2*XiXj statistics)
    arma::vec frequency;
    get_data_unique(unique_data, frequency, data);
    arma::uword N_unique = unique_data.n_rows; // For the pseudo posterior we need this sample size which is the number of unique rows in the data
    arma::mat unique_data_raw = unique_data.cols(0, P - 1); // used for the pseudo posterior gradient and logZ computation

    // --- Compute sufficient statistics for the observed data ---
    arma::vec obs_stats;
    sufficient_statistics_omrf(obs_stats, data.cols(0, P - 1), n_categories, P, n_thresholds, n_pars, 2.0);

    // --- Preallocate draws and sigma2_vec ---
    arma::mat draws(n_pars, n_keep, arma::fill::zeros);
    arma::vec sigma2_vec(nsim, arma::fill::zeros);

    // --- Parameters of adaptive proposal (FisherMALA) ---
    arma::mat R_n(n_pars, n_pars, arma::fill::eye);
    const arma::mat I_d(n_pars, n_pars, arma::fill::eye);   // diagonal matrix used in the adaptation of R_n
    double target_ar        = 0.574;   // optimal acceptance rate for MALA-type proposals (vs. 0.234 for RWM)
    double lambda           = 10.0;    // default value, Algorithm 1 in Fisher-adaptive MALA
    double learning_rate    = 0.015;   // default adaptation learning rate
    double sigma2_R         = sigma2; // this is  sigma2 / ((1/n_pars) * trace(R_n * R_n.t())), so sigma2_R = sigma2 at initial values
     
    // --- Parameters for updating the rescaling matrix ---
    arma::mat R_n_prev = R_n;
    double diff_norm = arma::norm(R_n - R_n_prev, "fro");
    double prev_norm = arma::norm(R_n_prev, "fro");
    double st = diff_norm / (prev_norm + 1e-12);
    double alpha_m = 0.05;
    double m = st;
    unsigned int counter_update = 0;
    double th_update = 3.0 / std::sqrt(static_cast<double>(data.n_rows));

    // --- Initialize counters ---
    arma::uword s           = 1,
                rejected    = 0,
                accepted    = 0;

    // --- Initialize structures used by compute_pseudo_gradient() ---

    // Lower triangular matrix indices excluding diagonal elements
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    // Utility [matrix 2 x P*(P-1)/2] where by column the indices (i,j) of the interaction effects 
    arma::umat matrix_indices_sigma;
    sigma_lower_tri_indices(matrix_indices_sigma, P);

    // Utility vectors indicating which stats and which category 
    arma::uvec which_stats;
    arma::vec category_stats;                           
    threshold_stat_index_map(which_stats, category_stats, n_categories, n_thresholds, P);

    // Utility vector indicating the offsets for each category in the thresholds vector
    arma::uvec category_offsets;
    category_offset_map(category_offsets, n_categories, P);

    // --- Setting starting values ---
    arma::vec eta_current = pars; // Original space
    sigma2_vec(s-1) = sigma2;

    // --- Preallocate matrices for Score and Hessian ---
    arma::mat Score(n_pars, n_pars), 
              hessian(n_pars, n_pars),
              interactions(P, P);

    // --- Preallocate matrices for the rescaling ---
    //
    // Three important relationships for the rescaling matrices:
    // 1. A = Gamma * Lt
    // 2. Ainv = invLt * invGamma --> A^{-1} = (L^T)^{-1} * Gamma^{-1}
    // 3. AinvT = invGamma.t() * invLt.t() --> A^{-T} = (A^{-1})^{T} = ( (L^T)^{-1} * Gamma^{-1} )^{T} = (Gamma^{-1})^{T} * ( (L^T)^{-1} )^{T}
    arma::mat invHW, HW, Gamma, invGamma, Lt, invLt, log_prior_curvature_mat;
    arma::vec log_prior_curvature;
    arma::vec pars_cumsum(n_pars, arma::fill::zeros);
    //             invHW(n_pars,n_pars,arma::fill::zeros),
    //             HW(n_pars,n_pars,arma::fill::zeros),
    //             Gamma(n_pars,n_pars,arma::fill::zeros),
    //             invGamma(n_pars,n_pars,arma::fill::zeros),
    //             Lt(n_pars,n_pars,arma::fill::zeros),
    //             invLt(n_pars,n_pars,arma::fill::zeros),
    //             log_prior_curvature_mat(n_pars,n_pars,arma::fill::zeros);
    // arma::vec log_prior_curvature(n_pars,arma::fill::zeros);

    update_target_rescaling(I_d,eta_current,unique_data,frequency,P,N_unique,n_pars,n_thresholds,n_categories,
                            lower_indices,interactions,matrix_indices_sigma,which_stats, category_stats,
                            category_offsets,Score,hessian,invHW,HW,Gamma,invGamma,Lt,invLt,log_prior_curvature,log_prior_curvature_mat);                       
    arma::vec beta_current = (Gamma * (Lt * (pars - pmles))) + pmles; // CoRe space                              

    // --- Calculate pseudo gradient and normalizing constant at current parameters, then trasnform to CoRe space gradient ---
    PseudoGradient core_current = compute_pseudo_gradient(obs_stats, unique_data_raw, frequency, eta_current, eta_current, 
                                                            n_categories, which_stats, lower_indices, category_offsets, P, N_unique, H, 
                                                            n_pars, n_thresholds, thresholds_alpha, thresholds_beta, 
                                                            interactions_location, interactions_scale);
    core_current.gradient =  invGamma.t() * (invLt.t() * core_current.gradient); // transform to CoRe space gradient

    // --- Start timer for the chain ---
    arma::wall_clock timer_chain;
    timer_chain.tic();

    // --- Progress bar setup (only if progress == true) ---
    arma::uword print_every = 0;
    Progress p(nsim, progress);
    if (progress) {
        print_every = std::max<arma::uword>(1, nsim / 100);
    }

    // --- Start sampling ---
    while(s < nsim){

        if (progress && Progress::check_abort()) {
            Rcpp::stop("Sampler aborted by user.");
        }

        // (0) Update cumulative means and adaptive A matrix 
        if (s < past_burnin) {
            pars_cumsum += eta_current;
            if(s > 1){ // for s = 1 we already initialized everything out of the loop
                diff_norm = arma::norm(R_n - R_n_prev, "fro");
                prev_norm = arma::norm(R_n_prev, "fro");
                st = diff_norm / (prev_norm + 1e-12);
                m = (1.0 - alpha_m) * m + alpha_m * st;
                if(m > th_update){
                    update_target_rescaling(I_d,pars_cumsum/static_cast<double>(s),unique_data,frequency,P,N_unique,n_pars,n_thresholds,n_categories,
                                            lower_indices,interactions,matrix_indices_sigma,which_stats,category_stats,category_offsets,Score,hessian,
                                            invHW,HW,Gamma,invGamma,Lt,invLt,log_prior_curvature,log_prior_curvature_mat);
                    counter_update++;
                    beta_current = (Gamma * (Lt * (eta_current - pmles))) + pmles; // CoRe space using the updated rescaling matrices
                }
                // --- Calculate pseudo gradient and normalizing constant at current parameters, then trasnform to CoRe space gradient ---
                core_current = compute_pseudo_gradient(obs_stats, unique_data_raw, frequency, eta_current, eta_current, 
                                                        n_categories, which_stats, lower_indices, category_offsets, P, N_unique, H, 
                                                        n_pars, n_thresholds, thresholds_alpha, thresholds_beta, 
                                                        interactions_location, interactions_scale);
                core_current.gradient =  invGamma.t() * (invLt.t() * core_current.gradient); // transform current parameters to CoRe space gradient (using the updated invGamma and invLt)
            }
        }

        // (1) Propose new parameters and transform to original space
        arma::vec beta_proposed = fishermala_propose(beta_current, core_current.gradient, R_n, sigma2_R); // Core space
        arma::vec eta_proposed = (invLt * (invGamma * (beta_proposed - pmles))) + pmles; // Original space

        // (2) Pseudo gradient and logZ at new parameters
        PseudoGradient core_proposed = compute_pseudo_gradient(obs_stats, unique_data_raw, frequency, eta_proposed, eta_current,
                                                                n_categories, which_stats, lower_indices, category_offsets, P, N_unique, H, 
                                                                n_pars, n_thresholds, thresholds_alpha, thresholds_beta, 
                                                                interactions_location, interactions_scale);
        core_proposed.gradient = invGamma.t() * (invLt.t() * core_proposed.gradient); // transform to CoRe space gradient

        // (3) Compute acceptance ratio
        double log_a = fishermala_core_log_acceptance_ratio(beta_current, beta_proposed, eta_current, eta_proposed, 
                                                            core_current.gradient, core_proposed.gradient, R_n, obs_stats, 
                                                            core_proposed.logZ_ratio, sigma2_R, n_thresholds, n_pars, thresholds_alpha, 
                                                            thresholds_beta, interactions_location, interactions_scale);

        // (4)  Adaptive step for R_n (FisherMALA)
        fishermala_update_preconditioner(R_n, log_a, core_current.gradient, core_proposed.gradient, I_d, lambda, s, adaptive_stage_n_iter);

        // (5) Adapt and normalize step size sigma2 and sigma2_R (FisherMALA)
        fishermala_update_step_size(sigma2, sigma2_R, R_n, n_pars,log_a, target_ar, learning_rate);
        // save adaptive sigma2 parameter
        sigma2_vec(s) = sigma2;   

        // (6) MH step: accept or reject the proposed parameters
        double u = R::runif(0.0, 1.0);

        // if(u < std::exp(log_a)) then set beta_current = beta_proposed update counter s++
        if(u < std::exp(log_a)){
            beta_current                        = beta_proposed;
            eta_current                         = eta_proposed;
            R_n_prev                            = R_n;
            if(s >= past_burnin){
                core_current.gradient           = core_proposed.gradient;
                draws.col(s - past_burnin)      = beta_proposed;
                accepted++;
            }
        }
        else{
            if(s >= past_burnin){
                draws.col(s - past_burnin) = beta_current;
                rejected++;
            }
        }

        // --- Print progress if verbose == true ---
        if (verbose && !progress) {
            Rcpp::Rcout << "(" << accepted << " of " << (accepted + rejected) << ") || iteration " << s
                        << " || acceptance = " << std::exp(log_a)
                        << " || sigma2 (" << sigma2 << ") with (de)increment = "
                        << learning_rate * (std::exp(log_a) - target_ar) << "\n";
        }

        if (s % print_every == 0) p.increment(print_every); // Update progress bar

        s++; // Increment counter
    }

    if (progress) {
        Rcpp::Rcout << "\n";
    }

    double time = timer_chain.toc(); // time elapsed (in seconds)

    // --- Calculate acceptance probability --- s / (s + rejected)
    double acceptance_probability = static_cast<double>(accepted)/(static_cast<double>(accepted) + static_cast<double>(rejected));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("draws") = draws.t(),
        Rcpp::Named("acceptance") = acceptance_probability,
        Rcpp::Named("seconds_elapsed") = time,
        Rcpp::Named("counter_update") = counter_update,
        Rcpp::Named("sigma2") = sigma2_vec
    );

    return out;
}

