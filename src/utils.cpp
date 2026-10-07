#include <string>
#include <map>
#include <vector>
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
arma::mat cpp_mvnrnd_arma(const arma::vec &mu, const arma::mat &Sigma, int n) {
  int p = mu.n_elem;
  arma::mat Z = rnorm_arma(p, n);   // each column is a sample
  arma::mat C = arma::chol(Sigma, "lower");
  arma::mat X = C * Z; // from standard to correlated normal
  X.each_col() += mu;
  return X;
}

// Utils Discrete MRF Pseudo-likelihood 
//
// pseudo-likelihood of a discrete mrf function value and derivatives value at specific parameters (used by the optimization algortihm)
// the function returns a list of named objects
// "loglik" is the negative pseudologlikelihood calculated at the parameters supplied via the input pars
// "gradient" is the negative gradient value at pars
// "hessian" is the negative hessian value (its inverse returns the matrix of variances and covariances of the model parameters)
// "HW" is the Godambe-Huber-White (GHW) covariance (it is already a matrix of variances and covariances)
// note: loglik, gradient and hessian are used by the trust algorithm to find the MPLEs
// [[Rcpp::export]]
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
    arma::uword n_pars = pars.n_elem; // this must be equal to n_thresholds + P * (P - 1) / 2
    arma::uword n_thresholds = arma::accu(n_categories - 1);
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

    // --- Positions of the parameters of the full conditional of each node p: its thresholds and its P - 1 interactions ---
    // category_offsets(p) is the position of the first threshold of node p; index_interaction(p, j) is the position of the
    // interaction between nodes p and j (j != p), with the interactions in lower-triangular column order
    arma::uvec category_offsets(P, arma::fill::zeros);
    for(arma::uword p = 1; p < P; p++){
        category_offsets(p) = category_offsets(p - 1) + n_categories(p - 1) - 1;
    }
    arma::umat index_interaction(P, P, arma::fill::zeros);
    for(arma::uword j = 0; j < P; j++){
        for(arma::uword i = j + 1; i < P; i++){
            arma::uword index_ij = n_thresholds + j * P - j * (j + 1) / 2 + (i - j - 1);
            index_interaction(i, j) = index_ij;
            index_interaction(j, i) = index_ij;
        }
    }

    // --- Unique response patterns and their frequencies: each pattern is processed once ---
    arma::mat data_unique;
    arma::vec frequency;
    get_data_unique(data_unique, frequency, data.cols(0, P - 1));
    arma::uword N_unique = data_unique.n_rows;

    // --- Split the patterns into a fixed number of chunks, independent of ncores, summed in order after the parallel loop ---
    // (the results do not depend on the number of threads; no BLAS routine is called inside the parallel loop, because BLAS
    // libraries with their own threads are not safe to call from several threads at once)
    arma::uword n_chunks = std::min<arma::uword>(16, N_unique);
    arma::uvec chunk_start(n_chunks + 1);
    for(arma::uword c = 0; c <= n_chunks; c++){
        chunk_start(c) = (c * N_unique) / n_chunks;
    }
    arma::vec loglik_chunk(n_chunks, arma::fill::zeros);
    arma::mat gradient_chunk(n_pars, n_chunks, arma::fill::zeros);
    arma::cube hessian_chunk(n_pars, n_pars, n_chunks, arma::fill::zeros);
    arma::mat gradient_patterns(n_pars, N_unique, arma::fill::zeros); // negative gradient of one observation of each pattern

    // --- Loop over the chunks, in parallel over ncores threads ---
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(static)
    #endif
    for(arma::uword c = 0; c < n_chunks; c++){
        arma::mat& hessian_c = hessian_chunk.slice(c);

        for(arma::uword n = chunk_start(c); n < chunk_start(c + 1); n++){
            arma::vec stats_n = data_unique.row(n).t(); // stats {X1, X2, ..., XP} of the n-th pattern
            double frequency_n = frequency(n);

            // --- Weighted sum of the other variables of each node: sum_{j != p} sigma_pj x_j (as a loop, no BLAS call) ---
            arma::vec xixj_sigma(P, arma::fill::zeros);
            for(arma::uword p = 0; p < P; p++){
                for(arma::uword j = 0; j < P; j++){
                    xixj_sigma(p) += interactions(j, p) * stats_n(j); // the diagonal of interactions is 0.0
                }
            }
            double loglik_n = 0.0;
            arma::vec gradient_n(n_pars, arma::fill::zeros);

            for(arma::uword p = 0; p < P; p++){
                arma::uword n_thresholds_p = n_categories(p) - 1;
                arma::uword x_p = static_cast<arma::uword>(stats_n(p));

                // --- Conditional probabilities P(Xp = h), h = 1, ..., m - 1 (log-sum-exp with the baseline category 0) ---
                arma::vec success_event(n_thresholds_p);
                for(arma::uword h = 1; h <= n_thresholds_p; h++){
                    success_event(h - 1) = thresholds(category_offsets(p) + h - 1) + static_cast<double>(h) * xixj_sigma(p);
                }
                double max_event = std::max(0.0, success_event.max());
                arma::vec probs_p = arma::exp(success_event - max_event);
                double denom_p = std::exp(-max_event) + arma::accu(probs_p);
                probs_p /= denom_p;
                double log_denom_p = max_event + std::log(denom_p);

                // --- Conditional expected value and variance of Xp ---
                double expected_X = 0.0;
                double expected_X_square = 0.0;
                for(arma::uword h = 1; h <= n_thresholds_p; h++){
                    expected_X += static_cast<double>(h) * probs_p(h - 1);
                    expected_X_square += static_cast<double>(h * h) * probs_p(h - 1);
                }
                double variance_X = expected_X_square - expected_X * expected_X;

                // --- Log pseudo-likelihood of node p: mu_{p, x_p} + x_p sum_{j != p} sigma_pj x_j - log(denom_p) ---
                if(x_p > 0){
                    loglik_n += thresholds(category_offsets(p) + x_p - 1);
                }
                loglik_n += stats_n(p) * xixj_sigma(p) - log_denom_p;

                // --- Gradient of the log pseudo-likelihood: I(Xp = h) - P(Xp = h) for the thresholds of node p ---
                for(arma::uword h = 1; h <= n_thresholds_p; h++){
                    gradient_n(category_offsets(p) + h - 1) += static_cast<double>(x_p == h) - probs_p(h - 1);
                }
                // --- and x_j (x_p - E[Xp]) for the interactions of node p ---
                for(arma::uword j = 0; j < P; j++){
                    if(j != p && stats_n(j) != 0.0){
                        gradient_n(index_interaction(p, j)) += stats_n(j) * (stats_n(p) - expected_X);
                    }
                }

                // --- Negative Hessian: frequency times the covariance of the sufficient statistics of node p under its full
                // conditional: Cov(mu_h, mu_k) = P(Xp = h) (I(h = k) - P(Xp = k)), Cov(mu_h, sigma_pj) = x_j P(Xp = h) (h - E[Xp]),
                // Cov(sigma_pj, sigma_pk) = x_j x_k Var(Xp) ---
                for(arma::uword h = 1; h <= n_thresholds_p; h++){
                    arma::uword index_threshold_h = category_offsets(p) + h - 1;
                    for(arma::uword k = 1; k <= n_thresholds_p; k++){
                        hessian_c(index_threshold_h, category_offsets(p) + k - 1) += frequency_n * probs_p(h - 1) * (static_cast<double>(h == k) - probs_p(k - 1));
                    }
                    double covariance_h = frequency_n * probs_p(h - 1) * (static_cast<double>(h) - expected_X);
                    for(arma::uword j = 0; j < P; j++){
                        if(j != p && stats_n(j) != 0.0){
                            hessian_c(index_threshold_h, index_interaction(p, j)) += covariance_h * stats_n(j);
                            hessian_c(index_interaction(p, j), index_threshold_h) += covariance_h * stats_n(j);
                        }
                    }
                }
                for(arma::uword j = 0; j < P; j++){
                    if(j == p || stats_n(j) == 0.0){
                        continue;
                    }
                    for(arma::uword k = 0; k < P; k++){
                        if(k == p || stats_n(k) == 0.0){
                            continue;
                        }
                        hessian_c(index_interaction(p, j), index_interaction(p, k)) += frequency_n * stats_n(j) * stats_n(k) * variance_X;
                    }
                }
            }

            loglik_chunk(c) -= frequency_n * loglik_n;  // (-=) because negative log pseudo-likelihood
            gradient_patterns.col(n) = -gradient_n;     // (-) because negative gradient
            gradient_chunk.col(c) += frequency_n * gradient_patterns.col(n);
        }
    }

    // --- Sum over the chunks, in order ---
    double loglik = arma::accu(loglik_chunk);
    arma::vec gradient = arma::sum(gradient_chunk, 1);
    arma::mat hessian(n_pars, n_pars, arma::fill::zeros);
    for(arma::uword c = 0; c < n_chunks; c++){
        hessian += hessian_chunk.slice(c);
    }

    // --- Covariance of the score: sum_n frequency_n g_n g_n', in one matrix product outside the parallel loop ---
    gradient_patterns.each_row() %= arma::sqrt(frequency).t();
    arma::mat square_score = gradient_patterns * gradient_patterns.t();

    // --- Godambe-Huber-White (GHW) covariance ---
    arma::mat inverse_negative_hessian = arma::inv_sympd(hessian);
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
double cpp_npseudologlik(
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

    // loop over people, in parallel over ncores threads: each person writes only to its own entry of loglik_vec
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores) if(ncores > 1) private(p, h)
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


// Deduplicate rows of `data`, returning unique rows and their frequencies
void get_data_unique(arma::mat& unique_data, arma::vec& frequency, const arma::mat& data) {
    arma::uword N = data.n_rows;
    arma::uword P = data.n_cols;

    // Count each distinct row with a std::map keyed by the row itself (as a vector of integers;
    // the values are small non-negative integers 0..m-1), keeping the order of first appearance.
    std::map<std::vector<int>, arma::uword> row_counts;
    std::vector<std::vector<int>> row_order;   // preserves first-seen order

    for (arma::uword n = 0; n < N; n++) {
        std::vector<int> row(P);
        for (arma::uword p = 0; p < P; p++) row[p] = static_cast<int>(data(n, p));

        auto it = row_counts.find(row);
        if (it == row_counts.end()) {
            row_counts[row] = 1;
            row_order.push_back(row);
        } else {
            it->second++;
        }
    }

    arma::uword n_unique = row_order.size();
    unique_data.set_size(n_unique, P);
    frequency.set_size(n_unique);

    for (arma::uword i = 0; i < n_unique; i++) {
        for (arma::uword p = 0; p < P; p++) unique_data(i, p) = static_cast<double>(row_order[i][p]);
        frequency(i) = static_cast<double>(row_counts[row_order[i]]);
    }
}


// function to calculate the negative pseudologlikelihood (with or without prior) of a discrete MRF model at many parameter
// vectors at once (one per column of pars_draws), used by the Savage-Dickey BSIR step of dmrfit(). Same value as
// cpp_npseudologlik() at each column, computed with matrix operations over the unique response patterns of the data
// (weighted by their frequencies) and a stabilized log-sum-exp in the normalizing constants of the full conditionals.
// [[Rcpp::export]]
arma::vec cpp_npseudologlik_draws(
    const arma::mat &pars_draws, // matrix of size [n_pars x M], one parameter vector [thresholds, interactions] per column
    const arma::mat &data, // matrix of size [N x P] with the data (categories coded 0, ..., m-1), without cross-products
    const arma::uvec &n_categories,
    const bool &with_prior,
    const double &thresholds_alpha,
    const double &thresholds_beta,
    const double &interactions_location,
    const double &interactions_scale)
{
    arma::uword P = n_categories.n_elem;
    arma::uword n_thresholds = arma::accu(n_categories - 1);
    arma::uword n_pars = pars_draws.n_rows;
    arma::uword M = pars_draws.n_cols;

    // --- Unique response patterns and their frequencies ---
    arma::mat X;
    arma::vec frequency;
    get_data_unique(X, frequency, data);
    arma::uword U = X.n_rows;

    // --- Position of the first threshold of each variable in the parameter vector ---
    arma::uvec offset(P, arma::fill::zeros);
    for (arma::uword p = 1; p < P; p++) offset(p) = offset(p - 1) + n_categories(p - 1) - 1;

    // --- Count of each observed threshold over the unique patterns, weighted by their frequencies ---
    arma::vec threshold_counts(n_thresholds, arma::fill::zeros);
    for (arma::uword u = 0; u < U; u++) {
        for (arma::uword p = 0; p < P; p++) {
            arma::uword x = static_cast<arma::uword>(X(u, p));
            if (x > 0) threshold_counts(offset(p) + x - 1) += frequency(u);
        }
    }

    // lower-triangle indices (column-major, as in cpp_npseudologlik) of the interaction parameters
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P, P), -1);

    arma::vec out(M);
    for (arma::uword m = 0; m < M; m++) {
        Rcpp::checkUserInterrupt();
        arma::vec thresholds = pars_draws(arma::span(0, n_thresholds - 1), m);
        arma::vec interactions_vec = pars_draws(arma::span(n_thresholds, n_pars - 1), m);
        arma::mat lower(P, P, arma::fill::zeros);
        lower(lower_indices) = interactions_vec;
        arma::mat interactions = arma::symmatl(lower); // symmetric, zero diagonal

        // --- Rest scores sum_{j != p} x_j sigma_pj of every pattern and variable ---
        arma::mat rest = X * interactions; // [U x P]

        // --- Numerator: observed thresholds and interactions (x' Sigma x = sum_{i<j} 2 x_i x_j sigma_ij) ---
        double loglik = arma::dot(threshold_counts, thresholds) + arma::dot(frequency, arma::sum(X % rest, 1));

        // --- Normalizing constants of the full conditionals: log(1 + sum_h exp(mu_ph + h * rest_p)), stabilized ---
        for (arma::uword p = 0; p < P; p++) {
            arma::uword H = n_categories(p) - 1;
            arma::mat eta(U, H + 1, arma::fill::zeros); // column 0: baseline category (eta = 0)
            for (arma::uword h = 1; h <= H; h++) eta.col(h) = thresholds(offset(p) + h - 1) + static_cast<double>(h) * rest.col(p);
            arma::vec eta_max = arma::max(eta, 1);
            arma::vec log_norm = eta_max + arma::log(arma::sum(arma::exp(eta.each_col() - eta_max), 1));
            loglik -= arma::dot(frequency, log_norm);
        }

        // --- Prior ---
        if (with_prior) {
            loglik += arma::accu(log_beta_prime(thresholds, thresholds_alpha, thresholds_beta)) +
                      arma::accu(log_dcauchy(interactions_vec, interactions_location, interactions_scale));
        }
        out(m) = -loglik;
    }
    return out;
}
