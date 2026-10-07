#' hdi (internal)
#' @description Highest posterior density interval of a set of draws: the shortest interval that contains
#'   \code{ceiling(prob * n)} of the \code{n} sorted draws (as \code{coda::HPDinterval}). Among intervals of equal
#'   length, the first one is returned. Missing values are dropped.
#' @param y numeric vector, the draws of one parameter
#' @param prob probability mass of the interval (default 0.95)
#' @return a vector with the lower and upper bounds (NA when there are fewer than 3 draws)
#' @noRd
.hdi <- function(y, prob = 0.95) {
    y <- sort.int(as.numeric(y), method = "quick") # sorting also drops NAs
    n <- length(y)
    window <- ceiling(prob * n) - 1 # an interval y[k], ..., y[k + window] contains window + 1 draws
    if (n < 3 || window < 1 || window >= n) return(c(NA_real_, NA_real_))

    # --- shortest interval among the n - window candidates ---
    lower <- y[seq_len(n - window)]
    upper <- y[(window + 1):n]
    k <- which.min(upper - lower)
    return(c(lower[k], upper[k]))
}


#' posterior_mode (internal)
#' @description Marginal posterior mode of a set of draws: the maximum of their kernel density estimate
#'   (stats::density with its default bandwidth), as in Arena and Marsman (2026). Missing values are dropped.
#' @param y numeric vector, the draws of one parameter
#' @return the mode (NA when there are fewer than 2 distinct draws)
#' @noRd
.posterior_mode <- function(y) {
    y <- y[!is.na(y)]
    if (length(unique(y)) < 2) return(NA_real_)
    d <- stats::density(y)
    return(d$x[which.max(d$y)])
}
