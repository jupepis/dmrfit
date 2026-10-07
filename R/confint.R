#' @title Confidence and credible intervals of a fitted discrete Markov random field
#'
#' @description Intervals for the parameters of a \code{dmrfit} or \code{dmrfit_bayes} fit. For a \code{dmrfit} fit,
#' the adjusted likelihood-ratio intervals computed with \code{dmrfit(..., lrt_intervals = TRUE)} (\code{method =
#' "lrt"}, see \code{\link{dmrfit}}), or Wald intervals from the Godambe-Huber-White (GHW) standard errors (\code{method = "wald"}). For a
#' \code{dmrfit_bayes} fit, highest posterior density intervals of the posterior draws.
#'
#' @param object a \code{dmrfit} or \code{dmrfit_bayes} object.
#' @param parm the parameters, by name (as in \code{names(object$argument)}) or position. By default, all free
#'   parameters (for a constrained fit, the thresholds and the interactions present in the structure).
#' @param level confidence or credible level. By default, the level of a \code{dmrfit} fit (its \code{level}
#'   argument), and 0.95 for a \code{dmrfit_bayes} fit. For \code{method = "lrt"}, it must be the level of the fit.
#' @param method for a \code{dmrfit} fit: \code{"lrt"} (default) or \code{"wald"}. Ignored for a \code{dmrfit_bayes}
#'   fit.
#' @param ... further arguments (currently unused).
#'
#' @return a matrix with one row per parameter and the lower and upper bounds as columns.
#'
#' @seealso \code{\link{dmrfit}}, \code{\link{dmrfit_bayes}}
#'
#' @examples
#' data(rads2)
#' dysphoria <- names(which(attr(rads2, "clusters") == "Dysphoria"))
#' fit <- dmrfit(rads2[, dysphoria], lrt_intervals = c("mu[1,3]", "theta[2,1]"))
#' confint(fit, parm = c("mu[1,3]", "theta[2,1]"))
#' confint(fit, parm = c("mu[1,3]", "theta[2,1]"), method = "wald")
#'
#' @method confint dmrfit
#' @export
#'
confint.dmrfit <- function(object, parm = NULL, level = NULL, method = c("lrt", "wald"), ...) {

    if (!inherits(object, "dmrfit")) {
        stop("object is not of class 'dmrfit'")
    }
    if (is.null(level)) {
        level <- if (inherits(object, "dmrfit_bayes")) 0.95 else .fit_level(object)
    }
    if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
        stop("level must be a number between 0 and 1.")
    }
    par_names <- names(object$argument)

    # --- parameters: by default the free ones (the absent edges of a constrained fit have no interval) ---
    n_thresholds <- sum(object$n_categories - 1)
    free <- if (isTRUE(object$structured) && !is.null(object$structure)) {
        c(rep(TRUE, n_thresholds), object$structure[lower.tri(object$structure)] == 1)
    } else rep(TRUE, length(par_names))
    if (is.null(parm)) {
        parm <- par_names[free]
    } else if (is.numeric(parm)) {
        if (any(parm < 1 | parm > length(par_names) | parm != round(parm))) {
            stop("parm must be positions between 1 and ", length(par_names), ".")
        }
        parm <- par_names[parm]
    } else if (!all(parm %in% par_names)) {
        stop("unknown parameter(s) in parm: ", paste(setdiff(parm, par_names), collapse = ", "), ". See names(object$argument).")
    }
    if (any(!free[match(parm, par_names)])) {
        stop("parm contains absent edges of the constrained fit: ", paste(parm[!free[match(parm, par_names)]], collapse = ", "), ".")
    }

    # --- bounds ---
    if (inherits(object, "dmrfit_bayes")) {
        bounds <- t(apply(object$draws[match(parm, par_names), , drop = FALSE], 1, .hdi, prob = level))
    } else {
        method <- match.arg(method)
        if (method == "lrt") {
            lrt <- object$lrt_intervals
            if (is.null(lrt)) {
                stop("No likelihood-ratio intervals in the fit: refit with dmrfit(..., lrt_intervals = TRUE), or use method = \"wald\".")
            }
            if (!isTRUE(all.equal(attr(lrt, "level"), level))) {
                stop("The likelihood-ratio intervals were computed at level ", attr(lrt, "level"),
                     ": refit with dmrfit(..., level = ", level, ") for another level.")
            }
            bounds <- as.matrix(lrt[parm, c("lower", "upper")])
        } else {
            se <- sqrt(diag(object$utils$HW))[match(parm, par_names)]
            z <- stats::qnorm((1 + level) / 2)
            bounds <- cbind(object$argument[parm] - z * se, object$argument[parm] + z * se)
        }
    }

    a <- (1 - level) / 2
    dimnames(bounds) <- list(parm, paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%"))
    return(bounds)
}
