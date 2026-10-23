# Edge selection (paper code): the summary of a sampler run, and the metrics of a run against the reference posterior
# and the truth, for the simulation study (es_project/). The summary is computed from the full chain inside the job
# (all kept iterations). Only the reference keeps a thinned set of draws, for the overlapping indices of the methods.

# Overlapping indices

#' Grid probabilities of a kernel density estimate
#'
#' @param x draws
#' @param from,to grid range
#' @param n_grid number of grid points
#' @return the density on the grid, floored at the smallest positive double (as in the CoRe paper) and normalized to
#'   sum to 1
#' @noRd
.paper_kde_grid <- function(x, from, to, n_grid) {
    y <- stats::density(x, n = n_grid, from = from, to = to)$y
    y[y < .Machine$double.xmin] <- .Machine$double.xmin
    y / sum(y)
}

#' Grid probabilities of a two-dimensional kernel density estimate
#'
#' @param x draws [draws x 2]
#' @param lims grid range, c(x1 range, x2 range)
#' @param n_grid number of grid points per dimension
#' @return the density on the n_grid x n_grid grid, floored and normalized as in \code{.paper_kde_grid()}
#' @noRd
.paper_kde2_grid <- function(x, lims, n_grid) {
    z <- MASS::kde2d(x[, 1], x[, 2], n = n_grid, lims = lims)$z
    z[z < .Machine$double.xmin] <- .Machine$double.xmin
    z / sum(z)
}

#' Integrals of the minimum and the maximum of one continuous component of two mixtures
#'
#' @param w1,w2 masses of the component
#' @param g1,g2 grid probabilities of its densities, or NULL when a side has too few draws to estimate it (the minimum
#'   is then taken as 0 and the maximum as w1 + w2, exact in the limit of a vanishing mass on that side)
#' @return c(min, max)
#' @noRd
.paper_component_minmax <- function(w1, w2, g1, g2) {
    if (is.null(g1) || is.null(g2)) return(c(min = 0, max = w1 + w2))
    c(min = sum(pmin(w1 * g1, w2 * g2)), max = sum(pmax(w1 * g1, w2 * g2)))
}

#' Marginal overlapping index of a spike-and-slab posterior
#'
#' The index of the CoRe paper, the integral of min(f1, f2) over the integral of max(f1, f2), for a posterior that is a
#' point mass at zero (probability 1 - p) plus a continuous part (probability p, density g):
#' \code{(min(1 - p1, 1 - p2) + sum min(p1 g1, p2 g2)) / (max(1 - p1, 1 - p2) + sum max(p1 g1, p2 g2))}. For a
#' parameter that is never zero (a threshold) it is the index of the CoRe paper.
#'
#' @param d1,d2 draws of one parameter (absent edge = exactly 0). d2 is the reference
#' @param pip1,pip2 inclusion probabilities (default: the fraction of nonzero draws)
#' @param n_grid grid size of the kernel estimates
#' @param min_draws smallest number of nonzero draws for a kernel estimate
#' @return \code{overlap} (point mass and slab), \code{slab} (the index of the two slabs, given inclusion, NA if a side
#'   has fewer than \code{min_draws} nonzero draws) and \code{inclusion}, \code{(1 - |p1 - p2|) / (1 + |p1 - p2|)}, the
#'   overlap if the two slabs were identical
#' @noRd
paper_marginal_overlap <- function(d1, d2, pip1 = mean(d1 != 0), pip2 = mean(d2 != 0), n_grid = 1024,
                                   min_draws = 10) {
    nz1 <- d1[d1 != 0]
    nz2 <- d2[d2 != 0]
    g1 <- g2 <- NULL
    if (length(nz1) >= min_draws || length(nz2) >= min_draws) {
        bounds <- range(c(nz1, nz2))
        if (length(nz1) >= min_draws) g1 <- .paper_kde_grid(nz1, bounds[1], bounds[2], n_grid)
        if (length(nz2) >= min_draws) g2 <- .paper_kde_grid(nz2, bounds[1], bounds[2], n_grid)
    }
    slab <- .paper_component_minmax(pip1, pip2, g1, g2)
    lo <- min(1 - pip1, 1 - pip2) + slab[["min"]]
    hi <- max(1 - pip1, 1 - pip2) + slab[["max"]]
    dp <- abs(pip1 - pip2)
    c(overlap = lo / hi,
      slab = if (!is.null(g1) && !is.null(g2)) sum(pmin(g1, g2)) / sum(pmax(g1, g2)) else NA_real_,
      inclusion = (1 - dp) / (1 + dp))
}

#' Bivariate overlapping index of a spike-and-slab posterior
#'
#' The integrals of \code{paper_marginal_overlap()} summed over the four inclusion patterns of a pair: (0, 0) a point
#' mass, (1, 0) and (0, 1) one-dimensional densities on a line, (1, 1) a two-dimensional density.
#'
#' @param d1,d2 draws of a pair of parameters [draws x 2] (absent edge = exactly 0). d2 is the reference
#' @param n_grid grid size per dimension of the two-dimensional kernel estimates
#' @param min_draws smallest number of draws in an inclusion pattern for a kernel estimate
#' @return \code{overlap} (all four patterns), \code{both_in} (the index of the (1, 1) pattern alone, NA if a side has
#'   fewer than \code{min_draws} draws there) and \code{coinclusion_tv}, the total variation between the two 2 x 2
#'   tables of the inclusion patterns
#' @noRd
paper_bivariate_overlap <- function(d1, d2, n_grid = 128, min_draws = 10) {
    d1 <- as.matrix(d1)
    d2 <- as.matrix(d2)
    pattern <- function(d) (d[, 1] != 0) + 2 * (d[, 2] != 0) # 0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (1, 1)
    c1 <- pattern(d1)
    c2 <- pattern(d2)
    w1 <- tabulate(c1 + 1, 4) / nrow(d1)
    w2 <- tabulate(c2 + 1, 4) / nrow(d2)
    lo <- min(w1[1], w2[1])
    hi <- max(w1[1], w2[1])
    for (c in 1:2) {
        x1 <- d1[c1 == c, c]
        x2 <- d2[c2 == c, c]
        g1 <- g2 <- NULL
        if (length(x1) >= min_draws || length(x2) >= min_draws) {
            bounds <- range(c(x1, x2))
            if (length(x1) >= min_draws) g1 <- .paper_kde_grid(x1, bounds[1], bounds[2], n_grid)
            if (length(x2) >= min_draws) g2 <- .paper_kde_grid(x2, bounds[1], bounds[2], n_grid)
        }
        mm <- .paper_component_minmax(w1[c + 1], w2[c + 1], g1, g2)
        lo <- lo + mm[["min"]]
        hi <- hi + mm[["max"]]
    }
    x1 <- d1[c1 == 3, , drop = FALSE]
    x2 <- d2[c2 == 3, , drop = FALSE]
    g1 <- g2 <- NULL
    if (nrow(x1) >= min_draws || nrow(x2) >= min_draws) {
        lims <- c(range(c(x1[, 1], x2[, 1])), range(c(x1[, 2], x2[, 2])))
        if (nrow(x1) >= min_draws) g1 <- .paper_kde2_grid(x1, lims, n_grid)
        if (nrow(x2) >= min_draws) g2 <- .paper_kde2_grid(x2, lims, n_grid)
    }
    mm <- .paper_component_minmax(w1[4], w2[4], g1, g2)
    lo <- lo + mm[["min"]]
    hi <- hi + mm[["max"]]
    c(overlap = lo / hi,
      both_in = if (!is.null(g1) && !is.null(g2)) sum(pmin(g1, g2)) / sum(pmax(g1, g2)) else NA_real_,
      coinclusion_tv = sum(abs(w1 - w2)) / 2)
}

#' Pairs of parameters for the bivariate overlapping index
#'
#' As in the CoRe paper: two interactions that are present in the truth and share a node (any two present
#' interactions if none share a node, any two interactions sharing a node if fewer than two are present), a threshold
#' of one node of a present interaction and that interaction, and two thresholds of different nodes. The pairs are
#' drawn with \code{seed} and the caller's random-number state is restored, so every method of a data set gets the
#' same pairs.
#'
#' @param truth a truth (\code{edges}, \code{n_categories})
#' @param seed seed of the draw
#' @return a list with the three pairs, each a vector of two parameter indices (1-based, thresholds then interactions)
#' @noRd
paper_edge_pairs <- function(truth, seed) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
        on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
    P <- length(truth$n_categories)
    n_thr <- sum(truth$n_categories - 1)
    ij <- which(lower.tri(diag(P)), arr.ind = TRUE) # dmrfit's interaction order
    thr_node <- rep(seq_len(P), truth$n_categories - 1)
    pick <- function(x) x[sample.int(length(x), 1)]
    shares <- function(a, b) length(intersect(ij[a, ], ij[b, ])) > 0
    present <- which(truth$edges == 1)
    all_pairs <- t(utils::combn(seq_len(nrow(ij)), 2))
    cand <- if (length(present) >= 2) t(utils::combn(present, 2)) else all_pairs
    cand_shared <- cand[apply(cand, 1, function(x) shares(x[1], x[2])), , drop = FALSE]
    if (nrow(cand_shared) == 0 && length(present) < 2) {
        cand_shared <- all_pairs[apply(all_pairs, 1, function(x) shares(x[1], x[2])), , drop = FALSE]
    }
    if (nrow(cand_shared) > 0) cand <- cand_shared
    ii <- cand[sample.int(nrow(cand), 1), ]
    edge <- if (length(present) >= 1) pick(present) else sample.int(nrow(ij), 1)
    t1 <- sample.int(n_thr, 1)
    list(interaction_interaction = n_thr + ii,
         threshold_interaction = c(pick(which(thr_node %in% ij[edge, ])), n_thr + edge),
         threshold_threshold = c(t1, pick(which(thr_node != thr_node[t1]))))
}

# Summary of a run

#' Effective sample size of a chain
#'
#' @param x draws
#' @return \code{coda::effectiveSize(x)}, or NA for a constant chain
#' @noRd
.paper_ess <- function(x) {
    if (stats::var(x) == 0) return(NA_real_)
    unname(coda::effectiveSize(x))
}

#' Summary of an edge-selection sampler run, from all kept iterations
#'
#' Computed inside the job, before any draws are dropped. Per interaction: the Rao-Blackwellized and the raw inclusion
#' probability, the model-averaged posterior mean and SD, the posterior mean and SD given inclusion, the ESS of the
#' indicator and of the draws, the switch rate of the indicator and its efficiency (the switch rate over 2 p (1 - p),
#' the rate of independent draws with the same inclusion probability p), the Monte Carlo standard error of the
#' inclusion probability, \code{sqrt(p (1 - p) / ESS)} (0 for an indicator that never switches), and the move counts.
#' Per threshold: the posterior mean, SD and ESS. For the run: the time, the acceptance rates and the number of kept
#' iterations.
#'
#' @param fit output of \code{paper_edge_sampler()}
#' @param seconds_total run time including the ingredients (default: the sampler's time)
#' @param n_store number of draws to keep (evenly spaced over the kept iterations), 0 for none
#' @return a list with \code{edges}, \code{thresholds}, \code{run} and, if \code{n_store > 0}, \code{draws}
#' @noRd
paper_edge_summary <- function(fit, seconds_total = fit$seconds_elapsed, n_store = 0) {
    D <- fit$draws
    n_thr <- fit$n_thresholds
    K <- ncol(D) - n_thr
    th <- D[, n_thr + seq_len(K), drop = FALSE]
    g <- fit$gamma
    T_kept <- nrow(D)
    pip <- fit$inclusion_rb
    raw <- colMeans(g)
    cond <- function(f) vapply(seq_len(K), function(k) if (sum(g[, k]) > 1) f(th[g[, k], k]) else NA_real_, 0)
    ess_gamma <- apply(g, 2, function(x) .paper_ess(as.numeric(x)))
    switches <- colSums(abs(diff(g))) / (T_kept - 1)
    ec <- fit$edge_counts
    edges <- data.frame(edge = seq_len(K), pip = pip, pip_raw = raw,
                        mean = colMeans(th), sd = apply(th, 2, stats::sd),
                        cmean = cond(mean), csd = cond(stats::sd),
                        ess_gamma = ess_gamma, ess_theta = apply(th, 2, .paper_ess),
                        switch_rate = switches, switch_efficiency = switches / (2 * raw * (1 - raw)),
                        mc_se = ifelse(is.na(ess_gamma), 0, sqrt(pip * (1 - pip) / ess_gamma)),
                        add_try = ec[, "add_try"], add_acc = ec[, "add_acc"], del_try = ec[, "del_try"],
                        del_acc = ec[, "del_acc"], within_try = ec[, "within_try"], within_acc = ec[, "within_acc"])
    edges$switch_efficiency[!is.finite(edges$switch_efficiency)] <- NA
    thr <- D[, seq_len(n_thr), drop = FALSE]
    thresholds <- data.frame(threshold = seq_len(n_thr), mean = colMeans(thr), sd = apply(thr, 2, stats::sd),
                             ess = apply(thr, 2, .paper_ess))
    run <- data.frame(method = fit$method, seconds = seconds_total, seconds_sampler = fit$seconds_elapsed,
                      kept_iterations = T_kept, block_acceptance = fit$block_acceptance,
                      add_acceptance = sum(ec[, "add_acc"]) / max(1, sum(ec[, "add_try"])),
                      delete_acceptance = sum(ec[, "del_acc"]) / max(1, sum(ec[, "del_try"])),
                      within_acceptance = sum(ec[, "within_acc"]) / max(1, sum(ec[, "within_try"])))
    out <- list(edges = edges, thresholds = thresholds, run = run)
    if (n_store > 0) out$draws <- D[unique(round(seq(1, T_kept, length.out = min(n_store, T_kept)))), , drop = FALSE]
    out
}

# Metrics
#   .compute_metrics() runs, in order, the metrics of the four groups decided for the study:
#   A. inclusion against the reference       .metrics_inclusion()
#   B. continuous part against the reference .metrics_continuous()
#   C. graph recovery against the truth      .metrics_graph()
#   E. efficiency                            .metrics_efficiency()
# (D, estimation against the truth, was left out.) The breakdowns by |theta| and by the reference PIP are done in the
# analysis, from the per-interaction table.

#' A. Inclusion against the reference posterior
#'
#' @param pip,pip_ref inclusion probabilities of the run and of the reference
#' @param present true indicators (logical)
#' @param uncertain logical, the uncertain interactions (reference PIP in the range of \code{.compute_metrics()})
#' @param clamp clamping of the PIPs for the log-odds
#' @return a list with the per-interaction table \code{edges} (absolute and squared PIP error, log-odds difference,
#'   agreement of the median probability graphs) and the one-row table \code{data_set} (A1 PIP mean absolute error,
#'   A2 squared error, A3 log-odds bias on null and on present edges, A4 agreement, A6 mean absolute error on the
#'   uncertain edges)
#' @noRd
.metrics_inclusion <- function(pip, pip_ref, present, uncertain, clamp) {
    lo <- function(p) stats::qlogis(pmin(pmax(p, clamp), 1 - clamp))
    edges <- data.frame(abs_error = abs(pip - pip_ref), sq_error = (pip - pip_ref)^2,
                        logodds_diff = lo(pip) - lo(pip_ref), agree = (pip > 0.5) == (pip_ref > 0.5))
    data_set <- data.frame(pip_mae = mean(edges$abs_error), pip_mse = mean(edges$sq_error),
                           logodds_bias_null = .mean_finite(edges$logodds_diff[!present]),
                           logodds_bias_present = .mean_finite(edges$logodds_diff[present]),
                           agreement = mean(edges$agree),
                           pip_mae_uncertain = .mean_finite(edges$abs_error[uncertain]), n_uncertain = sum(uncertain))
    list(edges = edges, data_set = data_set)
}

#' B. Continuous part against the reference posterior
#'
#' @param e,r per-interaction tables of the run and of the reference (\code{paper_edge_summary()})
#' @param ts,tr per-threshold tables of the run and of the reference
#' @param draws,reference_draws draws of the run (all kept iterations) and stored draws of the reference, or NULL
#' @param pairs output of \code{paper_edge_pairs()}, or NULL
#' @param n_thr number of thresholds
#' @param min_pip smallest PIP (run and reference) for the moments given inclusion
#' @return a list with \code{edges} (SD ratio and standardized mean difference given inclusion, marginal overlap of
#'   point mass and slab and of the slab alone), \code{thresholds} (standardized mean difference, SD ratio, overlap),
#'   \code{pairs} (bivariate overlap, of both interactions included alone, co-inclusion total variation) and the
#'   one-row \code{data_set} (B1 median SD ratio, B3 standardized mean difference, B4 thresholds, B5 overlap with B6
#'   the slab, B7 threshold overlap, B8 bivariate overlap with the both-included part)
#' @noRd
.metrics_continuous <- function(e, r, ts, tr, draws, reference_draws, pairs, n_thr, min_pip) {
    K <- nrow(e)
    with_draws <- !is.null(draws) && !is.null(reference_draws)
    cond_ok <- e$pip > min_pip & r$pip > min_pip & is.finite(e$csd) & is.finite(r$csd) & r$csd > 0
    edges <- data.frame(sd_ratio = ifelse(cond_ok, e$csd / r$csd, NA_real_),
                        mean_std = ifelse(cond_ok, (e$cmean - r$cmean) / r$csd, NA_real_),
                        overlap = NA_real_, overlap_slab = NA_real_)
    if (with_draws) {
        for (k in seq_len(K)) {
            ov <- paper_marginal_overlap(draws[, n_thr + k], reference_draws[, n_thr + k], e$pip[k], r$pip[k])
            edges$overlap[k] <- ov[["overlap"]]
            edges$overlap_slab[k] <- ov[["slab"]]
        }
    }
    thresholds <- data.frame(mean_std = (ts$mean - tr$mean) / tr$sd, sd_ratio = ts$sd / tr$sd, overlap = NA_real_)
    if (with_draws) {
        thresholds$overlap <- vapply(seq_len(n_thr), function(l)
            paper_marginal_overlap(draws[, l], reference_draws[, l])[["overlap"]], 0)
    }
    pair_tab <- NULL
    if (with_draws && !is.null(pairs)) {
        pair_tab <- do.call(rbind, lapply(names(pairs), function(nm) {
            ov <- paper_bivariate_overlap(draws[, pairs[[nm]]], reference_draws[, pairs[[nm]]])
            data.frame(pair = nm, first = pairs[[nm]][1], second = pairs[[nm]][2], overlap = ov[["overlap"]],
                       both_in = ov[["both_in"]], coinclusion_tv = ov[["coinclusion_tv"]])
        }))
    }
    pair_value <- function(nm, v) if (is.null(pair_tab)) NA_real_ else pair_tab[[v]][pair_tab$pair == nm]
    data_set <- data.frame(sd_ratio = .median_finite(edges$sd_ratio), mean_std = .mean_finite(edges$mean_std),
                           mean_std_abs = .mean_finite(abs(edges$mean_std)),
                           threshold_mean_std_abs = .mean_finite(abs(thresholds$mean_std)),
                           threshold_sd_ratio = .median_finite(thresholds$sd_ratio),
                           overlap = .mean_finite(edges$overlap), overlap_slab = .mean_finite(edges$overlap_slab),
                           threshold_overlap = .mean_finite(thresholds$overlap),
                           pair_ii = pair_value("interaction_interaction", "overlap"),
                           pair_ii_both_in = pair_value("interaction_interaction", "both_in"),
                           pair_ti = pair_value("threshold_interaction", "overlap"),
                           pair_ti_both_in = pair_value("threshold_interaction", "both_in"),
                           pair_tt = pair_value("threshold_threshold", "overlap"))
    list(edges = edges, thresholds = thresholds, pairs = pair_tab, data_set = data_set)
}

#' Graph recovery of one set of inclusion probabilities against the truth
#'
#' @param pip inclusion probabilities of the interactions
#' @param present true indicators (logical)
#' @return the AUC (threshold free: the area under the ROC curve over all PIP thresholds, from the ranks of the PIPs),
#'   the Brier score, and the sensitivity, specificity, precision and structural Hamming distance (false positives
#'   plus false negatives) of the median probability graph (PIP > 0.5)
#' @noRd
.graph_recovery <- function(pip, present) {
    selected <- pip > 0.5
    auc <- if (any(present) && any(!present)) {
        rk <- rank(pip)
        n1 <- sum(present)
        (sum(rk[present]) - n1 * (n1 + 1) / 2) / (n1 * sum(!present))
    } else NA_real_
    data.frame(auc = auc,
               sensitivity = if (any(present)) mean(selected[present]) else NA_real_,
               specificity = if (any(!present)) mean(!selected[!present]) else NA_real_,
               precision = if (any(selected)) mean(present[selected]) else NA_real_,
               brier = mean((pip - present)^2),
               shd = sum(selected != present))
}

#' C. Graph recovery against the truth, for the run and for the reference
#'
#' @param pip,pip_ref inclusion probabilities of the run and of the reference
#' @param present true indicators (logical)
#' @return a one-row table: \code{.graph_recovery()} of the run, and of the reference with the suffix \code{_ref}
#' @noRd
.metrics_graph <- function(pip, pip_ref, present) {
    run <- .graph_recovery(pip, present)
    ref <- .graph_recovery(pip_ref, present)
    names(ref) <- paste0(names(ref), "_ref")
    cbind(run, ref)
}

#' E. Efficiency
#'
#' @param e,r per-interaction tables of the run and of the reference
#' @param ts per-threshold table of the run
#' @param run run table of the run (\code{paper_edge_summary()})
#' @param uncertain logical, the uncertain interactions
#' @return a list with \code{edges} (ESS of the indicator and of the draws, switch efficiency, Monte Carlo floor of
#'   the PIP error, Monte Carlo standard error of the reference PIP) and the one-row \code{data_set} (E1 run time, E2
#'   ESS per second of the indicators of the uncertain edges, of the interactions and of the thresholds, E3 median
#'   switch efficiency of the uncertain edges, E4 acceptance rates, E5 Monte Carlo floor)
#' @noRd
.metrics_efficiency <- function(e, r, ts, run, uncertain) {
    edges <- data.frame(ess_gamma = e$ess_gamma, ess_theta = e$ess_theta, switch_efficiency = e$switch_efficiency,
                        mc_floor = sqrt(e$mc_se^2 + r$mc_se^2), mc_se_ref = r$mc_se)
    data_set <- data.frame(seconds = run$seconds,
                           ess_gamma_per_s = .mean_finite(edges$ess_gamma[uncertain]) / run$seconds,
                           ess_gamma_min_per_s = .min_finite(edges$ess_gamma[uncertain]) / run$seconds,
                           ess_theta_per_s = .mean_finite(edges$ess_theta) / run$seconds,
                           ess_threshold_per_s = .mean_finite(ts$ess) / run$seconds,
                           ess_threshold_min_per_s = .min_finite(ts$ess) / run$seconds,
                           switch_efficiency = .median_finite(edges$switch_efficiency[uncertain]),
                           block_acceptance = run$block_acceptance, add_acceptance = run$add_acceptance,
                           delete_acceptance = run$delete_acceptance, within_acceptance = run$within_acceptance,
                           mc_floor = .mean_finite(edges$mc_floor),
                           mc_floor_uncertain = .mean_finite(edges$mc_floor[uncertain]))
    list(edges = edges, data_set = data_set)
}

# mean, median and minimum of the finite values (NA if there are none)
.mean_finite <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
.median_finite <- function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
.min_finite <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_

#' Metrics of an edge-selection run against the reference posterior and the truth
#'
#' Runs, in order, the four groups of metrics: A inclusion and B continuous part against the reference, C graph
#' recovery against the truth, E efficiency (see their help).
#'
#' @param summary output of \code{paper_edge_summary()} for the run
#' @param reference output of \code{paper_edge_summary()} for the reference posterior
#' @param truth the truth of the data set (\code{true_pars}, \code{edges}, \code{n_categories})
#' @param draws draws of the run (all kept iterations), or NULL for no overlaps
#' @param reference_draws stored draws of the reference, or NULL
#' @param pairs output of \code{paper_edge_pairs()}, or NULL for no bivariate overlap
#' @param uncertain range of the reference PIP of the uncertain interactions
#' @param clamp clamping of the PIPs for the log-odds
#' @param min_pip smallest PIP (run and reference) for the moments given inclusion
#' @return a list with the per-interaction table \code{edges}, the per-threshold table \code{thresholds}, the per-pair
#'   table \code{pairs} and the one-row table \code{data_set} of the aggregated metrics of the four groups
#' @noRd
.compute_metrics <- function(summary, reference, truth, draws = NULL, reference_draws = NULL, pairs = NULL,
                             uncertain = c(0.05, 0.95), clamp = 0.005, min_pip = 0.1) {
    n_thr <- sum(truth$n_categories - 1)
    K <- length(truth$edges)
    P <- length(truth$n_categories)
    ij <- which(lower.tri(diag(P)), arr.ind = TRUE) # dmrfit's interaction order
    e <- summary$edges
    r <- reference$edges
    present <- truth$edges == 1
    unc <- r$pip > uncertain[1] & r$pip < uncertain[2]
    base <- data.frame(edge = seq_len(K), i = ij[, 1], j = ij[, 2], true_theta = truth$true_pars[n_thr + seq_len(K)],
                       present = present, pip_ref = r$pip, pip = e$pip, uncertain = unc)

    A <- .metrics_inclusion(e$pip, r$pip, present, unc, clamp)
    B <- .metrics_continuous(e, r, summary$thresholds, reference$thresholds, draws, reference_draws, pairs, n_thr,
                             min_pip)
    C <- .metrics_graph(e$pip, r$pip, present)
    E <- .metrics_efficiency(e, r, summary$thresholds, summary$run, unc)

    thresholds <- cbind(data.frame(threshold = seq_len(n_thr), true_mu = truth$true_pars[seq_len(n_thr)],
                                   ess = summary$thresholds$ess), B$thresholds)
    list(edges = cbind(base, A$edges, B$edges, E$edges), thresholds = thresholds, pairs = B$pairs,
         data_set = cbind(A$data_set, B$data_set, C, E$data_set))
}
