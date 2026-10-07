#' @title Plot a fitted discrete Markov random field
#'
#' @description Plots of a \code{dmrfit} or \code{dmrfit_bayes} fit, made with \pkg{ggplot2}. The interactions are
#' the estimates of \code{dmrfit()}, or a posterior summary of the draws of \code{dmrfit_bayes()} (see
#' \code{estimate}).
#' \describe{
#'   \item{\code{type = "network"}}{The nodes are the variables and the edges the pairwise interactions: the edge
#'     width is the absolute value of the interaction, the color and line type its sign. When the fit has
#'     Savage-Dickey Bayes factors, only the included edges are drawn (see Details).}
#'   \item{\code{type = "bf"}}{One circle per pair of variables (lower triangle): its color is the evidence for the
#'     edge given by its Savage-Dickey Bayes factor (see Details) and its area the absolute value of the interaction.
#'     Requires a fit with Bayes factors.}
#'   \item{\code{type = "trace"}}{The posterior draws (after burn-in) of the parameters in \code{pars} against the
#'     iteration, one panel per parameter. Only for a \code{dmrfit_bayes} fit.}
#'   \item{\code{type = "density"}}{The kernel density of the posterior draws of the parameters in \code{pars}, one
#'     panel per parameter, with a solid line at the posterior summary chosen by \code{estimate} and dashed lines at
#'     the highest posterior density interval of mass \code{prob}. Only for a \code{dmrfit_bayes} fit.}
#'   \item{\code{type = "intervals"}}{The estimate and the interval of every parameter in \code{pars}, the thresholds
#'     and the interactions in separate panels, each sorted by estimate. The intervals are the likelihood-ratio
#'     intervals of a \code{dmrfit} fit computed with \code{lrt_intervals = TRUE}, otherwise its Wald intervals, both
#'     at the level of the fit, or the highest posterior density intervals of mass \code{prob} of a
#'     \code{dmrfit_bayes} fit (see \code{\link{confint.dmrfit}}). The caption states which.}
#'   \item{\code{type = "centrality"}}{The expected influence of every variable, \eqn{EI_i = \sum_j \theta_{ij}}, the sum
#'     of its interactions, sorted by its estimate. For a \code{dmrfit_bayes} fit, the posterior summary chosen by
#'     \code{estimate} and the highest posterior density interval of mass \code{prob} of the expected influence
#'     computed on every draw, and next to it the posterior probability that the variable is the most central (the
#'     share of draws in which its expected influence is the largest). For a \code{dmrfit} fit, the estimate and its
#'     Wald interval at the level of the fit; the standard error of a sum of interactions follows from their sandwich
#'     covariance.}
#' }
#'
#' @param x a \code{dmrfit} or \code{dmrfit_bayes} object.
#' @param type the plot: \code{"network"} (default), \code{"bf"}, \code{"trace"}, \code{"density"},
#'   \code{"intervals"} or \code{"centrality"}.
#' @param estimate the posterior summary of the interactions of a \code{dmrfit_bayes} fit: \code{"mode"} (default) the
#'   marginal posterior mode, the maximum of the kernel density estimate of the draws (\code{stats::density}),
#'   \code{"mean"} or \code{"median"}. The mode is the default because the posterior can be asymmetric in small
#'   samples, which moves the mean and the median away from the mode; Arena and Marsman (2026) found this for the
#'   thresholds of rarely chosen response categories. Ignored for a \code{dmrfit} fit, which has point estimates.
#' @param pars (trace, density, intervals) the parameters to show (at most 9 for trace and density): their names, as in \code{names(x$argument)} (for
#'   instance \code{"sigma[2,1]"} for the interaction between variables 2 and 1, \code{"mu[1,2]"} for the second
#'   threshold of variable 1), or their positions in that vector. By default, the four interactions with the largest
#'   absolute posterior summary (see \code{estimate}) for trace and density, and all the free parameters for
#'   intervals. The parameters are labeled with the variable names.
#' @param prob (density, intervals, centrality) mass of the highest posterior density interval of a \code{dmrfit_bayes} fit
#'   (default 0.95). The intervals of a \code{dmrfit} fit are at its \code{level}.
#' @param all_edges (network) logical, whether to draw every estimated interaction instead of only the included ones
#'   (default FALSE). Every estimated interaction is also drawn for a fit without Bayes factors (only the free ones
#'   for a constrained fit).
#' @param groups (network) optional grouping of the nodes (for instance \code{attr(rads2, "clusters")}), a vector with
#'   one entry per variable, in the column order of the data or named by variable. At most four groups, shown by the
#'   color and shape of the nodes.
#' @param layout (network) the node positions: \code{"fr"} (default) the Fruchterman-Reingold layout of \pkg{igraph},
#'   which must be installed; \code{"circle"} the nodes on a circle, in the column order of the data; or a matrix of
#'   coordinates with two columns and one row per variable (in the column order of the data, or with the variable
#'   names as row names), rescaled by the same factor on both axes. The \code{"fr"} layout is computed from all
#'   estimated interactions, weighted by their absolute value, so that the nodes keep their positions when
#'   \code{all_edges} changes.
#' @param seed (network) random seed of the \code{"fr"} layout (default 30); other seeds give other arrangements. The
#'   caller's random number stream is restored on exit.
#' @param node_size (network) size of the nodes (default 10); the labels scale with it.
#' @param show_comments logical, whether to show the explanation below the plot (default TRUE): the rule for the drawn
#'   edges, or what the points, lines and probabilities show. With FALSE, only the plot and its legends are drawn.
#' @param ... further arguments (currently unused).
#'
#' @details The evidence for an edge is classified by its Savage-Dickey Bayes factor \eqn{BF_{01}}, the evidence for
#' the absence of the edge: included (\eqn{BF_{01} < 1/10}), weak included (\eqn{1/10 \le BF_{01} < 1/3}),
#' inconclusive (\eqn{1/3 \le BF_{01} < 3}), weak excluded (\eqn{3 \le BF_{01} < 10}) and excluded
#' (\eqn{BF_{01} \ge 10}), shown from orange (included) through gray to violet (excluded). The network plot draws the
#' included edges.
#'
#' The edge widths of the network plot are scaled to the range of the absolute interactions drawn in each plot. To
#' compare widths between plots, fix the range, for instance with
#' \code{+ ggplot2::scale_linewidth_continuous(limits = c(0, 1), range = c(0.3, 2.4))}.
#'
#' @return a \code{ggplot} object, which can be modified further with \pkg{ggplot2}.
#'
#' @references Arena, G. and Marsman, M. (2026). Bayesian inference for discrete Markov random fields through
#' coordinate rescaling. arXiv preprint. \doi{10.48550/arXiv.2601.17205}
#'
#' @seealso \code{\link{dmrfit}}, \code{\link{dmrfit_bayes}}
#'
#' @examples
#' data(rads2)
#' dysphoria <- names(which(attr(rads2, "clusters") == "Dysphoria"))
#' fit <- dmrfit(rads2[, dysphoria], with_prior = TRUE, savage_dickey = TRUE)
#' plot(fit)
#'
#' # all estimated interactions, with the nodes grouped and placed on a circle
#' groups <- c(rep("first", 3), rep("second", 4))
#' plot(fit, all_edges = TRUE, groups = groups, layout = "circle")
#'
#' # evidence for every edge from its Savage-Dickey Bayes factor
#' plot(fit, type = "bf")
#'
#' # trace and density of the posterior draws
#' fit_bayes <- dmrfit_bayes(rads2[, dysphoria], nsim = 1000, burnin = 500, progress = FALSE)
#' plot(fit_bayes, type = "trace")
#' plot(fit_bayes, type = "density", pars = c("sigma[2,1]", "mu[1,1]"))
#'
#' # estimates and intervals: Wald, likelihood-ratio and highest posterior density
#' plot(fit, type = "intervals")
#' fit_lrt <- dmrfit(rads2[, dysphoria], lrt_intervals = c("mu[1,3]", "mu[2,3]", "sigma[2,1]"))
#' plot(fit_lrt, type = "intervals", pars = c("mu[1,3]", "mu[2,3]", "sigma[2,1]"))
#' plot(fit_bayes, type = "intervals")
#'
#' # expected influence of every variable, without the explanation below the plot
#' plot(fit, type = "centrality", show_comments = FALSE)
#' plot(fit_bayes, type = "centrality")
#'
#' @method plot dmrfit
#' @export
#'
plot.dmrfit <- function(x, type = c("network", "bf", "trace", "density", "intervals", "centrality"), estimate = c("mode", "mean", "median"),
                        pars = NULL, prob = 0.95, all_edges = FALSE, groups = NULL, layout = c("fr", "circle"),
                        seed = 30, node_size = 10, show_comments = TRUE, ...) {

    type <- match.arg(type)
    if (!inherits(x, "dmrfit")) {
        stop("x is not of class 'dmrfit'")
    }
    if (!missing(estimate) && !inherits(x, "dmrfit_bayes")) {
        warning("estimate is ignored for a dmrfit fit, which has point estimates.")
    }
    estimate <- match.arg(estimate)
    if (type %in% c("trace", "density") && !inherits(x, "dmrfit_bayes")) {
        stop("type = \"", type, "\" needs posterior draws: fit the model with dmrfit_bayes().")
    }
    if (!is.numeric(prob) || length(prob) != 1 || prob <= 0 || prob >= 1) {
        stop("prob must be a number between 0 and 1.")
    }
    if (!is.logical(show_comments) || length(show_comments) != 1 || is.na(show_comments)) {
        stop("show_comments must be TRUE or FALSE.")
    }
    if (!is.logical(all_edges) || length(all_edges) != 1 || is.na(all_edges)) {
        stop("all_edges must be TRUE or FALSE.")
    }
    if (is.character(layout)) {
        layout <- match.arg(layout)
    }
    if (!is.numeric(seed) || length(seed) != 1) {
        stop("seed must be a single number.")
    }
    if (!is.numeric(node_size) || length(node_size) != 1 || node_size <= 0) {
        stop("node_size must be a positive number.")
    }

    p <- switch(type,
                network = .plot_network(x, estimate = estimate, all_edges = all_edges, groups = groups, layout = layout,
                                        seed = seed, node_size = node_size),
                bf = .plot_bf(x, estimate = estimate),
                trace = .plot_draws(x, type = "trace", estimate = estimate, pars = pars, prob = prob),
                density = .plot_draws(x, type = "density", estimate = estimate, pars = pars, prob = prob),
                intervals = .plot_intervals(x, estimate = estimate, pars = pars, prob = prob),
                centrality = .plot_centrality(x, estimate = estimate, prob = prob))

    # --- without the explanation below the plot (for instance, when it goes into the caption of a figure) ---
    if (!show_comments) {
        p <- p + ggplot2::labs(caption = NULL)
    }
    return(p)
}


# --- Colors (validated with the dataviz palette checker on the light surface) ---
# sign of the edges: blue (positive) and red (negative), with the line type as secondary encoding
.PLOT_SIGN_COLORS <- c(positive = "#2a78d6", negative = "#e34948")
# node groups: four hues that pass the all-pairs checks, with the node shape as secondary encoding
.PLOT_GROUP_COLORS <- c("#1baf7a", "#eda100", "#4a3aa7", "#e87ba4")
.PLOT_GROUP_SHAPES <- c(21, 22, 24, 23)
# --- Evidence for an edge from its Savage-Dickey Bayes factor BF_01 (evidence for exclusion) ---
# upper bounds of the BF_01 intervals; each interval runs from the previous bound (0 for the first) up to, but
# excluding, its own
.PLOT_BF_INTERVALS <- c("Included" = 1/10, "Weak included" = 1/3, "Inconclusive" = 3, "Weak excluded" = 10,
                        "Excluded" = Inf)
# from evidence for the presence of an edge (orange) through gray to its absence (violet); the weak classes are the
# two colors mixed half with white
.PLOT_EVIDENCE_COLORS <- c("Included" = "#eb6834", "Weak included" = "#f5b39a", "Inconclusive" = "gray65",
                           "Weak excluded" = "#a59dd3", "Excluded" = "#4a3aa7")
.PLOT_EVIDENCE_LABELS <- c("Included" = "Included (< 1/10)", "Weak included" = "Weak included (1/10 - 1/3)",
                           "Inconclusive" = "Inconclusive (1/3 - 3)", "Weak excluded" = "Weak excluded (3 - 10)",
                           "Excluded" = "Excluded (>= 10)")
.PLOT_SURFACE <- "#fcfcfb"
# trace and density plots: neutral colors (one series, no identity color)
.PLOT_TRACE <- "gray35"
.PLOT_DENSITY <- "#3d3d3a"
.PLOT_DENSITY_SHADE <- "#ebeae6"
.PLOT_INK <- "#1a1a19"


#' tint (internal)
#' @description Mixes a color with white.
#' @param col a color
#' @param amount share of white, between 0 (the color) and 1 (white)
#' @return the mixed color as a hex string
#' @noRd
.tint <- function(col, amount) {
    rgb_col <- grDevices::col2rgb(col) / 255
    mixed <- (1 - amount) * rgb_col + amount
    return(grDevices::rgb(mixed[1], mixed[2], mixed[3]))
}


#' interaction_estimates (internal)
#' @description Estimates of the pairwise interactions of a fit, one row per pair (i, j) with i > j: the estimates of
#'   dmrfit(), or a posterior summary of the draws of dmrfit_bayes().
#' @param x a dmrfit or dmrfit_bayes object
#' @param estimate the posterior summary for a dmrfit_bayes fit: "mode", "mean" or "median"
#' @return a data frame with the parameter name, the node indices i and j, and the estimate
#' @noRd
.interaction_estimates <- function(x, estimate = "mode") {
    n_thresholds <- sum(x$n_categories - 1)
    inter_idx <- (n_thresholds + 1):length(x$argument)
    par_names <- names(x$argument)[inter_idx]
    if (inherits(x, "dmrfit_bayes")) {
        summary_fun <- switch(estimate, mode = .posterior_mode, mean = mean, median = stats::median)
        estimate <- apply(x$draws[inter_idx, , drop = FALSE], 1, summary_fun)
    } else {
        estimate <- x$argument[inter_idx]
    }
    ij <- do.call(rbind, regmatches(par_names, gregexpr("[0-9]+", par_names)))
    return(data.frame(name = par_names, i = as.integer(ij[, 1]), j = as.integer(ij[, 2]), estimate = unname(estimate),
                      stringsAsFactors = FALSE))
}


#' network_layout (internal)
#' @description Node coordinates for the network plot: the Fruchterman-Reingold layout of igraph, computed from all
#'   estimated interactions (weighted by their absolute value, with a fixed seed and the caller's random number stream
#'   restored), or a circle, each axis rescaled to the interval from -1 to 1; or the coordinates given by the user, centered and rescaled
#'   by the same factor on both axes.
#' @param layout "fr", "circle", or a matrix with two columns and one row per variable
#' @param var_names the variable names
#' @param inter data frame of the estimated interactions (columns i, j and estimate)
#' @param seed (network) random seed of the "fr" layout
#' @return a P x 2 matrix with the variable names as row names
#' @noRd
.network_layout <- function(layout, var_names, inter, seed) {
    P <- length(var_names)
    rescale_axis <- function(z) if (diff(range(z)) > 0) 2 * (z - min(z)) / diff(range(z)) - 1 else 0 * z

    if (!is.character(layout)) {
        # --- coordinates given by the user ---
        layout <- as.matrix(layout)
        if (ncol(layout) != 2 || nrow(layout) != P || !is.numeric(layout)) {
            stop("layout must be \"fr\", \"circle\", or a numeric matrix with two columns and one row per variable (", P, ").")
        }
        if (!is.null(rownames(layout))) {
            if (!setequal(rownames(layout), var_names)) stop("the row names of layout must be the variable names.")
            layout <- layout[var_names, , drop = FALSE]
        }
        layout <- sweep(layout, 2, (apply(layout, 2, max) + apply(layout, 2, min)) / 2) # centered
        layout <- layout / max(abs(layout), .Machine$double.eps)                         # same factor on both axes
    } else if (layout == "fr") {
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("layout = \"fr\" requires the igraph package: install it with install.packages(\"igraph\"), or use ",
                 "layout = \"circle\".")
        }
        inter <- inter[which(inter$estimate != 0), , drop = FALSE] # absent edges of a constrained fit (and NA modes)
        g <- igraph::make_empty_graph(n = P, directed = FALSE)
        g <- igraph::add_edges(g, as.vector(t(as.matrix(inter[, c("i", "j")]))), weight = abs(inter$estimate))
        # fixed seed for a reproducible layout; the caller's random number stream is restored
        old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
        on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()))
        set.seed(seed)
        layout <- apply(igraph::layout_with_fr(g, weights = igraph::E(g)$weight), 2, rescale_axis)
    } else {
        angle <- pi / 2 - 2 * pi * (seq_len(P) - 1) / P # clockwise from the top
        layout <- cbind(cos(angle), sin(angle))
    }

    dimnames(layout) <- list(var_names, c("x", "y"))
    return(layout)
}


#' plot_network (internal)
#' @description Network plot of plot.dmrfit(type = "network"); see its documentation for the arguments.
#' @return a ggplot object
#' @noRd
.plot_network <- function(x, estimate, all_edges, groups, layout, seed, node_size) {
    P <- x$P
    var_names <- if (is.null(x$var_names)) paste0("V", seq_len(P)) else x$var_names
    inter <- .interaction_estimates(x, estimate)

    # --- edges to draw ---
    sd <- x$savage_dickey
    if (!all_edges && !is.null(sd)) {
        drawn <- names(sd$bf_01)[.evidence_class(sd$bf_01) == "Included"]
        rule <- "Edges: included interactions (BF01 < 1/10)"
    } else {
        if (!all_edges) {
            message("No Savage-Dickey Bayes factors in the fit (savage_dickey = FALSE): every estimated interaction is drawn.")
        }
        drawn <- inter$name
        if (isTRUE(x$structured) && !is.null(x$structure)) {
            drawn <- inter$name[x$structure[cbind(inter$i, inter$j)] == 1]
        }
        rule <- "Edges: all estimated interactions"
    }
    edges <- inter[which(inter$name %in% drawn & inter$estimate != 0), , drop = FALSE]

    # --- nodes ---
    xy <- .network_layout(layout, var_names, inter, seed)
    nodes <- data.frame(name = var_names, x = xy[, "x"], y = xy[, "y"], stringsAsFactors = FALSE)
    if (!is.null(groups)) {
        if (!is.null(names(groups))) {
            if (!all(var_names %in% names(groups))) stop("groups must have an entry for every variable.")
            groups <- groups[var_names]
        }
        if (length(groups) != P) stop("groups must have one entry per variable (", P, ").")
        nodes$group <- if (is.factor(groups)) droplevels(groups) else factor(groups, levels = unique(groups))
        if (nlevels(nodes$group) > length(.PLOT_GROUP_COLORS)) {
            stop("At most ", length(.PLOT_GROUP_COLORS), " groups can be shown.")
        }
    }

    # --- edges as curves between the node coordinates ---
    edges$sign <- factor(ifelse(edges$estimate > 0, "positive", "negative"), levels = c("positive", "negative"))
    edges$weight <- abs(edges$estimate)
    edges$x <- xy[edges$i, "x"]
    edges$y <- xy[edges$i, "y"]
    edges$xend <- xy[edges$j, "x"]
    edges$yend <- xy[edges$j, "y"]
    edges <- edges[order(edges$weight), , drop = FALSE] # strongest edges on top

    p <- ggplot2::ggplot() +
        ggplot2::geom_curve(data = edges, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                                                       colour = .data$sign, linetype = .data$sign, linewidth = .data$weight),
                            curvature = 0.15, lineend = "round", alpha = 0.85) +
        ggplot2::scale_colour_manual(values = .PLOT_SIGN_COLORS, name = NULL) +
        ggplot2::scale_linetype_manual(values = c(positive = "solid", negative = "22"), name = NULL) +
        ggplot2::scale_linewidth_continuous(range = c(0.3, 2.4), name = "|interaction|")

    if (is.null(groups)) {
        p <- p + ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = node_size,
                                     fill = .PLOT_SURFACE, colour = "gray45", stroke = 1.2)
    } else {
        # group color on the node outline, a light tint of it as fill (the labels stay dark on a light fill)
        group_cols <- stats::setNames(.PLOT_GROUP_COLORS[seq_len(nlevels(nodes$group))], levels(nodes$group))
        group_tints <- vapply(group_cols, .tint, character(1), amount = 0.75)
        p <- p + ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$group,
                                                                shape = .data$group),
                                     size = node_size, colour = unname(group_cols[as.character(nodes$group)]), stroke = 1.8) +
            ggplot2::scale_fill_manual(values = group_tints, name = NULL) +
            ggplot2::scale_shape_manual(values = stats::setNames(.PLOT_GROUP_SHAPES[seq_along(group_cols)], names(group_cols)),
                                        name = NULL) +
            ggplot2::guides(fill = ggplot2::guide_legend(order = 1, override.aes = list(colour = unname(group_cols), size = 5)),
                            shape = ggplot2::guide_legend(order = 1))
    }

    p <- p + ggplot2::geom_text(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
                                size = 0.28 * node_size, colour = .PLOT_INK) +
        ggplot2::coord_equal(clip = "off") +
        ggplot2::labs(caption = rule) +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                       plot.background = ggplot2::element_rect(fill = .PLOT_SURFACE, colour = NA),
                       plot.caption = ggplot2::element_text(colour = "gray35", hjust = 0.5),
                       plot.margin = ggplot2::margin(20, 20, 10, 20))
    return(p)
}


#' evidence_class (internal)
#' @description Evidence class of each Savage-Dickey Bayes factor BF_01, from the intervals in .PLOT_BF_INTERVALS.
#' @param bf01 Savage-Dickey Bayes factors BF_01
#' @return a factor with the names of .PLOT_BF_INTERVALS as levels
#' @noRd
.evidence_class <- function(bf01) {
    k <- findInterval(bf01, c(0, .PLOT_BF_INTERVALS[-length(.PLOT_BF_INTERVALS)]))
    return(factor(names(.PLOT_BF_INTERVALS)[k], levels = names(.PLOT_BF_INTERVALS)))
}


#' plot_bf (internal)
#' @description Bayes factor plot of plot.dmrfit(type = "bf"); see its documentation for the arguments.
#' @return a ggplot object
#' @noRd
.plot_bf <- function(x, estimate) {
    sd <- x$savage_dickey
    if (is.null(sd)) {
        stop("type = \"bf\" needs Savage-Dickey Bayes factors: refit with dmrfit(..., with_prior = TRUE, savage_dickey = TRUE).")
    }
    var_names <- if (is.null(x$var_names)) paste0("V", seq_len(x$P)) else x$var_names
    inter <- .interaction_estimates(x, estimate)

    bf01 <- rep(NA_real_, nrow(inter))
    bf01[match(names(sd$bf_01), inter$name)] <- unname(sd$bf_01)
    circles <- data.frame(row = factor(var_names[inter$i], levels = rev(var_names[-1])),
                          col = factor(var_names[inter$j], levels = var_names[-length(var_names)]),
                          evidence = .evidence_class(bf01),
                          weight = abs(inter$estimate), stringsAsFactors = FALSE)
    circles <- circles[!is.na(bf01), , drop = FALSE] # no circle for the absent edges of a constrained fit

    size_lab <- if (inherits(x, "dmrfit_bayes")) paste0("|posterior ", estimate, "|") else "|estimate|"
    p <- ggplot2::ggplot(circles, ggplot2::aes(x = .data$col, y = .data$row)) +
        # a thin outline keeps the light classes visible on the light background
        ggplot2::geom_point(ggplot2::aes(size = .data$weight, fill = .data$evidence), shape = 21, colour = "gray35",
                            stroke = 0.3, show.legend = TRUE) +
        ggplot2::scale_fill_manual(values = .PLOT_EVIDENCE_COLORS, labels = .PLOT_EVIDENCE_LABELS,
                                   name = expression("Evidence (" * BF["01"] * ")"), drop = FALSE) +
        ggplot2::scale_size_area(max_size = 9, name = size_lab) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::coord_equal() +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "gray92", linewidth = 0.4),
                       panel.grid.minor = ggplot2::element_blank(), legend.position = "right",
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.background = ggplot2::element_rect(fill = .PLOT_SURFACE, colour = NA)) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1, override.aes = list(size = 5)),
                        size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "gray65")))
    return(p)
}


#' par_labels (internal)
#' @description Readable labels of parameter names: the interaction sigma of variables i and j becomes
#'   "<var j>-<var i>", and threshold h of variable p (mu) becomes "<var p>: threshold h".
#' @param par_names parameter names, as in names(x$argument)
#' @param var_names the variable names
#' @return a character vector
#' @noRd
.par_labels <- function(par_names, var_names) {
    idx <- regmatches(par_names, gregexpr("[0-9]+", par_names))
    labels <- vapply(seq_along(par_names), function(k) {
        ij <- as.integer(idx[[k]])
        if (startsWith(par_names[k], "sigma")) paste0(var_names[ij[2]], "-", var_names[ij[1]])
        else paste0(var_names[ij[1]], ": threshold ", ij[2])
    }, character(1))
    return(labels)
}


#' bind_panels (internal)
#' @description Binds a list of data frames, one per panel, with the panel labels as an ordered factor.
#' @param frames list of data frames with a column parameter
#' @param labels the panel labels, in order
#' @return a data frame
#' @noRd
.bind_panels <- function(frames, labels) {
    out <- do.call(rbind, frames)
    out$parameter <- factor(out$parameter, levels = labels)
    return(out)
}


#' plot_draws (internal)
#' @description Trace and density plots of plot.dmrfit(type = "trace" or "density"); see its documentation for the
#'   arguments.
#' @param type "trace" or "density"
#' @return a ggplot object
#' @noRd
.plot_draws <- function(x, type, estimate, pars, prob) {
    var_names <- if (is.null(x$var_names)) paste0("V", seq_len(x$P)) else x$var_names
    par_names <- names(x$argument)
    summary_fun <- switch(estimate, mode = .posterior_mode, mean = mean, median = stats::median)

    # --- parameters to show: by default the four interactions with the largest absolute posterior summary ---
    if (is.null(pars)) {
        inter <- .interaction_estimates(x, estimate)
        pars <- inter$name[order(-abs(inter$estimate))][seq_len(min(4, nrow(inter)))]
    } else if (is.numeric(pars)) {
        if (any(pars < 1 | pars > length(par_names) | pars != round(pars))) {
            stop("pars must be positions between 1 and ", length(par_names), ".")
        }
        pars <- par_names[pars]
    } else if (!all(pars %in% par_names)) {
        stop("unknown parameter(s) in pars: ", paste(setdiff(pars, par_names), collapse = ", "),
             ". See names(x$argument).")
    }
    pars <- unique(pars)
    if (length(pars) > 9) {
        stop("At most 9 parameters can be shown; pars has ", length(pars), ".")
    }
    labels <- .par_labels(pars, var_names)
    draws <- x$draws[match(pars, par_names), , drop = FALSE]

    long <- data.frame(parameter = factor(rep(labels, each = ncol(draws)), levels = labels),
                       iteration = rep(seq_len(ncol(draws)), times = length(pars)),
                       value = as.vector(t(draws)))
    n_col <- if (length(pars) <= 3) length(pars) else if (length(pars) == 4) 2 else 3 # rows of up to three panels

    if (type == "trace") {
        p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$iteration, y = .data$value)) +
            ggplot2::geom_line(colour = .PLOT_TRACE, linewidth = 0.25) +
            ggplot2::labs(x = "Iteration (after burn-in)", y = NULL)
    } else {
        # --- kernel density of every parameter, its posterior summary and its HPD interval ---
        curves <- list()
        marks <- list()
        for (k in seq_along(pars)) {
            z <- draws[k, ]
            d <- stats::density(z)
            hpd <- .hdi(z, prob)
            curves[[k]] <- data.frame(parameter = labels[k], x = d$x, y = d$y)
            marks[[k]] <- data.frame(parameter = labels[k], center = summary_fun(z), lower = hpd[1], upper = hpd[2])
        }
        curves <- .bind_panels(curves, labels)
        marks <- .bind_panels(marks, labels)
        p <- ggplot2::ggplot(curves, ggplot2::aes(x = .data$x, y = .data$y)) +
            ggplot2::geom_area(fill = .PLOT_DENSITY_SHADE) +
            ggplot2::geom_line(colour = .PLOT_DENSITY, linewidth = 0.5) +
            ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = .data$center), colour = .PLOT_DENSITY,
                                linewidth = 0.5) +
            ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = .data$lower), colour = .PLOT_DENSITY,
                                linewidth = 0.4, linetype = "dashed") +
            ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = .data$upper), colour = .PLOT_DENSITY,
                                linewidth = 0.4, linetype = "dashed") +
            ggplot2::labs(x = NULL, y = "Density",
                          caption = paste0("Solid line: posterior ", estimate, ". Dashed lines: ", round(100 * prob),
                                           "% highest posterior density interval"))
    }

    p <- p + ggplot2::facet_wrap(~ parameter, ncol = n_col, scales = "free") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(face = "bold"),
                       plot.background = ggplot2::element_rect(fill = .PLOT_SURFACE, colour = NA),
                       plot.caption = ggplot2::element_text(colour = "gray35", hjust = 0.5))
    return(p)
}


#' plot_intervals (internal)
#' @description Interval plot of plot.dmrfit(type = "intervals"); see its documentation for the arguments.
#' @return a ggplot object
#' @noRd
.plot_intervals <- function(x, estimate, pars, prob) {
    var_names <- if (is.null(x$var_names)) paste0("V", seq_len(x$P)) else x$var_names
    par_names <- names(x$argument)
    if (is.numeric(pars)) pars <- par_names[pars]

    # --- intervals (confint() checks pars) and point estimates ---
    if (inherits(x, "dmrfit_bayes")) {
        bounds <- confint(x, parm = pars, level = prob)
        summary_fun <- switch(estimate, mode = .posterior_mode, mean = mean, median = stats::median)
        center <- apply(x$draws[match(rownames(bounds), par_names), , drop = FALSE], 1, summary_fun)
        rule <- paste0("Point: posterior ", estimate, ". Line: ", round(100 * prob), "% highest posterior density interval")
    } else {
        lrt <- !is.null(x$lrt_intervals)
        level <- .fit_level(x)
        bounds <- confint(x, parm = pars, level = level, method = if (lrt) "lrt" else "wald")
        center <- x$argument[rownames(bounds)]
        rule <- paste0("Point: estimate. Line: ", round(100 * level), "% ",
                       if (lrt) "likelihood-ratio interval" else "Wald interval (sandwich standard errors)")
    }

    is_inter <- startsWith(rownames(bounds), "sigma")
    iv <- data.frame(label = .par_labels(rownames(bounds), var_names), center = unname(center),
                     lower = bounds[, 1], upper = bounds[, 2],
                     type = factor(ifelse(is_inter, "Interactions", "Thresholds"), levels = c("Thresholds", "Interactions")))
    iv$label <- factor(iv$label, levels = iv$label[order(iv$type, iv$center)]) # sorted by estimate within each panel

    p <- ggplot2::ggplot(iv, ggplot2::aes(x = .data$center, y = .data$label)) +
        ggplot2::geom_vline(xintercept = 0, colour = "gray80", linewidth = 0.4) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), width = 0, orientation = "y",
                               colour = .PLOT_DENSITY, linewidth = 0.5) +
        ggplot2::geom_point(colour = .PLOT_DENSITY, size = 1.8) +
        ggplot2::facet_wrap(~ type, scales = "free", ncol = 2) + # own axes: thresholds and interactions differ in scale
        ggplot2::labs(x = NULL, y = NULL, caption = rule) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(face = "bold"),
                       plot.background = ggplot2::element_rect(fill = .PLOT_SURFACE, colour = NA),
                       plot.caption = ggplot2::element_text(colour = "gray35", hjust = 0.5))
    return(p)
}


#' expected_influence (internal)
#' @description Expected influence of every variable, EI_i = sum_j theta_ij: for a dmrfit_bayes fit, computed on every
#'   draw, summarized by the posterior summary, the HPD interval and the posterior probability of the largest expected
#'   influence; for a dmrfit fit, the estimate with its Wald interval, the standard error from the sandwich covariance
#'   of the interactions (sqrt(a' Sigma a), with a the indicator of the interactions of the variable).
#' @param x a dmrfit or dmrfit_bayes object
#' @param estimate the posterior summary for a dmrfit_bayes fit: "mode", "mean" or "median"
#' @param prob mass of the HPD interval of a dmrfit_bayes fit
#' @return a data frame with one row per variable: name, center, lower, upper and p_most_central (NA for a dmrfit fit)
#' @noRd
.expected_influence <- function(x, estimate, prob) {
    var_names <- if (is.null(x$var_names)) paste0("V", seq_len(x$P)) else x$var_names
    inter <- .interaction_estimates(x, "mean") # names and node indices of the interactions
    inter_idx <- match(inter$name, names(x$argument))

    # --- indicator of the interactions of each variable: EI = A theta ---
    A <- matrix(0, x$P, nrow(inter))
    A[cbind(inter$i, seq_len(nrow(inter)))] <- 1
    A[cbind(inter$j, seq_len(nrow(inter)))] <- 1

    if (inherits(x, "dmrfit_bayes")) {
        ei_draws <- A %*% x$draws[inter_idx, , drop = FALSE]
        summary_fun <- switch(estimate, mode = .posterior_mode, mean = mean, median = stats::median)
        hpd <- t(apply(ei_draws, 1, .hdi, prob = prob))
        most_central <- table(factor(apply(ei_draws, 2, which.max), levels = seq_len(x$P)))
        out <- data.frame(name = var_names, center = apply(ei_draws, 1, summary_fun), lower = hpd[, 1], upper = hpd[, 2],
                          p_most_central = as.vector(most_central) / ncol(ei_draws))
    } else {
        theta <- x$argument[inter_idx]
        Sigma <- x$utils$HW[inter_idx, inter_idx, drop = FALSE] # zero for the absent edges of a constrained fit
        se <- sqrt(diag(A %*% Sigma %*% t(A)))
        z <- stats::qnorm((1 + .fit_level(x)) / 2)
        center <- as.vector(A %*% theta)
        out <- data.frame(name = var_names, center = center, lower = center - z * se, upper = center + z * se,
                          p_most_central = NA_real_)
    }
    return(out)
}


#' plot_centrality (internal)
#' @description Centrality plot of plot.dmrfit(type = "centrality"); see its documentation for the arguments.
#' @return a ggplot object
#' @noRd
.plot_centrality <- function(x, estimate, prob) {
    ei <- .expected_influence(x, estimate, prob)
    ei$name <- factor(ei$name, levels = ei$name[order(ei$center)]) # sorted by expected influence
    bayes <- inherits(x, "dmrfit_bayes")
    rule <- if (bayes) {
        paste0("Point: posterior ", estimate, ". Line: ", round(100 * prob), "% highest posterior density interval\n",
               "Pr(most central): posterior probability of the largest expected influence")
    } else {
        paste0("Point: estimate. Line: ", round(100 * .fit_level(x)), "% Wald interval (sandwich standard errors)")
    }

    p <- ggplot2::ggplot(ei, ggplot2::aes(x = .data$center, y = .data$name)) +
        ggplot2::geom_vline(xintercept = 0, colour = "gray80", linewidth = 0.4) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), width = 0, orientation = "y",
                               colour = .PLOT_DENSITY, linewidth = 0.5) +
        ggplot2::geom_point(colour = .PLOT_DENSITY, size = 1.8) +
        ggplot2::labs(x = expression("Expected influence " * sum(theta[ij], j)), y = NULL, caption = rule) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
                       plot.background = ggplot2::element_rect(fill = .PLOT_SURFACE, colour = NA),
                       plot.caption = ggplot2::element_text(colour = "gray35", hjust = 0.5))
    if (bayes) {
        # posterior probability of the largest expected influence, in a column to the right of the panel
        p <- p + ggplot2::geom_text(ggplot2::aes(x = Inf, label = sprintf("%.2f", .data$p_most_central)), hjust = -0.3,
                                    size = 3, colour = .PLOT_INK) +
            ggplot2::annotate("text", x = Inf, y = Inf, label = "Pr(most\ncentral)", hjust = -0.1, vjust = -0.3, size = 3,
                              colour = "gray35", lineheight = 0.9) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::theme(plot.margin = ggplot2::margin(30, 60, 5.5, 5.5))
    }
    return(p)
}
