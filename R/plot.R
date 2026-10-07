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
#' }
#'
#' @param x a \code{dmrfit} or \code{dmrfit_bayes} object.
#' @param type the plot: \code{"network"} (default) or \code{"bf"}.
#' @param estimate the posterior summary of the interactions of a \code{dmrfit_bayes} fit: \code{"mode"} (default) the
#'   marginal posterior mode, the maximum of the kernel density estimate of the draws (\code{stats::density}),
#'   \code{"mean"} or \code{"median"}. The mode is the default because the posterior can be asymmetric in small
#'   samples, which moves the mean and the median away from the mode; Arena and Marsman (2026) found this for the
#'   thresholds of rarely chosen response categories. Ignored for a \code{dmrfit} fit, which has point estimates.
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
#' @method plot dmrfit
#' @export
#'
plot.dmrfit <- function(x, type = c("network", "bf"), estimate = c("mode", "mean", "median"), all_edges = FALSE,
                        groups = NULL, layout = c("fr", "circle"), seed = 30, node_size = 10, ...) {

    type <- match.arg(type)
    if (!inherits(x, "dmrfit")) {
        stop("x is not of class 'dmrfit'")
    }
    if (!missing(estimate) && !inherits(x, "dmrfit_bayes")) {
        warning("estimate is ignored for a dmrfit fit, which has point estimates.")
    }
    estimate <- match.arg(estimate)
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
                bf = .plot_bf(x, estimate = estimate))
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
