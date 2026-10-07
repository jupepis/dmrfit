# Tests for plot.dmrfit(): the edges drawn by the network plot (Savage-Dickey rule, constrained fits), the node
# groups and the layout, and the argument checks. The plots are built (ggplot_build) but not drawn.

data(rads2, package = "dmrfit")
X <- rads2[, c("D3", "D6", "D7", "D8", "D16")]
fit <- dmrfit(X, with_prior = TRUE, savage_dickey = TRUE, M = 500)
edge_layer <- function(p) ggplot2::layer_data(p, 1)  # the curves are the first layer

# --- the item names are stored and used as node labels
expect_equal(fit$var_names, colnames(X))
p <- plot(fit)
expect_inherits(p, "ggplot")
expect_silent(ggplot2::ggplot_build(p))
expect_equal(sort(ggplot2::layer_data(p, 3)$label), sort(colnames(X)))

# --- Savage-Dickey rule: the included edges (BF_01 < 1/10) are drawn, or every edge with all_edges = TRUE
sd <- fit$savage_dickey
expect_equal(nrow(edge_layer(plot(fit))), sum(sd$bf_01 < 1 / 10))
expect_equal(nrow(edge_layer(plot(fit, all_edges = TRUE))), choose(5, 2))

# --- without Bayes factors every estimated interaction is drawn, with a message; a constrained fit draws its free edges
fit_nobf <- dmrfit(X, with_prior = TRUE)
expect_message(p_nobf <- plot(fit_nobf), "every estimated interaction")
expect_equal(nrow(edge_layer(p_nobf)), choose(5, 2))
S <- matrix(1, 5, 5); diag(S) <- 0; S[1, 2] <- S[2, 1] <- 0; S[3, 5] <- S[5, 3] <- 0
fit_s <- dmrfit(X, structure = S)
expect_equal(nrow(edge_layer(plot(fit_s, all_edges = TRUE))), choose(5, 2) - 2)

# --- groups: in column order or named, at most four
p_g <- plot(fit, groups = c("a", "a", "b", "b", "c"))
expect_equal(nrow(ggplot2::layer_data(p_g, 2)), 5L)
expect_silent(ggplot2::ggplot_build(plot(fit, groups = setNames(c("a", "b", "a", "b", "a"), rev(colnames(X))))))
expect_error(plot(fit, groups = c("a", "b")), "one entry per variable")
expect_error(plot(fit, groups = letters[1:5]), "At most 4 groups")

# --- layout: the "fr" layout comes from all estimated interactions, so the nodes keep their positions when the
# threshold changes; a different seed gives a different arrangement; "circle" puts the nodes on the unit circle
node_xy <- function(p) ggplot2::layer_data(p, 2)[, c("x", "y")]
expect_equal(node_xy(plot(fit)), node_xy(plot(fit, all_edges = TRUE)))
expect_equal(node_xy(plot(fit)), node_xy(plot(fit, seed = 30)))
expect_false(isTRUE(all.equal(node_xy(plot(fit)), node_xy(plot(fit, seed = 20)))))
circle <- node_xy(plot(fit, layout = "circle"))
expect_equal(circle$x^2 + circle$y^2, rep(1, 5))

# --- layout given by the user: centered and rescaled by the same factor on both axes (proportions kept), in column
# order or by row names
L <- cbind(1:5, c(2, 4, 6, 8, 10))
expect_equal(node_xy(plot(fit, layout = L)), data.frame(x = seq(-0.5, 0.5, by = 0.25), y = seq(-1, 1, by = 0.5)))
L_named <- L[5:1, ]; rownames(L_named) <- rev(colnames(X))
expect_equal(node_xy(plot(fit, layout = L_named)), node_xy(plot(fit, layout = L)))
expect_error(plot(fit, layout = L[, 1, drop = FALSE]), "two columns")
expect_error(plot(fit, layout = "grid"))

# --- the groups follow the levels of a factor; node_size sets the node size and scales the labels
g <- factor(c("a", "a", "b", "b", "c"), levels = c("c", "b", "a"))
expect_equal(ggplot2::ggplot_build(plot(fit, groups = g))$plot$scales$get_scales("fill")$get_limits(), c("c", "b", "a"))
expect_equal(unique(ggplot2::layer_data(plot(fit, node_size = 14), 2)$size), 14)
expect_equal(unique(ggplot2::layer_data(plot(fit, node_size = 14), 3)$size), 0.28 * 14)
expect_error(plot(fit, node_size = 0), "positive number")
expect_error(plot(fit, seed = "a"), "single number")

# --- the caller's random number stream is left untouched by the igraph layout
set.seed(5); before <- .Random.seed
invisible(plot(fit))
expect_identical(.Random.seed, before)

# --- argument checks
expect_error(plot(fit, type = "other"))
expect_error(plot(fit, all_edges = NA), "TRUE or FALSE")

# --- dmrfit_bayes: the edges are the marginal posterior modes (default), means or medians of the draws
fit_b <- suppressWarnings(dmrfit_bayes(X, nsim = 400, burnin = 200, progress = FALSE))
inter_rows <- grep("^sigma", rownames(fit_b$draws))
est <- function(e) { d <- dmrfit:::.interaction_estimates(fit_b, e); setNames(d$estimate, d$name) }
expect_equal(est("mean"), apply(fit_b$draws[inter_rows, ], 1, mean))
expect_equal(est("median"), apply(fit_b$draws[inter_rows, ], 1, median))
expect_equal(est("mode"), apply(fit_b$draws[inter_rows, ], 1, function(z) { d <- density(z); d$x[which.max(d$y)] }))
expect_silent(ggplot2::ggplot_build(plot(fit_b, estimate = "median")))
expect_warning(plot(fit, estimate = "mean"), "ignored for a dmrfit fit")


# --- evidence classes of BF_01 (evidence for exclusion): each interval runs from the previous bound up to, but
# excluding, its own
cls <- dmrfit:::.evidence_class(c(0.05, 1/10, 0.2, 1/3, 1, 3, 5, 10, 20))
expect_equal(as.character(cls), c("Included", "Weak included", "Weak included", "Inconclusive", "Inconclusive",
                                  "Weak excluded", "Weak excluded", "Excluded", "Excluded"))

# --- Bayes factor plot: one circle per pair, colored by the evidence class, with area |estimate|
bf01 <- fit$savage_dickey$bf_01
p_bf <- plot(fit, type = "bf")
circ <- p_bf$data
expect_equal(nrow(circ), choose(5, 2))
ord <- match(paste0("sigma[", match(circ$row, colnames(X)), ",", match(circ$col, colnames(X)), "]"), names(fit$argument))
expect_equal(circ$weight, unname(abs(fit$argument[ord])))
expect_equal(as.character(circ$evidence), as.character(dmrfit:::.evidence_class(unname(bf01[names(fit$argument)[ord]]))))
expect_equal(levels(circ$evidence), c("Included", "Weak included", "Inconclusive", "Weak excluded", "Excluded"))
expect_silent(ggplot2::ggplot_build(p_bf))
# a constrained fit has no circle for its absent edges; checks
fit_s_bf <- dmrfit(X, structure = S, with_prior = TRUE, savage_dickey = TRUE, M = 500)
expect_equal(nrow(plot(fit_s_bf, type = "bf")$data), choose(5, 2) - 2)
expect_error(plot(fit_nobf, type = "bf"), "savage_dickey = TRUE")
expect_error(plot(fit, type = "bf_matrix"))
# for dmrfit_bayes, the circle area is the chosen posterior summary
expect_equal(sort(plot(fit_b, type = "bf", estimate = "mean")$data$weight), sort(unname(abs(est("mean")))))
