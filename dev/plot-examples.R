# Self-contained gallery of the plot() methods of dmrfit, on the RADS-2 data shipped with the package.
# Run from the package root: Rscript dev/plot-examples.R [output.pdf]
# Every plot is one page of the output pdf (default: dev/plot-examples.pdf).

library(dmrfit)
library(ggplot2)

out_file <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/plot-examples.pdf"

data(rads2)
clusters <- attr(rads2, "clusters")
items <- names(clusters)[clusters %in% c("Dysphoria", "Anhedonia/Negative Affect")] # 10 items, two clusters

# --- fits -------------------------------------------------------------------------------------------------------------
fit_bf   <- dmrfit(rads2[, items], with_prior = TRUE, savage_dickey = TRUE) # point estimates + SIR Bayes factors
fit_nobf <- dmrfit(rads2[, items], with_prior = TRUE)                       # point estimates only
S <- matrix(1, length(items), length(items), dimnames = list(items, items)); diag(S) <- 0
S[clusters[items] == "Dysphoria", clusters[items] != "Dysphoria"] <- 0      # no edges between the two clusters
S[clusters[items] != "Dysphoria", clusters[items] == "Dysphoria"] <- 0
fit_constrained <- dmrfit(rads2[, items], structure = S, with_prior = TRUE)
fit_lrt  <- dmrfit(rads2[, items], with_prior = TRUE, lrt_intervals = TRUE)  # likelihood-ratio intervals
fit_bayes <- dmrfit_bayes(rads2[, items], nsim = 2000, burnin = 1000, progress = FALSE) # CoRe posterior draws

# --- plots ------------------------------------------------------------------------------------------------------------
pages <- list(
    "1. dmrfit, included edges (BF01 < 1/10)"                  = plot(fit_bf),
    "2. dmrfit, all interactions (all_edges = TRUE)"           = plot(fit_bf, all_edges = TRUE),
    "3. dmrfit without Bayes factors (message, all edges)"     = plot(fit_nobf),
    "4. dmrfit, nodes grouped by cluster"                      = plot(fit_bf, groups = clusters[items]),
    "5. dmrfit, circular layout given by the user"             = plot(fit_bf, groups = clusters[items],
                                                                      layout = cbind(cos(2 * pi * seq_along(items) / length(items)),
                                                                                     sin(2 * pi * seq_along(items) / length(items)))),
    "6. dmrfit, constrained structure (no between-cluster edges)" = plot(fit_constrained, groups = clusters[items]),
    "7. dmrfit_bayes, posterior modes (default)"               = plot(fit_bayes, groups = clusters[items]),
    "8. dmrfit_bayes, posterior means (estimate = \"mean\")"  = plot(fit_bayes, groups = clusters[items], estimate = "mean"),
    "9. other arrangement (seed = 20 instead of 30)"          = plot(fit_bf, groups = clusters[items], seed = 20),
    "10. circular layout, larger nodes (node_size = 13)"      = plot(fit_bf, groups = clusters[items], layout = "circle",
                                                                      node_size = 13),
    "11. Bayes factors (dmrfit, SIR)"                         = plot(fit_bf, type = "bf"),
    "12. Bayes factors (dmrfit_bayes, posterior mode)"        = plot(fit_bayes, type = "bf"),
    "13. trace (default: the 4 interactions with the largest |posterior mode|)" = plot(fit_bayes, type = "trace"),
    "14. density, thresholds of D1 and interactions of D1 (9 panels)" = plot(fit_bayes, type = "density",
                                                                     pars = c(paste0("mu[1,", 1:3, "]"), paste0("sigma[", 2:7, ",1]"))),
    "15. density, posterior mean and 90% HPD interval"       = plot(fit_bayes, type = "density", estimate = "mean", prob = 0.9),
    "16. intervals (dmrfit, Wald)"                            = plot(fit_bf, type = "intervals"),
    "17. intervals (dmrfit, likelihood-ratio)"                = plot(fit_lrt, type = "intervals"),
    "18. intervals (dmrfit_bayes, HPD)"                       = plot(fit_bayes, type = "intervals"),
    "19. centrality (dmrfit, Wald)"                           = plot(fit_bf, type = "centrality"),
    "20. centrality (dmrfit_bayes, HPD and Pr(most central))"  = plot(fit_bayes, type = "centrality"),
    "21. the returned ggplot can be modified"                 = plot(fit_bf, groups = clusters[items]) +
                                                                  labs(title = "RADS-2: dysphoria and anhedonia items") +
                                                                  theme(legend.position = "right")
)

pdf(out_file, width = 7.5, height = 8)
for (k in seq_along(pages)) {
    print(pages[[k]] + labs(subtitle = names(pages)[k]))
}
invisible(dev.off())
cat("Wrote", length(pages), "plots to", out_file, "\n")
