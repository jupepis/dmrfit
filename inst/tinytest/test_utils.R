# Tests for the internal R helpers used by the plots.

# --- HPD interval (.hdi): 95% of the draws, the shortest such interval, shifted to the mode for skewed draws
set.seed(2)
y <- rnorm(1e5)
h <- dmrfit:::.hdi(y)
expect_equal(mean(y >= h[1] & y <= h[2]), 0.95, tolerance = 1e-4)
expect_equal(h, qnorm(c(0.025, 0.975)), tolerance = 0.02)
z <- rexp(1e5)
expect_true(dmrfit:::.hdi(z)[1] < quantile(z, 0.025))
expect_equal(dmrfit:::.hdi(1:2), c(NA_real_, NA_real_))

# --- marginal posterior mode (.posterior_mode): the peak of the kernel density estimate, NA for constant draws
set.seed(3)
expect_equal(dmrfit:::.posterior_mode(rgamma(1e5, shape = 3)), 2, tolerance = 0.05) # mode of Gamma(3, 1) is 2
expect_true(is.na(dmrfit:::.posterior_mode(rep(1, 10))))
