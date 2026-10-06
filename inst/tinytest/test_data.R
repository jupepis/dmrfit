# Tests for the rads2 dataset: dimensions, coding and the item clusters stored with it.
# The examples, the README and the help page rely on these properties.

data(rads2, package = "dmrfit")

# --- dimensions and coding
expect_inherits(rads2, "data.frame")
expect_equal(dim(rads2), c(917L, 25L))
expect_true(all(vapply(rads2, is.integer, logical(1))))
expect_false(anyNA(rads2))
expect_true(all(unlist(rads2) %in% 1:4))
# every item uses all four categories (the ordinal models assume 4 categories per item)
expect_true(all(vapply(rads2, function(v) length(unique(v)), integer(1)) == 4L))

# --- item numbering of the RADS-2 (items 2, 10, 23, 25 and 29 are not in the four-factor version)
expect_equal(names(rads2), paste0("D", c(1, 3:9, 11:22, 24, 26:28, 30)))

# --- item clusters
cl <- attr(rads2, "clusters")
expect_true(is.character(cl))
expect_equal(names(cl), names(rads2))
expect_equal(as.vector(table(cl)[c("Dysphoria", "Anhedonia/Negative Affect", "Negative Self-Assessment", "Somatic Complaint")]),
             c(7L, 3L, 8L, 7L))
expect_equal(names(which(cl == "Dysphoria")), c("D3", "D6", "D7", "D8", "D16", "D21", "D26"))
