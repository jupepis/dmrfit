# Precompile the vignette: run the code of vignettes/dmrfit.Rmd.orig and write vignettes/dmrfit.Rmd, with the printed
# output and the figures (vignettes/dmrfit-*.png). The package is then checked without running the vignette code.
# Run from the package root, with the current version of the package installed, after changes that affect the output:
#   Rscript vignettes/precompile.R

old_wd <- setwd("vignettes")
unlink(list.files(pattern = "^dmrfit-.*\\.png$")) # figures of the previous version
knitr::knit("dmrfit.Rmd.orig", output = "dmrfit.Rmd")
setwd(old_wd)
