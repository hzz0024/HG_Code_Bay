devtools::install_github("samyeaman/dgconstraint", build_vignettes = TRUE)
library(dgconstraint)
browseVignettes(package = "dgconstraint")
load("/Users/ryan/Downloads/dgconstraint-master/data/copper.RData")

copper[,1]

single_c_chisq(copper[,1], copper[,2], num_permute = 10000, na.rm = F)
single_p_chisq(copper[,1], copper[,2], num_permute = 10000, na.rm = F)

pairwise_c_chisq(copper, num_permute = 10000, na.rm = F)
allwise_p_chisq(copper, num_permute = 10000, na.rm = F)

single_c_hyper(copper[,1], copper[,2], na.rm = F)
single_p_hyper(copper[,1], copper[,2], na.rm = F)

pairwise_c_hyper(copper, na.rm = F)