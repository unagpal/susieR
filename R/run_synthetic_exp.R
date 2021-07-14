setwd("~/SuSiE-Ann/susieR/R/")


#install.packages("mvtnorm", repos="http://cran.us.r-project.org")
#install.packages("matrixcalc", repos="http://cran.us.r-project.org")

#library(mvtnorm)
#library(matrixcalc)

source("susie_ann_utils.R")

PVE <- 0.6
n <- 100
p <- 1000
L <- 1
G <- 3000
K <- 20
annot_sparsity <- 0.1
susie_ann_opt_itr <- 10
print("Now calling geneneration of synthetic data")
synth_data <- generate_synthetic_data_final(PVE, n, p, L, G, K, annot_sparsity)
fname <- ".Rdata"
fname = paste("synth_main_data", fname, sep="")
fname = paste("~/SuSiE-Ann/", fname)
save(synth_data, file=fname)
print("Saved synthetic data")
#synth_data_ld <- get(load("~/SuSiE-Ann/ synth_main_data.Rdata"))

# install.packages("rootSolve", repos="http://cran.us.r-project.org")
# install.packages("faux", repos="http://cran.us.r-project.org")
# install.packages("LaplacesDemon", repos="http://cran.us.r-project.org")
# install.packages("clusterGeneration", repos="http://cran.us.r-project.org")
# install.packages("SimDesign", repos="http://cran.us.r-project.org")
# install.packages("matrixStats", repos="http://cran.us.r-project.org")
# install.packages("qrnn", repos="http://cran.us.r-project.org")
# install.packages("numDeriv", repos="http://cran.us.r-project.org")
# install.packages("wordspace", repos="http://cran.us.r-project.org")
# install.packages("philentropy", repos="http://cran.us.r-project.org")
# library(rootSolve)
# library(faux)
# library(LaplacesDemon)
# library(clusterGeneration)
# library(SimDesign)
# library(matrixStats)
# library(qrnn)
# library(numDeriv)
# library(wordspace)
# library(philentropy)
source("susie_ann.R")

susie_ann_X <- synth_data$susie_ann_X
susie_ann_Y <- synth_data$susie_ann_Y
synth_A = synth_data$A
res <- susie_ann(susie_ann_X, susie_ann_Y, synth_A, G, FALSE, G, L=5, susie_ann_opt_itr = susie_ann_opt_itr,
                 elbo_optim_method="L-BFGS-B", elbo_opt_steps_per_itr = 100)
fname <- ".Rdata"
fname = paste("synth_main_res", fname, sep="")
fname = paste("~/SuSiE-Ann/", fname)
save(res, file=fname)
#analyze_susie_ann_fit(res, synth_data)
