source("SuSiE-Ann/susieR/R/susie_ann.R")
#generate_synthetic_data_basic(pve, n, p, L, num_loci, num_annotations)

#synth_data <- generate_synthetic_data_basic(0.4, 10, 2, 3, num_loci, 4)
#susie_ann_X <- list()
#susie_ann_Y <- list()
#for (loc in 1:num_loci){
#  susie_ann_X[[loc]] <- synth_data$output[[loc]]$X
#  susie_ann_Y[[loc]] <- synth_data$output[[loc]]$Y
#}
#synth_A = synth_data$A
#num_loci <- 5
susie_ann_opt_itr <- 50
num_trials <- 1
all_pi_1_true <- rep(0, num_trials)
all_pi_1_est <- rep(0, num_trials)
n <- 10
p <- 2
num_annotations <- 4
for (pve in c(0.4)){
  for (num_loci in c(5)){
    for (l in c(2)){
      print("PVE: ")
      print(pve)
      print("Num loci:")
      print(num_loci)
      print("L:")
      print(l)
      synth_data <- generate_synthetic_data_basic(pve, n, p, l, num_loci, num_annotations)
      susie_ann_X <- list()
      susie_ann_Y <- list()
      for (loc in 1:num_loci){
        susie_ann_X[[loc]] <- synth_data$output[[loc]]$X
        susie_ann_Y[[loc]] <- synth_data$output[[loc]]$Y
      }
      synth_A = synth_data$A
      res <- susie_ann(susie_ann_X, susie_ann_Y, synth_A, num_loci, FALSE, num_loci, L=l, susie_ann_opt_itr = susie_ann_opt_itr, elbo_opt_steps_per_itr = 10)
      print("---------------FINAL RESULTS-----------------")
      #print("Synthetic data:")
      #print(synth_data)
      #print("Final optimized pi:")
      #print(res$final_pi)
      all_pi_1_true[t] <- synth_data$pi[1]
      all_pi_1_est[t] <- res$final_pi[1]
      print("KL between pi:")
      print(KLD(res$final_pi, synth_data$pi)$sum.KLD.px.py)
      print("Predicted pi:")
      print(res$final_pi)
      print("True pi:")
      print(synth_data$pi)
      print("All iter pi:")
      print(res$all_iter_pi)
      
      #plot(res$final_pi[1], synth_data$pi[1])
      plot(susie_ann_Y[[1]],predict(res$susie_model[[1]]))
      print("Annotations matrix A:")
      print(synth_A)
      #res$susie_model
      #print("True pi:")
      #print(synth_data$pi)
      #print("True w")
      #print(synth_data$w)
      #print("Final optimized w:")
      #print(res$w)
      #print("Plotting ELBO")
      #plot(1:susie_ann_opt_itr, res$elbo_values)
    }
  }
}
#plot(all_pi_1_est, all_pi_1_true)
plot(1:length(res$all_iter_pi), res$all_iter_pi)
#plot (1:length(res$all_iter_alpha), res$all_iter_alpha)

# #Test Case 1: SuSiE-Ann example based on SuSiE vignette
# set.seed(1)
# n = 20
# p = 3
# b = rep(0,p)
# A = matrix(c(0,1,1,0,1,1), nrow=3, ncol=2)
# print("Defined matrix A")
# b[2] = 3
# X = matrix(rnorm(n*p),nrow=n,ncol=p)
# y = X %*% b + rnorm(n)
# print(dim(X))
# print(dim(y))
# #print(crossprod(y, X))
# print("Defined X and Y; now calling susie_ann")
# susie_ann_res <- susie_ann(X,y,A, FALSE, rho=0.1, L=1)
# res = susie_ann_res$susie_model
# #coef(res)
# plot(y,predict(res))
# plot(1:100, susie_ann_res$elbo_values)
#print(res)

# Test Case 2: Toy data for p=2
#(beta, corr, noise_ratio, annot_signal_ratio, n)
#toy_data_1 <- generate_toy_data(1.0, 0.9, 0.5, 5, 5)
#toy_data_2 <- generate_toy_data(1.0, 0.9, 0.5, 5, 5)
#susie_ann_res <- susie_ann(list(toy_data_1$X, toy_data_2$X), list(toy_data_1$Y, toy_data_2$Y), toy_data_1$A, 2, FALSE, 1, L=1, susie_ann_opt_itr=opt_itr)



# print("Printing toy data:")
# print(toy_data$X)
# print(toy_data$Y)
# print(toy_data$A)
# print(dim(toy_data$Y))
# opt_itr <- 30
# L <- 2
# susie_ann_res <- susie_ann(toy_data$X, toy_data$Y, toy_data$A, rho=0.1, L=L, susie_ann_opt_itr=opt_itr)
# optimized_susie_model <- susie_ann_res$susie_model
# print("Beta:")
# print(optimized_susie_model$beta)
# print("Pi:")
# print(susie_ann_res$final_pi)
# print("Alpha (Lxp):")
# print(optimized_susie_model$alpha)
# print("Mu (Lxp):")
# print(optimized_susie_model$mu)
# print("w:")
# print(susie_ann_res$w)
# print("rho:")
# print(susie_ann_res$all_rho[opt_itr])
# plot(1:opt_itr, susie_ann_res$all_rho, xlab="Alternating optimization iteration", ylab="Rho from optimized w")
# plot(1:opt_itr, susie_ann_res$elbo_values, xlab="Alternating optimization iteration", ylab="Optimized ELBO")
# 
# post_mean_coef = rep(0, ncol(toy_data$X))
# for (l in 1:L){
#   post_mean_coef = post_mean_coef + optimized_susie_model$beta[l] * optimized_susie_model$alpha * optimized_susie_model$mu
# }
# print("Posterior mean regression coefficients:")
# print(t(post_mean_coef))
# plot(toy_data$X %*% t(post_mean_coef), toy_data$Y, xlab="Average model predictions", ylab="Output")