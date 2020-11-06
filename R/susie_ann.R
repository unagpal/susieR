#' Optimization of SuSiE-Ann model via alternating calls to SuSiE-Ann IBSS
#' and gradient ascent optimization of ELBO with respect to annotation 
#' weights w

library(rootSolve)
library(faux)
library(LaplacesDemon)
library(clusterGeneration)
#source("SuSiE-Ann/susieR/R/set_X_attributes.R")
#source("SuSiE-Ann/susieR/R/initialize.R")
#source("SuSiE-Ann/susieR/R/update_each_effect.R")
#source("SuSiE-Ann/susieR/R/estimate_residual_variance.R")
#source("SuSiE-Ann/susieR/R/susie_utils.R")
source("SuSiE-Ann/susieR/R/susie.R")
source("SuSiE-Ann/susieR/R/susie_ann_elbo.R")

#' Sigmoid function
sigmoid <- function(vec){
  return (1/(1+exp(-vec)))
}

#' Obtaining pi from annotation weights w and annotation weights A
#' in basic SuSiE-Ann
pi_from_annotation_weights <- function(A, annotation_weights){
  return (exp(A %*% annotation_weights)/sum(exp(A %*% annotation_weights)))
}

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L) for extended SuSiE-Ann
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights[1:ncol(A)]))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#' Function whose root is taken in initializing annotation weights
difference_from_specified_rho <- function(A,equal_annotation_weight, L, specified_rho){
  all_differences <- rep(0, length(equal_annotation_weight))
  for (i in 1:length(equal_annotation_weight)){
    annotation_weights <- c(rep(equal_annotation_weight[i], ncol(A)-1), c(0))
    all_differences[i] = pi_rho_from_annotation_weights(A, annotation_weights, L)$rho - specified_rho
  }
  return (all_differences)
}

#' Given A and rho, solves for the initial annotation weights
#' where each annotation is weighted equally resulting in the
#' specified value of rho
init_annotation_weights_from_rho <- function(A, L, specified_rho){
  solved_annotation_weight <- uniroot.all(difference_from_specified_rho, interval=c(-100,100), A=A, L=L, specified_rho=specified_rho)[1]
  return (c(rep(solved_annotation_weight, ncol(A)-1), c(0)))
}

#' Taken from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/313138
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

#Function for generating an easy example with p=2
generate_toy_data <- function(beta, corr, noise_ratio, annot_signal_ratio, n){
  A <- matrix(c(annot_signal_ratio,1,1,annot_signal_ratio), nrow=2, ncol=2)
  x1 <- rnorm(n)
  x1 = x1 - mean(x1)
  x2 <- complement(x1, corr)
  x2 = x2 - mean(x2)
  if (var(x1) > var(x2))
    x2 = x2 * sqrt(var(x1)/var(x2))
  effect_inclusion <- rbinom(n, beta, size=1)
  y <- x1 * effect_inclusion + noise_ratio * rnorm(n)
  y = matrix(y, nrow=length(y), ncol=1)
  return (list(X=cbind(x1, x2), Y=y, A=A, effect_inclusions=effect_inclusion))
}

#Correlation absolute values ~ Beta(2,5)
#Correlation signs ~ Bernoulli(0.5)
generate_synthetic_data_basic <- function(pve, n, p, L, num_loci, num_annotations){
  output = list()
  w <- rnorm(num_annotations, mean = 0, sd = 1.5)
  A_lst <- list()
  pi <- list()
  for (l in 1:num_loci){
    A_l = matrix(0, nrow=p, ncol=num_annotations)
    for (i in 1:p)
      A_l[i,] <- rbinom(num_annotations, 1, 0.1)
    A_lst[[l]] <- A_l
    pi[[l]] <- softmax(A_lst[[l]] %*% w)
  }
  susie_ann_X <- list()
  susie_ann_Y <- list()
  all_selected_variables <- list()
  for (k in 1:num_loci){
    selected_locus_variables <- matrix(0, nrow=L, ncol=p)
    random_cov_matrix = genPositiveDefMat(dim=p, covMethod="eigen")$Sigma
    #random_cov_matrix <- rWishart(1,p,diag(p))
    random_corr_matrix <- matrix(0, nrow=p, ncol=p)
    for (i in 1:p){
      for (j in 1:p){
        random_corr_matrix[i,j] <- random_cov_matrix[i,j]/(sqrt(random_cov_matrix[i,i])*sqrt(random_cov_matrix[j,j]))
      }
    }
    random_corr_matrix = round(random_corr_matrix,3)
    X <- rnorm_multi(n = n,
                       mu = 0,
                       sd = 1,
                       r = random_corr_matrix,
                       varnames = c(1:p),
                       empirical = FALSE)
    b <- as.matrix(rep(0, p))
    for (l in 1:L){
      selected_variable <- rmultinom(1,1,pi[[k]])
      selected_var_ind <- match(1, selected_variable)
      selected_locus_variables[l,] <- selected_variable
      b[selected_var_ind] = rnorm(1)
    }
    all_selected_variables[[k]] <- selected_locus_variables
    var_xb <- var(as.matrix(X) %*% b)
    sigma <- sqrt((1/pve - 1)*var_xb)
    Y = as.matrix(X) %*% b + rnorm(n, sd=sigma)
    susie_ann_X[[k]] <- as.matrix(X)
    susie_ann_Y[[k]] <- as.matrix(Y)
    #print(Y)
    output[[k]] <- list(X=as.matrix(X), Y=as.matrix(Y), b=b, sigma=sigma)
  }
  return (list(susie_ann_X=susie_ann_X, susie_ann_Y=susie_ann_Y, A_lst=A_lst, w=w, pi=pi, output=output, selected_variables=all_selected_variables))
}

#' Function executing alternating optimization of SuSiE-Ann model
#' Currently gradient_opt_annotation_weights is a placeholder function
susie_ann <- function(X_lst,Y_lst, A_lst, num_loci, extended_model, batch_size,
                  annotation_weights=NULL, rho=0.2, 
                  L = min(10,ncol(X)),
                  scaled_prior_variance=0.2, residual_variance=NULL,
                  prior_weights=NULL, null_weight=NULL,
                  standardize=TRUE,intercept=TRUE,
                  estimate_residual_variance=TRUE,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim","EM","simple"),
                  check_null_threshold=0, prior_tol=1E-9,
                  residual_variance_upperbound = Inf,
                  s_init_lst = NULL,coverage=0.95,min_abs_corr=0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE, 
                  susie_ann_opt_itr=100,
                  elbo_opt_steps_per_itr=100,
                  max_iter=100,tol=1e-3,
                  verbose=FALSE,track_fit=FALSE) {
  num_loci <- length(X_lst)
  if (extended_model){
    for (l in 1:num_loci)
      A_lst[[l]] = cbind(A_lst[[l]], rep(1, nrow(A_lst[[l]])))
  }
  else
    all_iter_rho <- NULL
  all_initial_annotation_weight_elbo <- rep(0, susie_ann_opt_itr)
  all_opt_annotation_weight_elbo <-  rep(0, susie_ann_opt_itr)
  if (extended_model){
    all_iter_rho <- rep(0, susie_ann_opt_itr)
    if (!is.null(annotation_weights)){
      if (pi_rho_from_annotation_weights(A, annotation_weights[1:ncol(A)], L)$rho > 1){
        stop("Specified annotation weights are invalid; rho must be between 0 and 1")
      }
    }
    else{
      annotation_weights = init_annotation_weights_from_rho(A, L, rho)
      print("Initialized annotation weights from rho and they equal: ")
      print(annotation_weights)
    }
  }
  else{
    if (is.null(annotation_weights))
      annotation_weights = rep(0, ncol(A_lst[[1]]))
  }
  all_itr_pi <- list()
  all_itr_alpha <- list()
  for (susie_ann_itr in 1:susie_ann_opt_itr){
    if (!is.null(s_init_lst)){
      for (l in 1:num_loci){
        s_init_lst[[l]] <- susie(X_lst[[l]],Y_lst[[l]],extended_model, L=L,rho=rho,
                 scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
                 prior_weights=pi[[l]], null_weight=null_weight,
                 standardize=standardize,intercept=intercept,
                 estimate_residual_variance=estimate_residual_variance,
                 estimate_prior_variance = estimate_prior_variance,
                 estimate_prior_method = estimate_prior_method,
                 check_null_threshold=check_null_threshold, prior_tol=prior_tol,
                 residual_variance_upperbound = residual_variance_upperbound,
                 s_init = s_init_lst[[l]],coverage=coverage,min_abs_corr=min_abs_corr,
                 compute_univariate_zscore = compute_univariate_zscore,
                 na.rm = na.rm, max_iter=max_iter,tol=tol,
                 verbose=verbose,track_fit=track_fit)
      }
    }
    else{
      print("Obtaining pi from annotation weights:")
      if (extended_model)
        print(pi_rho_from_annotation_weights(A[[1]], annotation_weights,L)$pi)
      else
        print(pi_from_annotation_weights(A[[1]], annotation_weights))
      s_init_lst <- list()
      
      for (l in 1:num_loci){
        print(residual_variance)
        s_init_lst[[l]] <- susie(X_lst[[l]],Y_lst[[l]],extended_model, L=L,rho=rho,
                   #scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
                   prior_weights=prior_weights,
                   null_weight=null_weight,
                   standardize=standardize,intercept=intercept,
                   estimate_residual_variance=estimate_residual_variance,
                   estimate_prior_variance = estimate_prior_variance,
                   estimate_prior_method = estimate_prior_method,
                   check_null_threshold=check_null_threshold, prior_tol=prior_tol,
                   residual_variance_upperbound = residual_variance_upperbound,
                   coverage=coverage,min_abs_corr=min_abs_corr,
                   compute_univariate_zscore = compute_univariate_zscore,
                   na.rm = na.rm, max_iter=max_iter,tol=tol,
                   verbose=verbose,track_fit=track_fit)
      }
      print("Done initializing SuSiE in first iteration")
      #print(s_init$beta)
    } 
    #for (l in 1:num_loci){
      #print("Optimized SuSiE Model for Locus Number:")
      #print(l)
      #print(s_init_lst[[1]])
    #}
    print("starting to optimize annotation weights")
    if (susie_ann_itr != susie_ann_opt_itr)
      opt_annot_weights_results <- gradient_opt_annotation_weights(X_lst,Y_lst,A_lst,extended_model, annotation_weights,s_init_lst,batch_size,elbo_opt_steps_per_itr, run_until_convergence=TRUE)
    else
      opt_annot_weights_results <- gradient_opt_annotation_weights(X_lst,Y_lst,A_lst,extended_model, annotation_weights,s_init_lst,batch_size,elbo_opt_steps_per_itr, run_until_convergence=TRUE)
    print("finished optimizing annotation weights")
    all_initial_annotation_weight_elbo[susie_ann_itr] = opt_annot_weights_results$initial_elbo
    all_opt_annotation_weight_elbo[susie_ann_itr] = opt_annot_weights_results$elbo
    annotation_weights <- opt_annot_weights_results$annot_weights
    print("Returned annotation weights:")
    print(annotation_weights)
    pi <- opt_annot_weights_results$pi
    print("Returned pi:")
    print(pi)
    if (extended_model){
      rho <- updated_pi_rho$rho
      all_iter_rho[susie_ann_itr] = rho
      for (i in 1:num_loci){
        s_init_lst[[i]]$rho=rho
      }
    }
    all_itr_pi[[susie_ann_itr]] <- pi[[1]]
    #s_init$pi = opt_annot_weights_results$pi/sum(opt_annot_weights_results$pi)
    #s_init$rho = min(c(opt_annot_weights_results$rho, 1))
    #print("pi before gradient update:")
    #print(s_init$pi)
    for (i in 1:num_loci){
      s_init_lst[[i]]$pi=pi[[i]]
      s_init_lst[[i]]$alpha = opt_annot_weights_results$alpha[[i]]
      #print("alpha that is being set to alpha for susie fit:")
      #print(opt_annot_weights_results$alpha[[i]])
    }
    print("dim(alpha) for susie model 1 during susie-ann:")
    print(dim(s_init_lst[[1]]$alpha))
    all_itr_alpha[[l]] <- s_init_lst[[1]]$alpha[[1]]
    #print("all alphas:")
    #print(opt_annot_weights_results$alpha)
    #print("pi after gradient update:")
    #print('alpha after gradient update:')
    #print(s_init$alpha)
    print("Initial ELBO:")
    print(opt_annot_weights_results$initial_elbo)
    print("ELBO:")
    print(opt_annot_weights_results$elbo)
    print("One alternating optimization iteration complete")
    #The code below can be uncommented to compare the ELBO of the 'true' parameters
    #to that arrived at by SuSiE-Ann in synthetic experiments.
    # if (susie_ann_itr == susie_ann_opt_itr){
    #   s_final_lst <- list()
    #   for (l in 1:num_loci){
    #     s_final_lst[[l]] <- susie(X_lst[[l]],Y_lst[[l]],extended_model, L=L,rho=rho,
    #                               scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
    #                               prior_weights=synth_data$pi, null_weight=null_weight,
    #                               standardize=standardize,intercept=intercept,
    #                               estimate_residual_variance=estimate_residual_variance,
    #                               estimate_prior_variance = estimate_prior_variance,
    #                               estimate_prior_method = estimate_prior_method,
    #                               check_null_threshold=check_null_threshold, prior_tol=prior_tol,
    #                               residual_variance_upperbound = residual_variance_upperbound,
    #                               s_init = s_init_lst[[l]],coverage=coverage,min_abs_corr=min_abs_corr,
    #                               compute_univariate_zscore = compute_univariate_zscore,
    #                               na.rm = na.rm, max_iter=max_iter,tol=tol,
    #                               verbose=verbose,track_fit=track_fit)
    #   }
    #   correct_answer_opt_annot_weights_results <- gradient_opt_annotation_weights(X_lst,Y_lst,A_lst,extended_model, synth_data$w,s_final_lst,batch_size,0)
    #   print("Initial ELBO with right annot. weights")
    #   print(correct_answer_opt_annot_weights_results$initial_elbo)
    #   print("Final ELBO with right annot. weights")
    #   print(correct_answer_opt_annot_weights_results$elbo)
    # }
  }
  return(list(susie_models=s_init_lst, w=annotation_weights, final_pi=pi, all_iter_pi=all_itr_pi, all_iter_alpha=all_itr_alpha, initial_elbo_values=all_initial_annotation_weight_elbo, opt_elbo_values=all_opt_annotation_weight_elbo, all_rho=all_iter_rho))
}

analyze_susie_ann_fit <- function(res, synth_data) {
  #First plotting “initial” and “final” (i.e. before and after optimizing ELBO w.r.t w) ELBO as a function of iterations
  alternating_opt_itr <- length(res$initial_elbo_values)
  plot(1:alternating_opt_itr, res$initial_elbo_values, xlab="Alternating optimization iteration", ylab="SuSiE-Ann ELBO before optimizing w")
  plot(1:alternating_opt_itr, res$opt_elbo_values, xlab="Alternating optimization iteration", ylab="SuSiE-Ann ELBO after optimizing w")
  #Next obtaining predicted and actual pi for all loci
  actual_pi <- synth_data$pi
  pred_pi <- res$final_pi
  num_loci <- length(actual_pi)
  all_corr_pi <- rep(0, num_loci)
  all_pval_pi <- rep(0, num_loci)
  all_pi_mae <- rep(0, num_loci)
  all_kl_pi <- rep(0, num_loci)
  selected_loci_graphs <- sample(1:num_loci, 2)
  for (l in 1:num_loci){
    all_corr_pi[l] <- cor(pred_pi[[l]], actual_pi[[l]])
    fit = lm(actual_pi[[l]] ~ pred_pi[[l]])
    all_pval_pi[l] <- summary(fit)$coefficients[2,4]
    if (l %in% selected_loci_graphs){
      plot(pred_pi[[l]], actual_pi[[l]], xlab="Predicted pi", ylab="Actual pi")
    }
    all_pi_mae[l] <- mean(abs(pred_pi[[l]] - actual_pi[[l]]))
    all_kl_pi[l] <- KLD(pred_pi[[l]], actual_pi[[l]])$sum.KLD.px.py
  }
  print("All Locus R coefficients for Actual vs. Predicted Pi")
  print(all_corr_pi)
  hist(all_corr_pi, xlab="Locus R coefficients for Actual vs. Predicted Pi")
  print("All Locus P values for Actual vs. Predicted Pi Linear Regression")
  print(all_pval_pi)
  hist(all_pval_pi, xlab="Locus P-Values for Actual vs. Predicted Pi Linear Regression")
  print("All Locus MAE for Actual vs. Predicted Pi Linear Regression")
  print(all_pi_mae)
  hist(all_pi_mae, xlab="All Locus MAE for Actual vs. Predicted Pi Linear Regression")
  print("All locus D_{KL} (Predicted Pi, Actual Pi):")
  print(all_kl_pi)
  hist(all_kl_pi, xlab="All locus D_{KL} (Predicted Pi, Actual Pi)")
  #Getting predicted and actual annotation weights
  w_pred <- res$w
  w_actual <- synth_data$w
  plot(w_pred, w_actual, xlab="Predicted annotation weight", ylab="Actual annotation weight")
  print("Correlation coefficient between predicted and actual w:")
  print(cor(w_pred, w_actual))
  fit = lm(w_actual ~ w_pred)
  print("P value for linear regression between predicted and actual w:")
  print(summary(fit)$coefficients[2,4])
  #Next, comparing actual alpha (i.e. which variables are effects) with predicted alpha
  susie_models <- res$susie_models
  #all_pred_alpha[[i]] is nxp
  all_pred_alpha <- list()
  for (s in 1:length(susie_models)){
    all_pred_alpha[[s]] <- susie_models[[s]]$alpha
    print(dim(all_pred_alpha[[s]]))
  }
  print(dim(synth_data$all_selected_variables[[s]]))
  #synth_data$all_selected_variables[[s]]
}
