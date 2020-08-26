#' Optimization of SuSiE-Ann model via alternating calls to SuSiE-Ann IBSS
#' and gradient ascent optimization of ELBO with respect to annotation 
#' weights w

library(rootSolve)
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

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L)
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
  return (list(X=cbind(x1, x2), Y=y, A=A, effect_inclusions=effect_inclusion))
}

#' Function executing alternating optimization of SuSiE-Ann model
#' Currently gradient_opt_annotation_weights is a placeholder function
susie_ann <- function(X,Y, A,
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
                  s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE, 
                  susie_ann_opt_itr=100,
                  elbo_opt_steps_per_itr=100,
                  max_iter=100,tol=1e-3,
                  verbose=FALSE,track_fit=FALSE) {
  A = cbind(A, rep(1, nrow(A)))
  all_annotation_weight_elbo <- rep(0, susie_ann_opt_itr)
  all_iter_rho <- rep(0, susie_ann_opt_itr)
  if (!is.null(annotation_weights)){
    if (pi_rho_from_annotation_weights(A, annotation_weights[1:ncol(A)], L)$rho > 1){
      stop("Specified annotation weights are invalid; rho must be between 0 and 1")
    }
  }
  else
    annotation_weights = init_annotation_weights_from_rho(A, L, rho)
    print("Initialized annotation weights from rho and they equal: ")
    print(annotation_weights)
  for (susie_ann_itr in 1:susie_ann_opt_itr){
    if (susie_ann_itr > 1 | !is.null(s_init)){
      s_init <- susie(X,Y,L=L,rho=rho,
                 scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
                 prior_weights=pi, null_weight=null_weight,
                 standardize=standardize,intercept=intercept,
                 estimate_residual_variance=estimate_residual_variance,
                 estimate_prior_variance = estimate_prior_variance,
                 estimate_prior_method = estimate_prior_method,
                 check_null_threshold=check_null_threshold, prior_tol=prior_tol,
                 residual_variance_upperbound = residual_variance_upperbound,
                 s_init = s_init,coverage=coverage,min_abs_corr=min_abs_corr,
                 compute_univariate_zscore = compute_univariate_zscore,
                 na.rm = na.rm, max_iter=max_iter,tol=tol,
                 verbose=verbose,track_fit=track_fit)
    }
    else{
      print(pi_rho_from_annotation_weights(A, annotation_weights,L)$pi)
      s_init <- susie(X,Y,L=L,rho=rho,
                 scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
                 #prior_weights=pi_rho_from_annotation_weights(A, annotation_weights,L)$pi, 
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
      print(s_init$beta)
    } 
    opt_annot_weights_results <- gradient_opt_annotation_weights(X,Y,A,annotation_weights,s_init,elbo_opt_steps_per_itr)
    annotation_weights <- opt_annot_weights_results$annot_weights
    all_annotation_weight_elbo[susie_ann_itr] = opt_annot_weights_results$elbo
    updated_pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights,L)
    pi <- updated_pi_rho$pi
    rho <- updated_pi_rho$rho
    all_iter_rho[susie_ann_itr] = rho
    s_init$pi = opt_annot_weights_results$pi/sum(opt_annot_weights_results$pi)
    s_init$rho = min(c(opt_annot_weights_results$rho, 1))
    s_init$alpha = opt_annot_weights_results$alpha
    print("Finished one alternating optimization iteration")
  }
  return(list(susie_model=s_init, w=annotation_weights, final_pi=pi, elbo_values=all_annotation_weight_elbo, all_rho=all_iter_rho))
}

# Test Case 1: SuSiE-Ann example based on SuSiE vignette
# set.seed(1)
# n = 3
# p = 2
# b = rep(0,p)
# A = matrix(1:4, nrow=2, ncol=2)
# print("Defined matrix A")
# b[1] = 1
# X = matrix(rnorm(n*p),nrow=n,ncol=p)
# y = X %*% b + rnorm(n)
# print("Defined X and Y; now calling susie_ann")
# susie_ann_res <- susie_ann(X,y,A,rho=0.1, L=2)
# res = susie_ann_res$susie_model
# #coef(res)
# plot(y,predict(res))
# #print(res)

# Test Case 2: Toy data for p=2
toy_data <- generate_toy_data(0.7, 0.7, 0.1, 5, 20)
print("Printing toy data:")
print(toy_data$X)
print(toy_data$Y)
print(toy_data$A)
opt_itr <- 30
L <- 1
susie_ann_res <- susie_ann(toy_data$X, toy_data$Y, toy_data$A, rho=0.1, L=L, susie_ann_opt_itr=opt_itr)
optimized_susie_model <- susie_ann_res$susie_model
print("Beta:")
print(optimized_susie_model$beta)
print("Pi:")
print(susie_ann_res$final_pi)
print("Alpha (Lxp):")
print(optimized_susie_model$alpha)
print("Mu (Lxp):")
print(optimized_susie_model$mu)
print("w:")
print(susie_ann_res$w)
print("rho:")
print(susie_ann_res$all_rho[opt_itr])
plot(1:opt_itr, susie_ann_res$all_rho, xlab="Alternating optimization iteration", ylab="Rho from optimized w")
plot(1:opt_itr, susie_ann_res$elbo_values, xlab="Alternating optimization iteration", ylab="Optimized ELBO")

post_mean_coef = rep(0, ncol(toy_data$X))
for (l in 1:L){
  post_mean_coef = post_mean_coef + optimized_susie_model$beta[l] * optimized_susie_model$alpha * optimized_susie_model$mu
}
print("Posterior mean regression coefficients:")
print(t(post_mean_coef))
plot(toy_data$X %*% t(post_mean_coef), toy_data$Y, xlab="Average model predictions", ylab="Output")