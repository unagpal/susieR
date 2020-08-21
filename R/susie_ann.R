#' Optimization of SuSiE-Ann model via alternating calls to SuSiE-Ann IBSS
#' and gradient ascent optimization of ELBO with respect to annotation 
#' weights w

source("SuSiE-Ann/susieR/R/set_X_attributes.R")
source("SuSiE-Ann/susieR/R/initialize.R")
source("SuSiE-Ann/susieR/R/update_each_effect.R")
source("SuSiE-Ann/susieR/R/estimate_residual_variance.R")
source("SuSiE-Ann/susieR/R/susie_utils.R")
source("SuSiE-Ann/susieR/R/susie.R")
source("SuSiE-Ann/susieR/R/susie_ann_elbo.R")

#' Sigmoid function
sigmoid <- function(vec){
  return (1/(1+exp(-vec)))
}

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L)
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#' Function executing alternating optimization of SuSiE-Ann model
#' Currently gradient_opt_annotation_weights is a placeholder function
susie_ann <- function(X,Y, A, annotation_weights=rep(1/ncol(A),nrow(A)),
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
  init_pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights, L)
  pi <- init_pi_rho$pi
  rho <- init_pi_rho$rho
  for (susie_ann_itr in 1:susie_ann_itr){
    if (susie_ann_itr > 1){
      s_init = s
    }
    s <- susie(X,Y,L=L,rho=rho,
               scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
               prior_weights=pi, null_weight=null_weight,
               standardize=standardize,intercept=intercept,
               estimate_residual_variance=estimate_residual_variance,
               estimate_prior_variance = estimate_prior_variance,
               estimate_prior_method = estimate_prior_method,
               check_null_threshold=check_null_threshold, prior_tol=prior_toll,
               residual_variance_upperbound = residual_variance_upperbound,
               s_init = s_init,coverage=coverage,min_abs_corr=min_abs_corr,
               compute_univariate_zscore = compute_univariate_zscore,
               na.rm = na.rm, max_iter=max_itr,tol=tol,
               verbose=verbose,track_fit=track_fit)
    annotation_weights <- gradient_opt_annotation_weights(X,Y,A,annotation_weights,s,elbo_opt_steps_per_itr)
    updated_pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights,L)
    pi <- updated_pi_rho$pi
    rho <- updated_pi_rho$rho
  }
  return(list(susie_model=s, w=annotation_weights))
}