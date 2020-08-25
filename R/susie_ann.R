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

#' Function whose root is taken in initializing annotation weights
difference_from_specified_rho <- function(A,equal_annotation_weight, L, specified_rho){
  annotation_weights <- c(rep(equal_annotation_weight, ncol(A)-1), c(0))
  return (pi_rho_from_annotation_weights(A, annotation_weights, L)$rho - specified_rho)
}

#' Given A and rho, solves for the initial annotation weights
#' where each annotation is weighted equally resulting in the
#' specified value of rho
init_annotation_weights_from_rho <- function(A, L, specified_rho){
  solved_annotation_weight <- uniroot.all(difference_from_specified_rho, interval=c(-100,100), A=A, L=L, specified_rho=specified_rho)[1]
  return (c(rep(solved_annotation_weight, ncol(A)-1), c(0)))
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
  if (annotation_weights != NULL){
    if (pi_rho_from_annotation_weights(A, annotation_weights, L)$rho > 1){
      stop("Specified annotation weights are invalid; rho must be between 0 and 1")
    }
  }
  else
    annotation_weights = init_annotation_weights_from_rho(A, L, rho)
  for (susie_ann_itr in 1:susie_ann_opt_itr){
    if (susie_ann_itr > 1 | s_init != NULL){
      s_init <- susie(X,Y,L=L,rho=rho,
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
    }
    else{
      s_init <- susie(X,Y,L=L,rho=rho,
                 scaled_prior_variance=scaled_prior_variance, residual_variance=residual_variance,
                 prior_weights=prior_weights, null_weight=null_weight,
                 standardize=standardize,intercept=intercept,
                 estimate_residual_variance=estimate_residual_variance,
                 estimate_prior_variance = estimate_prior_variance,
                 estimate_prior_method = estimate_prior_method,
                 check_null_threshold=check_null_threshold, prior_tol=prior_toll,
                 residual_variance_upperbound = residual_variance_upperbound,
                 coverage=coverage,min_abs_corr=min_abs_corr,
                 compute_univariate_zscore = compute_univariate_zscore,
                 na.rm = na.rm, max_iter=max_itr,tol=tol,
                 verbose=verbose,track_fit=track_fit)
    } 
    annotation_weights <- gradient_opt_annotation_weights(X,Y,A,annotation_weights,s_init,elbo_opt_steps_per_itr)
    updated_pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights,L)
    pi <- updated_pi_rho$pi
    rho <- updated_pi_rho$rho
    print("Finished one alternating optimization iteration")
  }
  return(list(susie_model=s, w=annotation_weights))
}

set.seed(1)
n = 10
p = 5
b = rep(0,p)
A = matrix(1:25, nrow=5, ncol=5)
print("Defined matrix A")
b[1:3] = 1
X = matrix(rnorm(n*p),nrow=n,ncol=p)
y = X %*% b + rnorm(n)
print("Defined X and Y; now calling susie_ann")
susie_ann_res <- susie_ann(X,y,A, L=3)
res = susie_ann_res$susie_model
#coef(res)
plot(y,predict(res))
#print(res)