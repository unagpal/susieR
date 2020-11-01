#' @title update each effect once
#' @param X an n by p matrix of regressor variables
#' @param Y an n vector of response variable
#' @param s a SuSiE fit
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
#' @param check_null_threshold float a threshold on the log scale to compare likelihood between current estimate and zero the null
#' @keywords internal
source("SuSiE-Ann/susieR/R/single_effect_regression.R")
source("SuSiE-Ann/susieR/R/sparse_multiplication.R")
source("SuSiE-Ann/susieR/R/elbo.R")
update_each_effect <- function (X, Y, s, estimate_prior_variance=FALSE,
                                estimate_prior_method="optim",check_null_threshold) {
  if(estimate_prior_variance==FALSE) estimate_prior_method="none"

  # Repeat for each effect to update
  L = nrow(s$alpha)
  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      if (s$extended_model)
        s$Xr = s$Xr - compute_Xb(X, (s$beta[l] * s$alpha[l,] * s$mu[l,]))
      else
        s$Xr = s$Xr - compute_Xb(X, (s$alpha[l,] * s$mu[l,]))
      #compute residuals
      R = Y - s$Xr
      res <- single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                      estimate_prior_method,check_null_threshold)
      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$lbf[l] <- res$lbf_model
      s$var_lbf[l,] <- exp(res$lbf)
      s$KL[l] <- -res$loglik + SER_posterior_e_loglik(X,R,s$sigma2,res$alpha*res$mu,res$alpha*res$mu2)
      
      
      if (s$extended_model){
        s$beta[l] <- s$rho*activated_effect_susie_ann_likelihood(X, Y, s, l) / ((1-s$rho)*deactivated_effect_susie_ann_likelihood(X, Y, s, l) + s$rho * activated_effect_susie_ann_likelihood(X, Y, s))
        s$Xr <- s$Xr + compute_Xb(X, (s$beta[l] * s$alpha[l,] * s$mu[l,]))
      }
      else
        s$Xr <- s$Xr + compute_Xb(X, (s$alpha[l,] * s$mu[l,]))
    }
  }

  return(s)
}

#' @title Computes p(y | X, b, beta, prior_variance, residual_variance) approximately where effect l is always activated, that is, beta[l]=1
activated_effect_susie_ann_likelihood <- function(X, Y, s, l){
  activated_effect_susie_model <- s
  activated_effect_susie_model$beta[l] <- 1
  return(susie_ann_likelihood(X, Y, activated_effect_susie_model))
}

#' @title Computes p(y | X, b, beta, prior_variance, residual_variance) approximately where effect l is always deactivated, that is, beta[l]=0
deactivated_effect_susie_ann_likelihood <- function(X, Y, s, l){
  deactivated_effect_susie_model <- s
  deactivated_effect_susie_model$beta[l] <- 0
  return(susie_ann_likelihood(X, Y, deactivated_effect_susie_model))
}

#' @title Computes p(y | X, b, beta, prior_variance, residual_variance) approximately by drawing from the model posterior
susie_ann_likelihood <- function(X, Y, s, posterior_draws=500){
  p <- ncol(X)
  L <- nrow(s$alpha)
  all_likelihoods <- rep(0, posterior_draws)
  for (post_draw in 1:posterior_draws){
    post_sample_b <- rep(0, p)
    sampled_beta <- rbinom(n=L, size=1, prob=min(c(s$beta, 1.0)))
    for (l in 1:L){
      sampled_gamma <- rmultinom(1, size=1, prob=s$alpha[l,])
      sampled_coef <- rnorm(L, mean=s$mu[l,], sqrt(s$mu2[l,] - s$mu[l,]^2))
      post_sample_b = post_sample_b + s$mu[l,] * sampled_gamma * sampled_beta[l]
    }
    residuals <- Y - compute_Xb(X,post_sample_b)
    all_likelihoods[post_draw] = exp(sum(dnorm(residuals,mean=0, sd=sqrt(s$sigma2), log=TRUE)))     
  }
  return(mean(all_likelihoods))
}
