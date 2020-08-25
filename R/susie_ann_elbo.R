#'Code for calculating the SuSiE-Ann ELBO and for optimizing
#'ELBO with respect to annotation weights w via gradient ascent
#'pi_rho_from_annotation_weights is used to calculate pi from w and A

library(numDeriv)
#' source("SuSiE-Ann/susieR/R/susie_ann.R")

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L)
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#Function executing gradient-based optimization of ELBO w.r.t 
#annotation weights w. Returns updated annotation weights.
gradient_opt_annotation_weights <- function(X,Y,A,annotation_weights,s,elbo_opt_steps_per_itr, step_size=0.001){
  for (itr in 1:elbo_opt_steps_per_itr){
    elbo_gradient <- grad(elbo, x=annotation_weights, X=X, Y=Y, A=A, s=s, elbo_opt_steps_per_itr=elbo_opt_steps_per_itr)
    annotation_weights = annotation_weights + step_size * elbo_gradient
  }
  return (annotation_weights)
}

#Calculating ELBO; used in differentiation w.r.t. annotation 
#weights w
elbo <- function(X,Y,A,annotation_weights,s,elbo_opt_steps_per_itr){
  L <- nrow(s$var_lbf)
  p <- ncol(s$var_lbf)
  #Calculating pi and alpha from annotation weights and annotations
  pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights, length(s$beta))
  pi <- pi_rho$pi
  rho <- pi_rho$rho
  alpha <- alpha_from_pi_bf(pi, s$var_lbf)
  #Calculating first term of ELBO
  post_avg_b <- rep(0, p)
  for (l in 1:L){
    post_avg_b = post_avg_b + s$beta[l] * alpha[l,] * s$mu[l,]
  }
  elbo_first_term <- norm(Y - compute_Xb(X, post_avg_b), type="2")^2
  #Calculating second term of ELBO
  elbo_second_term <- 0
  for (l in 1:L){
    for (i in 1:n){
      expected_product <- 0
      for (j in 1:p){
        expected_product = expected_product + s$beta[l] * alpha[l,j] * X[i,j] * s$mu[l,j]
      }
      deactivated_contrib <- (1-s$beta[l]) * expected_product^2
      activated_contrib <- 0
      for (j in 1:p){
        sigma2_lj <- s$mu2[l,j] - s$mu[l,j]^2
        activated_contrib = activated_contrib + s$beta[l] * alpha[l,j] * ((expected_product - X[i,j]*s$mu[l,j])^2 + (X[i,j]^2 * sigma2_lj))
      }
      elbo_second_term = elbo_second_term + activated_contrib + deactivated_contrib
    }
  }
  #Calculating third term of ELBO
  elbo_third_term <- 0
  for (l in 1:L){
    beta_l <- s$beta[l]
    for (j in 1:p){
      alpha_lj <- alpha[l,j]
      elbo_third_term = elbo_third_term + beta_l * alpha_lj * log((rho * pi[j])/(beta_l * alpha_lj))
    }
  }
  return (elbo_first_term + elbo_second_term + elbo_third_term)
}

#Calculates Lxp matrix alpha based on p-vector pi and 
#Lxp matrix bf containing Bayes factors for each effect
alpha_from_pi_bf <- function(pi, bf){
  p <- length(pi)
  L <- nrow(bf)
  alpha <- matrix(1/p, nrow=L, ncol=p)
  for (i in 1:L){
    alpha[i,] = (pi*bf[i,])/sum(pi*bf[i,])
  }
  return (alpha)
}