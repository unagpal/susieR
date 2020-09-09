#'Code for calculating the SuSiE-Ann ELBO and for optimizing
#'ELBO with respect to annotation weights w via gradient ascent
#'pi_rho_from_annotation_weights is used to calculate pi from w and A

library(numDeriv)
#' source("SuSiE-Ann/susieR/R/susie_ann.R")

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L)
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights[1:ncol(A)]))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#Function executing gradient-based optimization of ELBO w.r.t 
#annotation weights w. Returns updated annotation weights.
gradient_opt_annotation_weights <- function(X,Y,A,annotation_weights,s,elbo_opt_steps_per_itr, step_size=0.001){
  for (itr in 1:elbo_opt_steps_per_itr){
    #print("current annotation weights w:")
    #print(annotation_weights)
    elbo_gradient <- grad(elbo, x=annotation_weights, X=X, Y=Y, A=A, susie_fit=s)
    #print("elbo gradient:")
    #print(elbo_gradient)
    #print("current elbo:")
    #print(elbo(X, Y, A, annotation_weights, s))
    annotation_weights = annotation_weights + step_size * elbo_gradient
  }
  optimized_elbo <- elbo(X, Y, A, annotation_weights, s)
  #Calculating pi and alpha from annotation weights and annotations
  pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights, length(s$beta))
  pi <- pi_rho$pi
  rho <- pi_rho$rho
  alpha <- alpha_from_pi_bf(pi, s$var_lbf)
  return (list(annot_weights=annotation_weights, pi=pi, rho=rho, alpha=alpha, elbo=optimized_elbo))
}

#Calculating ELBO; used in differentiation w.r.t. annotation 
#weights w
elbo <- function(X,Y,A,annotation_weights,susie_fit){
  L <- nrow(susie_fit$var_lbf)
  p <- ncol(susie_fit$var_lbf)
  #Calculating pi and alpha from annotation weights and annotations
  pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights, length(susie_fit$beta))
  pi <- pi_rho$pi
  rho <- pi_rho$rho
  alpha <- alpha_from_pi_bf(pi, susie_fit$var_lbf)
  #Calculating first term of ELBO
  post_avg_b <- rep(0, p)
  for (l in 1:L){
    post_avg_b = post_avg_b + susie_fit$beta[l] * alpha[l,] * susie_fit$mu[l,]
  }
  post_avg_b = matrix(post_avg_b, nrow=length(post_avg_b), ncol=1)
  elbo_first_term <- norm(Y - X %*% post_avg_b, type="2")^2
  #Calculating second term of ELBO
  elbo_second_term <- 0
  for (l in 1:L){
    for (i in 1:n){
      expected_product <- 0
      for (j in 1:p){
        expected_product = expected_product + susie_fit$beta[l] * alpha[l,j] * X[i,j] * susie_fit$mu[l,j]
      }
      deactivated_contrib <- (1-susie_fit$beta[l]) * expected_product^2
      activated_contrib <- 0
      for (j in 1:p){
        sigma2_lj <- susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2
        activated_contrib = activated_contrib + susie_fit$beta[l] * alpha[l,j] * ((expected_product - X[i,j]*susie_fit$mu[l,j])^2 + (X[i,j]^2 * sigma2_lj))
      }
      elbo_second_term = elbo_second_term + activated_contrib + deactivated_contrib
    }
  }
  #Calculating third term of ELBO (note: this term is undefined if pi[j]=alpha_lj=0)
  elbo_third_term <- 0
  for (l in 1:L){
    beta_l <- susie_fit$beta[l]
    for (j in 1:p){
      alpha_lj <- alpha[l,j]
      #This conditional technically should not be here according to the formula
      #but avoids the ELBO being undefined
      if (alpha_lj != 0){
        elbo_third_term = elbo_third_term + beta_l * alpha_lj * log((rho * pi[j])/(beta_l * alpha_lj))
      }
      else
        print("Error: third term of ELBO undefined")
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