#'Code for calculating the SuSiE-Ann ELBO and for optimizing
#'ELBO with respect to annotation weights w via gradient ascent
#'pi_rho_from_annotation_weights is used to calculate pi from w and A

library(numDeriv)
#' source("SuSiE-Ann/susieR/R/susie_ann.R")

#logsumexp and softmax taken from https://gist.github.com/aufrank/83572
#note: numerical stability is important for softmax
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

#' Obtaining pi from annotation weights w and annotation weights A
#' in basic SuSiE-Ann
pi_from_annotation_weights <- function(A, annotation_weights){
  return (softmax(A %*% annotation_weights))
}

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L) in SuSiE-Ann extended
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights[1:ncol(A)]))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#Calculates Lxp matrix alpha based on p-vector pi and 
#Lxp matrix bf containing Bayes factors for each effect
alpha_from_pi_bf <- function(pi, bf){
  #print('printingproblematic bf in alpha from pi bf')
  p <- length(pi)
  L <- nrow(bf)
  #print(bf)
  #print(pi)
  #print(p)
  alpha <- matrix(1/p, nrow=L, ncol=p)
  #print(alpha)
  for (i in 1:L){
    alpha[i,] = (pi*bf[i,])/sum(pi*bf[i,])
  }
  #print('ending alpha from pi bf')
  return (alpha)
}

#Function executing gradient-based optimization of ELBO w.r.t 
#annotation weights w. Returns updated annotation weights.
gradient_opt_annotation_weights <- function(X_lst,Y_lst,A,is_extended,annotation_weights,susie_fits, batch_size, elbo_opt_steps_per_itr, run_until_convergence=FALSE, step_size=0.05){
  alpha <- list()
  initial_elbo <- cross_locus_elbo(X_lst,Y_lst,A,annotation_weights, susie_fits, is_extended, batch_size)
  if (run_until_convergence==TRUE){
    convergence <- FALSE
    convergence_epsilon <- 0.0005
    itr <- 1
    while (convergence==FALSE || itr <= elbo_opt_steps_per_itr){
    #for (itr in 1:elbo_opt_steps_per_itr){
      #print("current annotation weights w:")
      #print(annotation_weights)
      elbo_gradient <- grad(cross_locus_elbo, x=annotation_weights, X_lst=X_lst, Y_lst=Y_lst, A=A, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size)
      previous_annotation_weights <- annotation_weights
      annotation_weights <- annotation_weights + step_size * elbo_gradient
      itr = itr + 1
      if (norm(annotation_weights - previous_annotation_weights, type="2") < convergence_epsilon){
        convergence = TRUE
      }
    }
  }
  else{
    for (itr in 1:elbo_opt_steps_per_itr){
      elbo_gradient <- grad(cross_locus_elbo, x=annotation_weights, X_lst=X_lst, Y_lst=Y_lst, A=A, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size)
      previous_annotation_weights <- annotation_weights
      annotation_weights <- annotation_weights + step_size * elbo_gradient
      itr = itr + 1
    }
  }
  optimized_elbo <- cross_locus_elbo(X_lst,Y_lst,A,annotation_weights, susie_fits, is_extended, batch_size)
  if (is_extended){
    #Calculating pi and alpha from annotation weights and annotations
    pi_rho <- pi_rho_from_annotation_weights(A, annotation_weights, length(susie_fits[1]$beta))
    pi <- pi_rho$pi
    rho <- pi_rho$rho
  }
  else{
    #Calculating pi from annotation weights and annotations
    pi <- as.vector(softmax(A %*% annotation_weights))
    #print("Pi:")
    #print(pi)
  }
  for (i in 1:length(X_lst)){
    alpha[[i]] <- alpha_from_pi_bf(pi, susie_fits[[i]]$var_lbf)    
  }
  if (is_extended)
    return (list(annot_weights=annotation_weights, pi=pi, rho=rho, alpha=alpha, elbo=optimized_elbo, initial_elbo=initial_elbo))
  else
    return (list(annot_weights=annotation_weights, pi=pi, alpha=alpha, elbo=optimized_elbo, initial_elbo=initial_elbo))
}

#batch_size is number of randomly-selected loci for gradient descent step
cross_locus_elbo <- function(X_lst,Y_lst,A,annotation_weights, susie_fits, is_extended, batch_size){
  elbo <- 0
  selected_locus_indices <- sample(1:length(X_lst), batch_size)
  for (locus_index in selected_locus_indices){
    if (!is_extended){
      elbo = elbo + elbo_basic(X_lst[[locus_index]], Y_lst[[locus_index]], A, annotation_weights, susie_fits[[locus_index]])
      #print("current elbo while being summed across models:")
      #print(elbo)
    }
    else
      elbo = elbo + elbo_extended(X_lst[[locus_index]], Y_lst[[locus_index]], A, annotation_weights, susie_fits[[locus_index]])
  }
  return (elbo)
}

#Calculating ELBO for extended SuSiE-Ann; 
#used in differentiation w.r.t. annotation weights w
elbo_basic <- function(X,Y,A,annotation_weights,susie_fit){
  L <- nrow(susie_fit$var_lbf)
  p <- ncol(susie_fit$var_lbf)
  #Calculating pi and alpha from annotation weights and annotations
  pi <- softmax(A %*% annotation_weights)
  #print('this is where the problem is')
  #print(susie_fit$var_lbf)
  #print(pi)
  #print(annotation_weights)
  #print(A)
  alpha <- alpha_from_pi_bf(pi, susie_fit$var_lbf)
  #Calculating first term of ELBO
  post_avg_b <- rep(0, p)
  for (l in 1:L){
    post_avg_b = post_avg_b + alpha[l,] * susie_fit$mu[l,]
  }
  post_avg_b = matrix(post_avg_b, nrow=length(post_avg_b), ncol=1)
  residual_variance <- susie_fit$sigma2
  #print("posterior average predictions:")
  #print(alpha)
  #print(post_avg_b)
  #print(X %*% post_avg_b)
  #print("true output:")
  #print(Y)
  elbo_first_term <- -1/(2*residual_variance) * norm(Y - X %*% post_avg_b, type="2")^2
  #print("done with elbo first term =")
  #print(elbo_first_term)
  #Calculating second term of ELBO
  scaled_elbo_second_term <- 0
  for (l in 1:L){
    post_avg_b_l = alpha[l,] * susie_fit$mu[l,]
    scaled_elbo_second_term = scaled_elbo_second_term - norm(X %*% post_avg_b_l, type="2")^2
    for (j in 1:p){
      mu2_lj <- susie_fit$mu[l,j]^2
      sigma2_lj <- susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2
      scaled_elbo_second_term = scaled_elbo_second_term + norm(X[,j], type="2")^2 * alpha[l,j] * (mu2_lj + sigma2_lj)
    }
  }
  elbo_second_term <- -1/(2*residual_variance) * scaled_elbo_second_term
  #print('elbo second term=')
  #print(elbo_second_term)
  #Calculating third term of ELBO (note: this term is undefined if pi[j]=alpha_lj=0)
  elbo_third_term <- 0
  for (l in 1:L){
    for (j in 1:p){
      alpha_lj <- alpha[l,j]
      if (susie_fit$V[l] != 0)
        elbo_third_term = elbo_third_term + alpha_lj/2 * (1 + log((susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2)/susie_fit$V[l]) - susie_fit$mu2[l,j]/susie_fit$V[l])
      else
        elbo_third_term = elbo_third_term + alpha_lj/2
      #This conditional technically should not be here according to the formula
      #but avoids the ELBO being undefined
      if (alpha_lj != 0){
        elbo_third_term = elbo_third_term + alpha_lj * log(pi[j]/alpha_lj)
      }
      else
        print("Error: third term of ELBO undefined due to alpha_lj=0")
    }
    #print("elbo third term after one loop")
    #print(elbo_third_term)
  }
  #print('elbo third term:')
  #print(elbo_third_term)
  #print('susie V:')
  #print(susie_fit$V)
  return (elbo_first_term + elbo_second_term + elbo_third_term)
}

#Calculating ELBO for extended SuSiE-Ann; 
#used in differentiation w.r.t. annotation weights w
elbo_extended <- function(X,Y,A,annotation_weights,susie_fit){
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
  residual_variance <- susie_fit$sigma2
  elbo_first_term <- -1/(2*residual_variance) * norm(Y - X %*% post_avg_b, type="2")^2
  #Calculating second term of ELBO
  scaled_elbo_second_term <- 0
  for (l in 1:L){
    for (i in 1:n){
      expected_product <- 0
      for (j in 1:p){
        expected_product = expected_product + susie_fit$beta[l] * alpha[l,j] * X[i,j] * susie_fit$mu[l,j]
      }
      deactivated_contrib <- expected_product^2
      activated_contrib <- 0
      for (j in 1:p){
        sigma2_lj <- susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2
        activated_contrib = activated_contrib +  alpha[l,j] * ((expected_product - X[i,j]*susie_fit$mu[l,j])^2 + (X[i,j]^2 * sigma2_lj))
      }
      scaled_elbo_second_term = scaled_elbo_second_term + susie_fit$beta[l]*activated_contrib + (1-susie_fit$beta[l]) * deactivated_contrib
    }
  }
  elbo_second_term <- -1/(2*residual_variance) * scaled_elbo_second_term
  #Calculating third term of ELBO (note: this term is undefined if pi[j]=alpha_lj=0)
  elbo_third_term <- 0
  #The beta_l and alpha_lj conditionals technically should not be here according to the formula
  #but prevent the ELBO from being undefined
  for (l in 1:L){
    beta_l <- susie_fit$beta[l]
    if (beta_l != 0 && beta_l != 1)
      elbo_third_term = elbo_third_term + (1-beta_l) * log((1-rho)/(1-beta_l)) + beta_l * log(rho/beta_l)
    else
      print("Error: third term of ELBO undefined due to beta_l=0 or 1")
    for (j in 1:p){
      alpha_lj <- alpha[l,j]
      sigma2_lj <- susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2
      if (alpha_lj != 0)
        elbo_third_term = elbo_third_term + beta_l * alpha_lj * (log(pi[j]/alpha_lj) + 1/2 + log(sigma2_lj/susie_fit$V) - susie_fit$mu2[l,j]/susie_fit$V)
      else
        print("Error: third term of ELBO undefined due to alpha_lj=0")
    }
  }
  return (elbo_first_term + elbo_second_term + elbo_third_term)
}

