#'Code for calculating the SuSiE-Ann ELBO and for optimizing
#'ELBO with respect to annotation weights w via gradient ascent or L-BFGS-B

library(qrnn)
library(numDeriv)
library(wordspace)
library(philentropy)
library(SimDesign)
library(matrixStats)
#' source("SuSiE-Ann/susieR/R/susie_ann.R")

#logsumexp and softmax taken from https://gist.github.com/aufrank/83572
#note: numerical stability is important for softmax
logsumexp <- function (x) {
  y = max(x)
  return (y + log(sum(exp(x - y))))
}

#numerically stable softmax function
softmax <- function (x) {
  return (exp(x - logsumexp(x)))
}

#' Obtaining pi from annotation weights w and annotation matrix A
#' in basic SuSiE-Ann/Proposal 0
pi_from_annotation_weights <- function(A, annotation_weights){
  return (softmax(A %*% annotation_weights))
}

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L) in SuSiE-Ann extended/Proposal 2
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights[1:ncol(A)]))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#Calculates Lxp matrix alpha based on p-vector pi and 
#Lxp matrix bf containing Bayes factors for each effect
alpha_from_pi_bf <- function(pi, bf){
  p <- length(pi)
  L <- nrow(bf)
  alpha <- matrix(1/p, nrow=L, ncol=p)
  for (i in 1:L){
    if (sum(pi*bf[i,]) > 0)
      alpha[i,] = (pi*bf[i,])/sum(pi*bf[i,])
    else
      alpha[i,] = 1/p
  }
  return (alpha)
}

#Key function executing gradient-based optimization of ELBO w.r.t 
#annotation weights w. Returns updated annotation weights. Parameters:
#batch_size: batch size of loci
#optim_method: optimization method, "SGD" or "L-BFGS-B"
gradient_opt_annotation_weights <- function(X_lst,Y_lst,A_lst,is_extended, optim_method, annotation_weights,susie_fits, batch_size, elbo_opt_steps_per_itr=100, run_until_convergence=TRUE, step_size=0.1){
  print("Starting gradient optimization of w")
  alpha <- list()
  #Calculating ELBO prior to optimization of w
  initial_elbo <- cross_locus_elbo(annotation_weights,X_lst,Y_lst,A_lst,susie_fits, is_extended, batch_size)
  num_loci <- length(X_lst)
  #Stochastic gradient based ELBO ascent w.r.t. w
  if (optim_method=="SGD"){
    #Running to convergence, typically done in last alternating optimization iteration
    if (run_until_convergence==TRUE){
      convergence <- FALSE
      convergence_epsilon <- 0.005
      itr <- 1
      while (convergence==FALSE && itr<10*elbo_opt_steps_per_itr){
        print("Running until convergence, starting new itr")
        
        elbo_gradient <- grad(cross_locus_elbo, x=annotation_weights, X_lst=X_lst, Y_lst=Y_lst, A_lst=A_lst, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size)
        previous_annotation_weights <- annotation_weights
        
        #Can toggle between moving in normalized gradient direction or unnormalized gradient direction
        annotation_weights <- annotation_weights +  elbo_gradient/(100*norm(elbo_gradient, type="2"))
        #annotation_weights <- annotation_weights + step_size * elbo_gradient
        
        #elbo_improvement <- cross_locus_elbo(X_lst,Y_lst,A_lst,annotation_weights, susie_fits, is_extended, batch_size) - cross_locus_elbo(X_lst,Y_lst,A_lst,previous_annotation_weights, susie_fits, is_extended, batch_size)
        
        #Seeing/printing improvement in ELBO with each gradient step
        elbo_improvement <- cross_locus_elbo(annotation_weights,X_lst,Y_lst,A_lst, susie_fits, is_extended, batch_size) - cross_locus_elbo(previous_annotation_weights,X_lst,Y_lst,A_lst, susie_fits, is_extended, batch_size)
        print("Improvement in ELBO (running till max convergence itr or convergence):")
        print(elbo_improvement)
        itr = itr + 1
        
        previous_annotation_weight_norm <- norm(previous_annotation_weights, type="2")
        #Can toggle between two convergence criteria: 1) convergence of ELBO 2) convergence of w vector
        if (elbo_improvement < convergence_epsilon){
        #if (norm(annotation_weights - previous_annotation_weights, type="2") < convergence_epsilon*previous_annotation_weight_norm){
          convergence = TRUE
        }
      }
    }
    #Optimizing for set number of gradient steps
    else{
      for (itr in 1:elbo_opt_steps_per_itr){
        elbo_gradient <- grad(cross_locus_elbo, x=annotation_weights, X_lst=X_lst, Y_lst=Y_lst, A_lst=A_lst, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size)
        previous_annotation_weights <- annotation_weights +  elbo_gradient/(10*norm(elbo_gradient, type="2"))
        annotation_weights <- annotation_weights + step_size * elbo_gradient
        itr = itr + 1
      }
    }
  }
  #Optimizing via L-BFGS-B
  else if (optim_method=="L-BFGS-B"){
    optim_res <- optim(par=annotation_weights, fn=cross_locus_elbo, gr = NULL, X_lst=X_lst,Y_lst=Y_lst,A_lst=A_lst, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size, method=optim_method, lower = -Inf, upper = Inf,
          control = list(fnscale=-1))
    annotation_weights = optim_res$par
    print("Optimized ELBO value:")
    print(optim_res$value)
  }
  
  #Adam optimization not yet implemented as I am looking for an R implementation
  #enabling optimization of a user-specified non-neural-network function
  #else if (optim_method=="ADAM"){
    #adam_res <- adam(f=(-1) * cross_locus_elbo, p=annotation_weights, X_lst=X_lst,Y_lst=Y_lst,A_lst=A_lst, susie_fits=susie_fits, is_extended=is_extended, batch_size=batch_size, print.level=0)
    #annotation_weihgts = adam_res$estimate
    #print("Optimized ELBO value:")
    #print(adam_res$minimum)
  #}
  optimized_elbo <- cross_locus_elbo(annotation_weights,X_lst,Y_lst,A_lst, susie_fits, is_extended, batch_size)
  pi <- list()
  #Updating pi and rho based on optimized w for each locus, for extended SuSiE-Ann/Proposal 2
  if (is_extended){
    #Calculating pi and alpha from annotation weights and annotations
    pi_rho <- pi_rho_from_annotation_weights(A_lst, annotation_weights, length(susie_fits[1]$beta))
    pi <- pi_rho$pi
    rho <- pi_rho$rho
  }
  #Updating pi based on optimized w for each locus, for basic SuSiE-Ann/Proposal 0
  else{
    for (l in 1:num_loci){
      #Calculating pi from annotation weights and annotations
      pi[[l]] <- as.vector(softmax(A_lst[[l]] %*% annotation_weights))
    }
  }
  #Calculating alpha based on optimized w for each locus
  for (l in 1:num_loci){
    log_pi <- as.vector(log(pi[[l]]))
    log_pi_bf <- log_pi + susie_fits[[l]]$var_lbf
    log_alpha <- log_pi_bf - rowLogSumExps(log_pi_bf)
    alpha[[l]] <- exp(log_alpha)
  }
  print("Returned alpha from SuSiE-Ann:")
  print(alpha)
  if (is_extended)
    return (list(annot_weights=annotation_weights, pi=pi, rho=rho, alpha=alpha, elbo=optimized_elbo, initial_elbo=initial_elbo))
  else
    return (list(annot_weights=annotation_weights, pi=pi, alpha=alpha, elbo=optimized_elbo, initial_elbo=initial_elbo))
}

#The SuSiE-Ann ELBO function, based on random subsampling of the loci, that is optimized. Parameters:
#batch_size: number of randomly-selected loci for each gradient descent step
#is_extended: True indicates SuSiE-Ann extended/Proposal 2, False indicates basic SuSiE-Ann/Proposal 0
cross_locus_elbo <- function(annotation_weights,X_lst,Y_lst,A_lst, susie_fits, is_extended, batch_size){
  elbo <- 0
  #Sampling batch_size loci for gradient-based optimization of w
  selected_locus_indices <- sample(1:length(X_lst), batch_size)
  for (locus_index in selected_locus_indices){
    if (!is_extended){
      elbo = elbo + elbo_basic(X_lst[[locus_index]], Y_lst[[locus_index]], A_lst[[locus_index]], annotation_weights, susie_fits[[locus_index]])
    }
    else
      elbo = elbo + elbo_extended(X_lst[[locus_index]], Y_lst[[locus_index]], A_lst[[locus_index]], annotation_weights, susie_fits[[locus_index]])
  }
  return (elbo)
}

#Non-vectorized ELBO for basic SuSiE-Ann (one locus)
#used in differentiation w.r.t. annotation weights w
# elbo_basic <- function(X,Y,A,annotation_weights,susie_fit){
#   L <- nrow(susie_fit$var_lbf)
#   p <- ncol(susie_fit$var_lbf)
#   #Calculating pi and alpha from annotation weights and annotations
#   pi <- softmax(A %*% annotation_weights)
#   #print('this is where the problem is')
#   #print(susie_fit$var_lbf)
#   #print(pi)
#   #print(annotation_weights)
#   #print(A)
#   alpha <- alpha_from_pi_bf(pi, susie_fit$var_lbf)
#   #Calculating first term of ELBO
#   post_avg_b <- rep(0, p)
#   for (l in 1:L){
#     post_avg_b = post_avg_b + alpha[l,] * susie_fit$mu[l,]
#   }
#   post_avg_b = matrix(post_avg_b, nrow=length(post_avg_b), ncol=1)
#   residual_variance <- susie_fit$sigma2
#   #print("posterior average predictions:")
#   #print(alpha)
#   #print(post_avg_b)
#   #print(X %*% post_avg_b)
#   #print("true output:")
#   #print(Y)
#   elbo_first_term <- -1/(2*residual_variance) * norm(Y - X %*% post_avg_b, type="2")^2
#   #print("done with elbo first term =")
#   #print(elbo_first_term)
#   #Calculating second term of ELBO
#   scaled_elbo_second_term <- 0
#   for (l in 1:L){
#     post_avg_b_l = alpha[l,] * susie_fit$mu[l,]
#     scaled_elbo_second_term = scaled_elbo_second_term - norm(X %*% post_avg_b_l, type="2")^2
#     for (j in 1:p){
#       mu2_lj <- susie_fit$mu[l,j]^2
#       sigma2_lj <- susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2
#       scaled_elbo_second_term = scaled_elbo_second_term + norm(X[,j], type="2")^2 * alpha[l,j] * (mu2_lj + sigma2_lj)
#     }
#   }
#   elbo_second_term <- -1/(2*residual_variance) * scaled_elbo_second_term
#   #print('elbo second term=')
#   #print(elbo_second_term)
#   #Calculating third term of ELBO (note: this term is undefined if pi[j]=alpha_lj=0)
#   elbo_third_term <- 0
#   for (l in 1:L){
#     for (j in 1:p){
#       alpha_lj <- alpha[l,j]
#       if (susie_fit$V[l] != 0)
#         elbo_third_term = elbo_third_term + alpha_lj/2 * (1 + log((susie_fit$mu2[l,j] - susie_fit$mu[l,j]^2)/susie_fit$V[l]) - susie_fit$mu2[l,j]/susie_fit$V[l])
#       else
#         elbo_third_term = elbo_third_term + alpha_lj/2
#       #This conditional technically should not be here according to the formula
#       #but avoids the ELBO being undefined
#       if (alpha_lj != 0){
#         elbo_third_term = elbo_third_term + alpha_lj * log(pi[j]/alpha_lj)
#       }
#       else
#         print("Error: third term of ELBO undefined due to alpha_lj=0")
#     }
#     #print("elbo third term after one loop")
#     #print(elbo_third_term)
#   }
#   #print('elbo third term:')
#   #print(elbo_third_term)
#   #print('susie V:')
#   #print(susie_fit$V)
#   return (elbo_first_term + elbo_second_term + elbo_third_term)
# }

#Computes vectorized ELBO for basic SuSiE-Ann (one locus)
elbo_basic <- function(X,Y,A,annotation_weights,susie_fit){
  L <- nrow(susie_fit$var_lbf)
  p <- ncol(susie_fit$var_lbf)
  #Calculating pi and alpha from annotation weights and annotations
  #in a numerically stable way (logsumexp)
  pi <- softmax(A %*% annotation_weights)
  log_pi <- as.vector(log(pi))
  log_pi_bf <- log_pi + susie_fit$var_lbf
  log_alpha <- log_pi_bf - rowLogSumExps(log_pi_bf)
  alpha <- exp(log_alpha)
  
  #Calculating first term of ELBO
  post_avg_coeff <- alpha * susie_fit$mu
  post_avg_b <- colSums(post_avg_coeff)
  post_avg_b = matrix(post_avg_b, nrow=length(post_avg_b), ncol=1)
  residual_variance <- susie_fit$sigma2
  elbo_first_term <- -1/(2*residual_variance) * norm(Y - X %*% post_avg_b, type="2")^2

  #Calculating second term of ELBO
  post_average_pred_by_effect <- X %*% t(post_avg_coeff)
  total_pred_by_effect_norm <- sum(colNorms(post_average_pred_by_effect)^2)
  scaled_elbo_second_term <- (-1) * total_pred_by_effect_norm
  X_colnorms_2 = colNorms(X)^2
  scaled_elbo_second_term = scaled_elbo_second_term + sum(X_colnorms_2 * (alpha * susie_fit$mu2))
  elbo_second_term <- -1/(2*residual_variance) * scaled_elbo_second_term
  
  #Calculating third term of ELBO
  moment_ratios <- susie_fit$mu2/susie_fit$V
  log_variance_ratios <- log((susie_fit$mu2 - susie_fit$mu^2)/susie_fit$V)
  elbo_third_term <- sum((alpha/2) * (log_variance_ratios - moment_ratios + 1))
  for (l in 1:L){
    elbo_third_term = elbo_third_term - quiet(KL(rbind(alpha[l,], as.vector(pi)), unit="log"))
  }
  #print("Overall ELBO:")
  #print(elbo_first_term + elbo_second_term + elbo_third_term)
  
  #Ensuring each term of ELBO is computed properly and is defined
  if (is.na(elbo_first_term))
    print("First term of ELBO is NA")
  if (is.na(elbo_second_term))
    print("Second term of ELBO is NA")
  if (is.na(elbo_third_term))
    print("Third term of ELBO is NA")
  return (elbo_first_term + elbo_second_term + elbo_third_term)
}

#Calculating ELBO for extended SuSiE-Ann (not vectorized); 
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
  #Calculating third term of ELBO (note: the code here is being updated)
  elbo_third_term <- 0
  #The beta_l and alpha_lj conditionals technically should not be here according to the formula
  #but prevent the ELBO from being undefined; this will be fixed by calling a KL divergence function
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

