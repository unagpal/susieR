#' Optimization of SuSiE-Ann model via alternating calls to SuSiE-Ann IBSS
#' and gradient ascent optimization of ELBO with respect to annotation 
#' weights w. G will denote the number of loci.
#' SuSiE-Ann basic and Proposal 0 will be used interchangably
#' Extended SuSiE-Ann and Proposal 2 will be used interchangably

library(rootSolve)
#library(faux)
library(LaplacesDemon)
library(clusterGeneration)
#library(SimDesign)

#source("SuSiE-Ann/susieR/R/set_X_attributes.R")
#source("SuSiE-Ann/susieR/R/initialize.R")
#source("SuSiE-Ann/susieR/R/update_each_effect.R")
#source("SuSiE-Ann/susieR/R/estimate_residual_variance.R")
#source("SuSiE-Ann/susieR/R/susie_utils.R")
source("~/SuSiE-Ann/susieR/R/susie.R")
source("~/SuSiE-Ann/susieR/R/susie_ann_elbo.R")
source("~/SuSiE-Ann/susieR/R/susie_ann_utils.R")

#' Sigmoid function
sigmoid <- function(vec){
  return (1/(1+exp(-vec)))
}

#' Obtaining pi from annotation weights w and annotation weights A
#' in Proposal 0 (SuSiE-Ann Basic)
pi_from_annotation_weights <- function(A, annotation_weights){
  return (exp(A %*% annotation_weights)/sum(exp(A %*% annotation_weights)))
}

#' Obtaining pi and rho from annotation weights w and annotations A
#' (and number of effects L) for extended SuSiE-Ann (Proposal 2)
pi_rho_from_annotation_weights <- function(A, annotation_weights, L){
  tilda_pi = 1 - sigmoid((-A %*% annotation_weights[1:ncol(A)]))^(1/L)
  rho = sum(tilda_pi)
  pi = tilda_pi/rho
  return (list(pi=pi, rho=rho))
}

#' Function whose root/zero is taken in initializing annotation weights for 
#' SuSiE-Ann Proposal 2 such that the resulting pi is uniform
difference_from_specified_rho <- function(A,equal_annotation_weight, L, specified_rho){
  all_differences <- rep(0, length(equal_annotation_weight))
  for (i in 1:length(equal_annotation_weight)){
    annotation_weights <- c(rep(equal_annotation_weight[i], ncol(A)-1), c(0))
    all_differences[i] = pi_rho_from_annotation_weights(A, annotation_weights, L)$rho - specified_rho
  }
  return (all_differences)
}

#' Given A and rho (in SuSiE-Ann proposal 2), solves for the initial annotation weights
#' where each annotation is weighted equally resulting in the specified value of rho
init_annotation_weights_from_rho <- function(A, L, specified_rho){
  solved_annotation_weight <- uniroot.all(difference_from_specified_rho, interval=c(-100,100), A=A, L=L, specified_rho=specified_rho)[1]
  return (c(rep(solved_annotation_weight, ncol(A)-1), c(0)))
}

#' Generates a random variable with a specified correlation to an existing variable
#' Taken from https://stats.stackexchange.com/questions/15011/
#' generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/313138
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

#Function for generating an easy synthetic dataset with p=2 and G=1
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

#n: number of observations
#G: number of loci
#Data-generating function to examine changes in SuSiE-Ann inference
#quality as a function of number of loci G and annot_sparsity, the Bernoulli probability
#determining which of the two annotations is informative about variable
#inclusion.
motivating_example <- function(n, G, annot_sparsity){
  p <- 11
  all_active_snp_pairs <- list()
  snp_pair_ctr <- 1
  for (i in 1:(p-1)){
    for (j in 1:(p-1)){
      if (i != j){
        all_active_snp_pairs[[snp_pair_ctr]] <- c(i,j)
        snp_pair_ctr = snp_pair_ctr + 1
      }
    }
  }
  print(length(all_active_snp_pairs))
  locus_active_snp_pair_ind <- sample(x=1:length(all_active_snp_pairs), size=G, replace = FALSE)
  locus_active_snp_pairs <- all_active_snp_pairs[locus_active_snp_pair_ind]
  print(locus_active_snp_pairs)
  A_lst <- list()
  for (g in 1:G){
    active_snp_pairs <- locus_active_snp_pairs[g][[1]]
    A_g <- matrix(-20, nrow=p, ncol=2)
    informative_annot <- sample(1:2, size=1, prob=c(annot_sparsity, 1-annot_sparsity))
    if (informative_annot == 1){
      uninformative_annot <- 2
    }
    else{
      uninformative_annot <- 1
    }
    A_g[active_snp_pairs[1], informative_annot] = 2
    A_g[active_snp_pairs[1], uninformative_annot] = 0
    A_g[active_snp_pairs[2], informative_annot] = 1
    A_g[active_snp_pairs[2], uninformative_annot] = 0
    A_lst[[g]] <- A_g
  }
  w <- c(1,1)
  pi <- list()
  for (l in 1:G){
    pi[[l]] <- softmax(A_lst[[l]] %*% w)
  }
  L <- 1
  num_loci <- 2
  susie_ann_X <- list()
  susie_ann_Y <- list()
  output <- list()
  all_selected_variables <- list()
  for (k in 1:G){
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
    X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
    #X <- rnorm_multi(n = n,
    #                 mu = 0,
    #                 sd = 1,
    #                 r = random_corr_matrix,
    #                 varnames = c(1:p),
    #                 empirical = FALSE)
    b <- as.matrix(rep(0, p))
    for (l in 1:L){
      selected_variable <- rmultinom(1,1,pi[[k]])
      selected_var_ind <- match(1, selected_variable)
      selected_locus_variables <- matrix(selected_variable, nrow=1, ncol=length(selected_variable))
      #selected_locus_variables[l,] <- selected_variable
      b[selected_var_ind] = rnorm(1)
    }
    all_selected_variables[[k]] <- selected_locus_variables
    #var_xb <- var(as.matrix(X) %*% b)
    #sigma <- sqrt((1/pve - 1)*var_xb)
    #print(as.matrix(X))
    #print(b)
    Y = as.matrix(X) %*% b + rnorm(n, sd=1)
    susie_ann_X[[k]] <- as.matrix(X)
    susie_ann_Y[[k]] <- as.matrix(Y)
    #print(Y)
    output[[k]] <- list(X=as.matrix(X), Y=as.matrix(Y), b=b, sigma=sigma)
  }
  return (list(susie_ann_X=susie_ann_X, susie_ann_Y=susie_ann_Y, A_lst=A_lst, w=w, pi=pi, output=output, selected_variables=all_selected_variables))
}

#Smaller-scale version of motivating_example function above for up to 6 loci
#Used in generating 'Motivating Example' in Figure 2 of Graphical Models report
#and SuSiE-Ann lab meeting presentation
contrived_example <- function(n, num_loci){
  #num_loci <- 2
  all_active_snp_pairs <- list()
  snp_pair_ctr <- 1
  for (i in 1:3){
    for (j in 1:3){
      if (i != j){
        all_active_snp_pairs[[snp_pair_ctr]] <- c(i,j)
        snp_pair_ctr = snp_pair_ctr + 1
      }
    }
  }
  all_locus_snp_pairs <- all_active_snp_pairs[1:num_loci]
  A_lst <- list()
  for (g in 1:num_loci){
    active_snp_pairs <- all_locus_snp_pairs[g][[1]]
    A_g <- matrix(-20, nrow=4, ncol=2)
    if (g %% 2){
      informative_annot <- 1
      uninformative_annot <- 2
    }
    else{
      informative_annot <- 2
      uninformative_annot <- 1
    }
    A_g[active_snp_pairs[1], informative_annot] = 2
    A_g[active_snp_pairs[1], uninformative_annot] = 0
    A_g[active_snp_pairs[2], informative_annot] = 1
    A_g[active_snp_pairs[2], uninformative_annot] = 0
    A_lst[[g]] <- A_g
  }
  #A_1 <- matrix(c(2,1,-20,-20, 0,0,-20,-20), nrow=4, ncol=2)
  #A_2 <- matrix(c(-20,0,0,-20,-20,1,2,-20), nrow=4, ncol=2)
  #A_lst <- list()
  #A_lst[[1]] <- A_1
  #A_lst[[2]] <- A_2
  w <- c(1,1)
  pi <- list()
  for (l in 1:num_loci){
    pi[[l]] <- softmax(A_lst[[l]] %*% w)
  }
  L <- 1
  p <- 4
  susie_ann_X <- list()
  susie_ann_Y <- list()
  output <- list()
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
    X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
    #X <- rnorm_multi(n = n,
    #                 mu = 0,
    #                 sd = 1,
    #                 r = random_corr_matrix,
    #                 varnames = c(1:p),
    #                 empirical = FALSE)
    b <- as.matrix(rep(0, p))
    for (l in 1:L){
      selected_variable <- rmultinom(1,1,pi[[k]])
      selected_var_ind <- match(1, selected_variable)
      selected_locus_variables <- matrix(selected_variable, nrow=1, ncol=length(selected_variable))
      #selected_locus_variables[l,] <- selected_variable
      b[selected_var_ind] = rnorm(1)
    }
    all_selected_variables[[k]] <- selected_locus_variables
    #var_xb <- var(as.matrix(X) %*% b)
    #sigma <- sqrt((1/pve - 1)*var_xb)
    #print(as.matrix(X))
    #print(b)
    Y = as.matrix(X) %*% b + rnorm(n, sd=1)
    susie_ann_X[[k]] <- as.matrix(X)
    susie_ann_Y[[k]] <- as.matrix(Y)
    #print(Y)
    output[[k]] <- list(X=as.matrix(X), Y=as.matrix(Y), b=b, sigma=sigma)
  }
  return (list(susie_ann_X=susie_ann_X, susie_ann_Y=susie_ann_Y, A_lst=A_lst, w=w, pi=pi, output=output, selected_variables=all_selected_variables))
}

#Function generating synthetic data satisfying a given PVE, n, p, L, G, and 
#number of annotations K. Currently joint correlation matrix between SNPs is determined
#correlation matrix with random (positive) Eigenvalues. One alternative could be:
#Correlation absolute values ~ Beta(2,5)
#Correlation signs ~ Bernoulli(0.5)
generate_synthetic_data <- function(pve, n, p, L, num_loci, num_annotations){
  output = list()
  w <- rnorm(num_annotations, mean = 0, sd = 1.5)
  
  A_lst <- list()
  pi <- list()
  z_score_cutoff <- (1-0.1)/1.28
  annotations_mean <- rep(0.1, num_annotations)
  for (l in 1:num_loci){
    random_cov_matrix <- genPositiveDefMat(dim=num_annotations, covMethod="eigen")$Sigma
    random_cov_matrix = random_cov_matrix*sqrt(z_score_cutoff)/max(random_cov_matrix)
    A_l = matrix(0, nrow=p, ncol=num_annotations)
    for (i in 1:p)
      annot_row <- rmvnorm(n=1, mean=annotations_mean, sigma=random_cov_matrix)
      for (k in 1:num_annotations){
        if (annot_row[k] >= 1){
          annot_row[k] = 1
        }
        else{
          annot_row[k] = 0
        }
      }
      A_l[i,] <- annot_row
    A_lst[[l]] <- A_l
    pi[[l]] <- softmax(A_lst[[l]] %*% w)
  }
  print("Done creating A")
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
        random_corr_matrix[i,j] = random_cov_matrix[i,j]/(sqrt(random_cov_matrix[i,i])*sqrt(random_cov_matrix[j,j]))
      }
    }
    random_corr_matrix = round(random_corr_matrix,3)
    X <- round(rmvnorm(n=n, mean=rep(1, p), sigma=random_corr_matrix ))
    print(dim(X))
    for (obs in 1:n){
      for (var in 1:p){
        if (X[obs, var]<0){
          X[obs, var] = 0
        }
        else if (X[obs, var]>2){
          X[obs, var]= 2
        }
      }
    }
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
    print("Done creating X, Y, b for one locus")
  }
  return (list(susie_ann_X=susie_ann_X, susie_ann_Y=susie_ann_Y, A_lst=A_lst, w=w, pi=pi, output=output, selected_variables=all_selected_variables))
}

#Function generating synthetic data satisfying a given PVE, n, p, L, G, and 
#number of annotations K. Currently joint correlation matrix between SNPs is determined
#correlation matrix with random (positive) Eigenvalues. One alternative could be:
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

#' Key function executing alternating optimization of SuSiE-Ann model. Optimization of ELBO
#' with respect to w is in susie_ann_elbo.R. Parameters:
#' extended_model: FALSE indicates basic SuSiE-Ann/Proposal 0; TRUE indicates extended SuSiE-Ann/Proposal 2
#' elbo_optim_method: "SGD", "L-BFGS-B", or "Adam" can be used to optimize ELBO w.r.t. w
#' opt_alpha: whether optimization of annotation weights w incorporates its effect on alpha
susie_ann <- function(X_lst,Y_lst, A_lst, num_loci, extended_model, batch_size,
                  annotation_weights=NULL, rho=0.2, 
                  L = min(10,ncol(X)),
                  scaled_prior_variance=1.0, residual_variance=NULL,
                  prior_weights=NULL, null_weight=NULL,
                  standardize=TRUE,intercept=TRUE,
                  estimate_residual_variance=TRUE,
                  estimate_prior_variance = FALSE,
                  estimate_prior_method = c("optim","EM","simple"),
                  check_null_threshold=0, prior_tol=1E-9,
                  residual_variance_upperbound = Inf,
                  s_init_lst = NULL,coverage=0.95,min_abs_corr=0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE, 
                  susie_ann_opt_itr=100,
                  elbo_optim_method="L-BFGS-B", optim_alpha=TRUE,
                  elbo_opt_steps_per_itr=100,
                  reg_type="None",
                  max_iter=100,tol=1e-3,
                  verbose=FALSE,track_fit=FALSE) {
  if (reg_type == "L1"){
    l1_sigma_2 <- 4
    l2_sigma_2 <- 0
  }
  else if (reg_type=="L2"){
    l1_sigma_2 <- 0
    l2_sigma_2 <- 4
  }
  else{
    l1_sigma_2 <- 0
    l2_sigma_2 <- 0
  }
  
  num_loci <- length(X_lst)
  #Adds an intercept term for w through an extra annotation of value 1
  #for SuSiE-Ann proposal 2
  if (intercept){
    for (l in 1:num_loci)
      Y_lst[[l]] <- Y_lst[[l]] - mean(Y_lst[[l]])
  }
  if (extended_model){
    for (l in 1:num_loci)
      A_lst[[l]] = cbind(A_lst[[l]], rep(1, nrow(A_lst[[l]])))
  }
  else
    all_iter_rho <- NULL
  all_initial_annotation_weight_elbo <- rep(0, susie_ann_opt_itr)
  all_opt_annotation_weight_elbo <-  rep(0, susie_ann_opt_itr)
  #Initialization of annotation weights for SuSiE-Ann proposal 2 such that
  #the resulting pi is uniform
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
  #Initialization of annotation weights to zero vector for SuSiE-Ann proposal 0
  #if they are null
  else{
    if (is.null(annotation_weights)){
      print("ncol Alst")
      print(ncol(A_lst[[1]]))
      annotation_weights = rep(0, ncol(A_lst[[1]]))
    }
  }
  all_itr_pi <- list()
  all_itr_alpha <- list()
  #Loop over SuSiE-Ann alternating optimization iterations
  for (susie_ann_itr in 1:susie_ann_opt_itr){
    #Training SuSiE models from initialization in case partially trained
    #SuSiE models are passed in, or in alternating SuSiE-Ann optimization
    #starting the second iteration when the previous iteration's SuSiE
    #models are available
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
                 s_init = s_init_lst[[l]], coverage=coverage,min_abs_corr=min_abs_corr,
                 compute_univariate_zscore = compute_univariate_zscore,
                 na.rm = na.rm, max_iter=max_iter,tol=tol,
                 verbose=verbose,track_fit=track_fit)
        print("SuSiE ELBO during IBSS:")
        print(s_init_lst[[l]]$elbo)
      }
    }
    #Training SuSiE models from scratch in case no pretrained SuSiE models
    #are passed in and it is the first alternating optimization iteration
    else{
      s_init_lst <- list()
      
      for (l in 1:num_loci){
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
        print("SuSiE ELBO during IBSS:")
        print(s_init_lst[[l]]$elbo)
      }
      print("Done initializing SuSiE in first iteration")
    }
    #Saving raw SuSiE results as a baseline for SuSiE-Ann
    if (susie_ann_itr == 1)
      raw_susie_results <- s_init_lst
    #Gradient-based optimization of annotation weights
    #If desired, can run until convergence only in the last iteration
    print("Optimizing annotation weights")
    if (susie_ann_itr != susie_ann_opt_itr)
      opt_annot_weights_results <- gradient_opt_annotation_weights(X_lst,Y_lst,A_lst,extended_model,elbo_optim_method, optim_alpha, annotation_weights,s_init_lst,batch_size, l1_sigma_2, l2_sigma_2, elbo_opt_steps_per_itr, run_until_convergence=TRUE)
    else
      opt_annot_weights_results <- gradient_opt_annotation_weights(X_lst,Y_lst,A_lst,extended_model, elbo_optim_method, optim_alpha, annotation_weights,s_init_lst,batch_size, l1_sigma_2, l2_sigma_2, elbo_opt_steps_per_itr, run_until_convergence=TRUE)
    print("finished optimizing annotation weights")
    #Retrieving initial ELBO, final ELBO, optimized w, and optimized pi
    #from optimization of ELBO w.r.t. w
    all_initial_annotation_weight_elbo[susie_ann_itr] = opt_annot_weights_results$initial_elbo
    all_opt_annotation_weight_elbo[susie_ann_itr] = opt_annot_weights_results$elbo
    annotation_weights <- opt_annot_weights_results$annot_weights
    print("Returned annotation weights:")
    print(annotation_weights)
    pi <- opt_annot_weights_results$pi
    #Updating rho if extended SuSiE-Ann/proposal 2 is being used
    if (extended_model){
      rho <- updated_pi_rho$rho
      all_iter_rho[susie_ann_itr] = rho
      for (i in 1:num_loci){
        s_init_lst[[i]]$rho=rho
      }
    }
    #Tracking pi for the first locus across alternating optimization iterations
    all_itr_pi[[susie_ann_itr]] <- pi[[1]]
    
    #s_init$pi = opt_annot_weights_results$pi/sum(opt_annot_weights_results$pi)
    #s_init$rho = min(c(opt_annot_weights_results$rho, 1))
    #print("pi before gradient update:")
    #print(s_init$pi)
    
    #Updating pi and alpha of SuSiE models based on optimization of ELBO w.r.t. w
    for (i in 1:num_loci){
      s_init_lst[[i]]$pi=pi[[i]]
      s_init_lst[[i]]$alpha = opt_annot_weights_results$alpha[[i]]
    }
    
    #Tracking alpha for the first locus across alternating optimization iterations
    all_itr_alpha[[l]] <- s_init_lst[[1]]$alpha[[1]]
    #Comparing ELBO prior to optimization of ELBO w.r.t. w to 
    #ELBO after optimization of ELBO w.r.t. w
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
  return(list(susie_models=s_init_lst, raw_susie_results=raw_susie_results, w=annotation_weights, final_pi=pi, all_iter_pi=all_itr_pi, all_iter_alpha=all_itr_alpha, initial_elbo_values=all_initial_annotation_weight_elbo, opt_elbo_values=all_opt_annotation_weight_elbo, all_rho=all_iter_rho))
}

#Analyzes and plots figures illustrating the quality of a SuSiE-Ann fit to synthetic data, including:
#a) plots of SuSiE-Ann ELBO over alternating optimization iterations
#b) Comparison of actual to predicted pi (e.g. locus R coefficients, locus regression p value, KL div)
#c) Comparison of actual to predicted w (plot, correlation coefficient, and correlation coefficient p value)
#d) Comparison of predicted to 'true' alpha (based on averaging one-hot vectors across effects):
#e.g. correlation, average SuSiE-Ann alpha of 'true' effect variables
#e) Comparison of SuSiE and SuSiE-Ann predictions of Y across loci
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
  selected_loci_graphs <- sample(1:num_loci, 4)
  for (l in 1:num_loci){
    all_corr_pi[l] <- cor(pred_pi[[l]], actual_pi[[l]])
    fit <- lm(actual_pi[[l]] ~ pred_pi[[l]])
    #print(summary(fit)$coefficients)
    if (var(pred_pi[[l]])==0)
      all_pval_pi[l] <- 1
    else
      all_pval_pi[l] <- summary(fit)$coefficients[2,4]
    #print("Now indiv")
    #print(summary(fit)$coefficients[2,4])
    #all_pval_pi[l] <- summary(fit)$coefficients[2,4]
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
  all_pred_alpha <- list()
  raw_susie_pred_alpha <- list()
  all_true_alpha_dist <- list()
  all_true_effect_vars <- list()
  average_alpha_corr <- rep(0, num_loci)
  average_alpha_pval <- rep(0, num_loci)
  raw_susie_average_alpha_corr <- rep(0, num_loci)
  raw_susie_average_alpha_pval <- rep(0, num_loci)
  L <- nrow(synth_data$selected_variables[[1]])
  activated_variable_alpha <- rep(0, num_loci * L)
  raw_susie_activated_variable_alpha <- rep(0, num_loci * L)
  raw_susie_results <- res$raw_susie_results
  for (s in 1:length(susie_models)){
    all_pred_alpha[[s]] <- susie_models[[s]]$alpha
    raw_susie_pred_alpha[[s]] <- raw_susie_results[[s]]$alpha
    all_true_alpha_dist[[s]] <- synth_data$selected_variables[[s]]
    locus_true_effect_vars <- rep(0, L)
    remaining_effects <- c(1:L)
    raw_susie_remaining_effects <- c(1:L)
    for (r in 1:L){
      locus_true_effect_vars[r] <- match(1, all_true_alpha_dist[[s]][r,])
      remaining_effect_alpha <- rep(0, length(remaining_effects))
      remaining_effect_ctr <- 1
      for (l in remaining_effects){
        remaining_effect_alpha[remaining_effect_ctr] <- all_pred_alpha[[s]][l,locus_true_effect_vars[r]]
        remaining_effect_ctr = remaining_effect_ctr + 1
        max_remaining_effect_alpha <- max(remaining_effect_alpha)
        #argmax_remaining_effect_alpha <- match(max_remaining_effect_alpha, remaining_effect_alpha)
        remaining_effects <- remaining_effects[!remaining_effects == max_remaining_effect_alpha]
        activated_variable_alpha[r*(l-1) + l] <- max_remaining_effect_alpha
      }
      raw_susie_remaining_effect_alpha <- rep(0, length(raw_susie_remaining_effects))
      raw_susie_remaining_effect_ctr <- 1
      for (l in raw_susie_remaining_effects){
        raw_susie_remaining_effect_alpha[raw_susie_remaining_effect_ctr] <- raw_susie_pred_alpha[[s]][l,locus_true_effect_vars[r]]
        raw_susie_remaining_effect_ctr = raw_susie_remaining_effect_ctr+1
        raw_susie_max_remaining_effect_alpha <- max(raw_susie_remaining_effect_alpha)
        raw_susie_remaining_effects <- raw_susie_remaining_effects[!raw_susie_remaining_effects==raw_susie_max_remaining_effect_alpha]
        raw_susie_activated_variable_alpha[r*(l-1) + l] <- raw_susie_max_remaining_effect_alpha
      }
    }
    all_true_effect_vars[[s]] <- locus_true_effect_vars
    average_true_alpha <- colMeans(all_true_alpha_dist[[s]])
    average_pred_alpha <- colMeans(all_pred_alpha[[s]])
    average_raw_susie_pred_alpha <- colMeans(raw_susie_pred_alpha[[s]])
    average_alpha_corr[s] <- cor(average_pred_alpha, average_true_alpha)
    fit <- lm(average_true_alpha ~ average_pred_alpha)
    if (var(average_pred_alpha) > 0)
      average_alpha_pval[s] <- summary(fit)$coefficients[2,4]
    else
      average_alpha_pval[s] <- 1
    # print("summary fit coefficients:")
    # print(summary(fit)$coefficients)
    # print("Var of average pred alpha")
    # print(var(average_pred_alpha))
    #average_alpha_pval[s] <- summary(fit)$coefficients[2,4]
    raw_susie_average_alpha_corr[s] <- cor(average_raw_susie_pred_alpha, average_true_alpha)
    fit <- lm(average_true_alpha ~ average_raw_susie_pred_alpha)
    # print("summary fit coefficients: ")
    # print(summary(fit)$coefficients)
    # print("made it here 1")
    raw_susie_average_alpha_pval[s] <- summary(fit)$coefficients[2,4]
    # print("made it here 2")
    # print(raw_susie_average_alpha_pval[s])
    if (s %in% selected_loci_graphs){
      plot(average_pred_alpha, average_true_alpha, xlab="Average SuSiE-Ann predicted alpha for SNP across effects", ylab="True alpha for SNP")
      plot(average_raw_susie_pred_alpha, average_true_alpha, xlab="Average SuSiE predicted alpha for SNP across effects", ylab="True alpha for SNP")
    }
  }
  print("Correlation coefficients between predicted average alpha and true average alpha for each locus:")
  print(average_alpha_corr)
  print("These linear regressions have p values of:")
  print(average_alpha_pval)
  print("Average SuSiE Ann Alpha for Activated Variables (Generous Interpretation*):")
  print(mean(activated_variable_alpha))
  hist(activated_variable_alpha, xlab="SuSiE-Ann Alpha for Activated Variables (Generous Interpretation*)")
  #Next comparing results to SuSiE: correlations with true output, inference of alpha
  raw_susie_results <- res$raw_susie_results
  hist(raw_susie_activated_variable_alpha, xlab="SuSiE Alpha for Activated Variables (Generous Interpretation*)")
  print("Raw SuSiE Correlation coefficients between predicted and true average alpha for each locus:")
  print(raw_susie_average_alpha_corr)
  print("These linear regressions have p values of:")
  print(raw_susie_average_alpha_pval)
  print("Average raw SuSiE Alpha of Activated Variables (Generous Interpretation*):")
  print(mean(raw_susie_activated_variable_alpha))
  
  all_prediction_pval <- rep(0, num_loci)
  all_prediction_r <- rep(0, num_loci)
  raw_susie_prediction_pval <- rep(0, num_loci)
  raw_susie_prediction_r <- rep(0, num_loci)
  susie_ann_Y <- synth_data$susie_ann_Y
  for (l in 1:num_loci){
    true_y <- susie_ann_Y[[l]]
    susie_ann_pred_y <- predict(susie_models[[l]])
    susie_pred_y <- predict(raw_susie_results[[l]])
    all_prediction_r[l] <- cor(susie_ann_pred_y, true_y)
    fit = lm(true_y ~ susie_ann_pred_y)
    all_prediction_pval[s] <- summary(fit)$coefficients[2,4]
    raw_susie_prediction_r[l] <- cor(susie_pred_y, true_y)
    fit = lm(true_y ~ susie_pred_y)
    raw_susie_prediction_pval[s] <- summary(fit)$coefficients[2,4]
    if (l %in% selected_loci_graphs){
      plot(susie_ann_pred_y, true_y, xlab="SuSiE-Ann Y hat for one locus", ylab="True Y")
      plot(susie_pred_y, true_y, xlab="SuSiE Y hat for one locus", ylab="True Y")
    }
  }
  print("SuSiE-Ann Y vs. Y hat Correlation Coefficient (R) by Locus:")
  print(all_prediction_r)
  print("SuSiE-Ann Y vs. Y hat Correlation p-value by Locus:")
  print(all_prediction_pval)
  print("SuSiE Y vs. Y hat Correlation Coefficient (R) by Locus:")
  print(raw_susie_prediction_r)
  print("SuSiE Y vs. Y hat Correlation p-value by Locus")
  print(raw_susie_prediction_pval)
}
