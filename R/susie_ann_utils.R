#logsumexp and softmax taken from https://gist.github.com/aufrank/83572
#note: numerical stability is important for softmax
library(mvtnorm)
library(matrixcalc)

logsumexp <- function (x) {
  y = max(x)
  return (y + log(sum(exp(x - y))))
}

#numerically stable softmax function
softmax <- function (x) {
  return (exp(x - logsumexp(x)))
}

generate_synthetic_data_final <- function(pve, n, p, L, G, K, annot_sparsity, method="rounding"){
  output = list()
  A_lst <- list()
  pi <- list()
  susie_ann_X <- list()
  susie_ann_Y <- list()
  all_selected_variables <- list()
  b <- list()
  all_sigma <- list()
  #w is drawn with iid N(0,2) components
  w <- rnorm(K, mean = 0, sd = 2)
  
  #Setting Wishart scale for annotations and SNPs
  #based on https://en.wikipedia.org/wiki/Wishart_distribution "Choice of Parameters"
  all_annotation_covariances <- rWishart(1, K, diag(K)/K)
  all_snp_covariances = rWishart(G, p, diag(p)/p)
  annot_zscore_1 <- qnorm(1-annot_sparsity,mean=0,sd=1)
  if (method == "rounding"){
    marginal_sigma <- sqrt(diag(all_annotation_covariances[, , 1]))
    mu <- 1 - annot_zscore_1 * marginal_sigma
    for (g in 1:G){
      A_g <- rmvnorm(p, mean=mu, sigma=all_annotation_covariances[, , 1])
      A_g[A_g<1]=0
      A_g[A_g>1]=1
      A_lst[[g]]<- A_g
      pi[[g]] <- softmax(A_g %*% w)
    }
    for (g in 1:G){
      #Sampling X_g
      snp_freq = runif(p, min = 0.1, max = 0.5)
      prob_2 = snp_freq^2
      prob_1 = 2 * snp_freq * (1-snp_freq)
      prob_0 = 1 - prob_2 - prob_1
      z_0 = qnorm(prob_0,mean=0,sd=1)
      z_1 = qnorm(prob_0 + prob_1, mean=0,sd=1)
      sigma = 1/(z_1 - z_0)
      mu = (-z_0)/(z_1 - z_0)
      snp_corr_matrix <- cov2cor(all_snp_covariances[, , g])
      #Covariance from correlation matrix: https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
      snp_cov_matrix <- diag(sigma) %*% snp_corr_matrix %*% diag(sigma)
      X_g <- rmvnorm(n, mean=mu, sigma=snp_cov_matrix)
      X_g[X_g<0] = 0
      X_g[abs(X_g-0.5)<0.5] = 1
      X_g[X_g>1] = 2
      susie_ann_X[[g]] <- as.matrix(X_g)
      
      #Sampling b_g
      selected_locus_variables <- rep(0, L)
      b_g <- as.matrix(rep(0, p))
      for (l in 1:L){
        selected_variable <- rmultinom(1,1,pi[[g]])
        selected_var_ind <- match(1, selected_variable)
        selected_locus_variables[l] <- selected_var_ind
        b_g[selected_var_ind] = rnorm(1)
      }
      b[[g]] <- b_g
      all_selected_variables[[g]] <- selected_locus_variables
      var_xb <- var(as.matrix(X_g) %*% b_g)
      sigma_g <- sqrt((1/pve - 1)*var_xb)
      all_sigma[[g]] <- sigma_g
      Y_g = as.matrix(X_g) %*% b_g + rnorm(n, sd=sigma_g)
      susie_ann_Y[[g]] <- as.matrix(Y_g)
      if (g %% 10 == 0){
        print("Finished generating 10 loci")
      }
    }
  }
  else if (method=="exact"){
    #Calculating Gaussian covariance matrix for sampling binary annotations
    large_quantity <- 100
    annot_gaussian_var <- diag(K)
    cutoff_vec <- annot_zscore_1 * rep(1, K)
    binary_search_itr <- 7
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        annot_corr <- (all_annotation_covariances[i, j, 1]/(sqrt(all_annotation_covariances[i, i, 1]*all_annotation_covariances[j, j, 1])))[1]
        lhs <- annot_sparsity^2 + annot_sparsity * (1-annot_sparsity) * annot_corr
        #Approximately solve for gaussian_corr inducing annot_corr via binary search
        annot_corr_mat <- diag(2)
        low <- -1.0
        high <- 1.0
        for (itr in 1:binary_search_itr){
          mid <- (low + high)/2
          annot_corr_mat[1,2] = mid
          annot_corr_mat[2,1] = mid
          rhs <- pmvnorm(mean=rep(0,2), corr=annot_corr_mat, lower=rep(annot_zscore_1, 2), upper=rep(Inf, 2))[1]
          if (lhs > rhs){
            low = mid
          }
          else if (lhs < rhs){
            high = mid
          }
        }
        gaussian_corr = mid
        annot_gaussian_var[i,j] = gaussian_corr
        annot_gaussian_var[j,i] = gaussian_corr
      }
    }
    #Sampling binary annotations
    for (g in 1:G){
      A_g <- rmvnorm(p, mean=rep(0, K), sigma=annot_gaussian_var)
      A_g[A_g < annot_zscore_1] = 0
      A_g[A_g >= annot_zscore_1] = 1
      A_lst[[g]]<- A_g
      pi[[g]] <- softmax(A_g %*% w)
    }
    for (g in 1:G){
      #Sampling X_g
      X_g_gaussian_var <- diag(p)
      snp_freq <- runif(p, min = 0.1, max = 0.5)
      prob_2 <- snp_freq^2
      prob_1 <- 2 * snp_freq * (1-snp_freq)
      prob_0 <- 1 - prob_2 - prob_1
      #Calculating cutoff values for locus g
      c_1 <- qnorm(prob_0,mean=0,sd=1)
      c_2 <- qnorm(prob_0 + prob_1, mean=0,sd=1)
      #Calculating mean & stdev of each SNP value
      snp_mu_g <- prob_1 * 1 + prob_2 * 2
      snp_var_g <- prob_0 * (rep(0,p) - snp_mu_g)^2 + prob_1 * (rep(1,p) - snp_mu_g)^2 + prob_2 * (rep(2,p) - snp_mu_g)^2
      snp_stdev_g <- sqrt(snp_var_g)
      print(snp_stdev_g)
      search_itr <- 200
      for (i in 1:(p-1)){
        for (j in (i+1):p){
          snp_corr <- all_snp_covariances[i,j,g]/(sqrt(all_snp_covariances[i,i,g]*all_snp_covariances[j,j,g]))
          lhs <- snp_mu_g[i] * snp_mu_g[j] + snp_stdev_g[i]*snp_stdev_g[j]*snp_corr
          #Approximately solve for snp_gaussian_corr inducing snp_corr via grid search
          annot_corr_mat <- diag(2)
          itr <- 1
          all_rhs <- rep(0, search_itr)
          for (snp_gaussian_corr in seq(from = -1, to = 1, length.out = search_itr)){
            annot_corr_mat[1,2] = snp_gaussian_corr
            annot_corr_mat[2,1] = snp_gaussian_corr
            annot_corr_mat <- as.matrix(annot_corr_mat)
            if (is.positive.semi.definite(annot_corr_mat)){
              p_11 <- pmvnorm(mean=rep(0,2), corr=annot_corr_mat, lower=c(c_1[i], c_1[j]), upper=c(c_2[i], c_2[j]))
              p_21 <- pmvnorm(mean=rep(0,2), corr=annot_corr_mat, lower=c(c_2[i], c_1[j]), upper=c(Inf, c_2[j]))
              p_12 <- pmvnorm(mean=rep(0,2), corr=annot_corr_mat, lower=c(c_1[i], c_2[j]), upper=c(c_2[i], Inf))
              p_22 <- pmvnorm(mean=rep(0,2), corr=annot_corr_mat, lower=c(c_2[i], c_2[j]), upper=rep(Inf,2))
              all_rhs[itr] <- p_11 + 2 * p_21 + 2 * p_12 + 4 * p_22
            }
            else{
              all_rhs[itr] <- Inf
            }
            itr = itr + 1
          }
          print(all_rhs)
          snp_gaussian_corr <- seq(from = -1, to = 1, length.out = search_itr)[which.min(abs(lhs - all_rhs))]
          X_g_gaussian_var[i,j] = snp_gaussian_corr
          X_g_gaussian_var[j,i] = snp_gaussian_corr
        }
      }
      print(X_g_gaussian_var)
      #Gaussian sampling & discretizing to get X_g
      X_g <- rmvnorm(n, mean=snp_mu_g, sigma=X_g_gaussian_var)
      for (j in 1:p){
        X_g_col_j <- X_g[,j]
        X_g_col_j[X_g_col_j<c_1[j]] = 0
        X_g_col_j[abs(X_g_col_j - (c_1[j] + c_2[j])/2) < (c_2[j]-c_1[j])/2] = 1
        X_g_col_j[X_g_col_j>=c_2[j]] = 2
        X_g[,j] = X_g_col_j
      }
      susie_ann_X[[g]] <- as.matrix(X_g)
      
      #Sampling b_g
      selected_locus_variables <- rep(0, L)
      b_g <- as.matrix(rep(0, p))
      for (l in 1:L){
        selected_variable <- rmultinom(1,1,pi[[g]])
        selected_var_ind <- match(1, selected_variable)
        selected_locus_variables[l] <- selected_var_ind
        b_g[selected_var_ind] = rnorm(1)
      }
      b[[g]] <- b_g
      all_selected_variables[[g]] <- selected_locus_variables
      var_xb <- var(as.matrix(X_g) %*% b_g)
      sigma_g <- sqrt((1/pve - 1)*var_xb)
      all_sigma[[g]] <- sigma_g
      Y_g = as.matrix(X_g) %*% b_g + rnorm(n, sd=sigma_g)
      susie_ann_Y[[g]] <- as.matrix(Y_g)
    }
  }
  return (list(susie_ann_X=susie_ann_X, susie_ann_Y=susie_ann_Y, A_lst=A_lst,w=w, pi=pi, b=b, sigma=all_sigma, selected_variables=all_selected_variables))
}