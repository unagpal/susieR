#' Consistent implementation of diag
fix_diag=function(x) {
  if(length(x)==1) matrix(x) else diag(x)
}

unscale = function(x) {
  x = sweep(x, 1, attr(x, "scaled:scale"), "*")
  sweep(x, 1, attr(x, "scaled:center"), "+")
}

#' Simple SVD based imputation of missing genotypes
#'
#' @param geno [samples] x [SNPs] genotype matrix (0/1/2)
#' @param prop_var Proportion of variance that the PCs should explain
#' @return Complete genotype matrix.
easy_impute=function(geno, prop_var=0.95) {
  temp=geno
  temp=t(scale(t(geno)))
  temp[is.na(temp)]=0
  s=svd(temp)
  v=s$d^2/sum(s$d^2)
  to_use=cumsum(v)<prop_var
  s$d[!to_use]=0.0
  recon=s$u %*% fix_diag(s$d) %*% t(s$v)
  temp[is.na(geno)]=recon[is.na(geno)]
  temp=unscale(temp)
  stopifnot(max(abs(temp[!is.na(geno)]-geno[!is.na(geno)]))<1e-10)
  temp=round(temp)
  class(temp)="integer"
  temp
}
