% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie_ss.R
\name{check_semi_pd}
\alias{check_semi_pd}
\title{Check whether A is positive semidefinite}
\usage{
check_semi_pd(A, tol)
}
\arguments{
\item{A}{a symmetric matrix}
}
\value{
a list of result: \cr
\item{matrix}{The matrix with eigen decomposition}
\item{status}{whether A is positive semidefinite}
\item{eigenvalues}{eigenvalues of A truncated by r_tol}
}
\description{
Check whether A is positive semidefinite
}
