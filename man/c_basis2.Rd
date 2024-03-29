% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{c_basis2}
\alias{c_basis2}
\title{CB construction}
\usage{
c_basis2(A, eps_limit = 1e-10)
}
\arguments{
\item{A}{[nxk matrix] must have rank k}

\item{eps_limit}{[double] used as a limit of small value}
}
\value{
T (n x n) the basis matrix
}
\description{
The SVD-based constraint basis construction for non-overlapping
subsets of constraints. Algorithm 1 from the reference.
Creating a basis from the matrix A
}
