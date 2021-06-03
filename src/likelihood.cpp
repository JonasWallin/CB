
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]


//[[Rcpp::export]]
Eigen::SparseMatrix<double>   likelihood(Eigen::SparseMatrix<double>  Q) {
  Eigen::AMDOrdering<int> ordering;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm;
  
  //ordering(Q, perm); // Call AMD
  ordering(Q.selfadjointView<Eigen::Lower>(), perm);
  return(Q);
}