#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

//Sample multivariate normal variates ARMA
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmvnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}

//Inverse matrix with arma
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
arma::mat armaSolve(arma::mat A) {
  arma::mat AA(A);
  arma::mat Ainv = arma::inv_sympd(AA);
  return Ainv;
}