// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::colvec armaRowSums(const arma::mat& x) {
  return arma::sum(x, 1);
}


// [[Rcpp::export]]
arma::mat eigenMatMult(const arma::mat& A, const arma::mat& B){
    return A * B;
}

