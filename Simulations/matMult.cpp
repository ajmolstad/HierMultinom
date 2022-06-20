// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

double Exp(double x) // the functor we want to apply
{
    return std::exp(x);
}

// [[Rcpp::export]]
SEXP eigenMatMult_e(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;
    C.unaryExpr(&Exp);
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
arma::colvec armaRowSums(const arma::mat& x) {
  return arma::sum(x, 1);
}
