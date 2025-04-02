#define ARMA_USE_SUPERLU 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List rcpp_sparse_eigen(arma::sp_mat X, int rank) {
    arma::mat eigvec;
    arma::vec eigval;
    arma::eigs_sym(eigval, eigvec, X, rank);
    return Rcpp::List::create(Rcpp::Named("s")=eigval, Rcpp::Named("v")=eigvec);
}

// [[Rcpp::export]]
//arma::vec rcpp_solve(arma::sp_mat A, arma::vec b) {
    //return arma::spsolve(A, b);
//}
