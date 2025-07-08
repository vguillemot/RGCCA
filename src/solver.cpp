// [[Rcpp::depends(RcppEigen)]]

#include "RGCCA_types.h"

//============================================================================
//Sparse solver
//============================================================================

// [[Rcpp::export]]
Rcpp::XPtr<EigenSparseSolver> EigenSparseSolver_new() {
    Rcpp::XPtr<EigenSparseSolver> ptr (new EigenSparseSolver, true);
    return ptr;
}

// [[Rcpp::export]]
void EigenSparseSolver_compute(Rcpp::XPtr<EigenSparseSolver> ptr,
        const EigenSparseMatrix& L) {
    ptr->compute(L);
}

// [[Rcpp::export]]
Eigen::VectorXd EigenSparseSolver_solve(Rcpp::XPtr<EigenSparseSolver> ptr, 
        const Eigen::VectorXd& b) {
    return ptr->solve<Eigen::VectorXd>(b);
}

//============================================================================
//Dense solver
//============================================================================

// [[Rcpp::export]]
Rcpp::XPtr<EigenDenseSolver> EigenDenseSolver_new() {
    Rcpp::XPtr<EigenDenseSolver> ptr (new EigenDenseSolver, true);
    return ptr;
}

// [[Rcpp::export]]
void EigenDenseSolver_compute(Rcpp::XPtr<EigenDenseSolver> ptr,
        const Eigen::MatrixXd& L) {
    ptr->compute(L);
}

// [[Rcpp::export]]
Eigen::VectorXd EigenDenseSolver_solve(Rcpp::XPtr<EigenDenseSolver> ptr,
        const Eigen::VectorXd& b) {
    return ptr->solve<Eigen::VectorXd>(b);
}
