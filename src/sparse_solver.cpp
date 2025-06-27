#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// find how to use them to simplify code
typedef Eigen::SparseMatrix<int> laplacian;
typedef Eigen::SimplicialLDLT<laplacian> LDLTsolver;

//============================================================================
//Sparse solver
//============================================================================

// [[Rcpp::export]]
Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> LDLTsparse_new() {
    Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
        ptr (new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(), true);
    return ptr;
}

// [[Rcpp::export]]
void LDLTsparse_compute(Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> ptr,
        const Eigen::SparseMatrix<double>& L) {
    ptr->compute(L);
}

// [[Rcpp::export]]
Eigen::VectorXd LDLTsparse_solve(Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> ptr,
        const Eigen::VectorXd& b) {
    return ptr->solve<Eigen::VectorXd>(b);
}

//============================================================================
//Dense solver
//============================================================================

// [[Rcpp::export]]
Rcpp::XPtr<Eigen::LDLT<Eigen::MatrixXd>> LDLTdense_new() {
    Rcpp::XPtr<Eigen::LDLT<Eigen::MatrixXd>>
        ptr (new Eigen::LDLT<Eigen::MatrixXd>(), true);
    return ptr;
}

// [[Rcpp::export]]
void LDLTdense_compute(Rcpp::XPtr<Eigen::LDLT<Eigen::MatrixXd>> ptr,
        const Eigen::MatrixXd& L) {
    ptr->compute(L);
}

// [[Rcpp::export]]
Eigen::VectorXd LDLTdense_solve(Rcpp::XPtr<Eigen::LDLT<Eigen::MatrixXd>> ptr,
        const Eigen::VectorXd& b) {
    return ptr->solve<Eigen::VectorXd>(b);
}
