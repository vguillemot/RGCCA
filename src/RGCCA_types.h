// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// find how to use them to simplify code
typedef Eigen::SparseMatrix<double> EigenSparseMatrix;
typedef Eigen::SimplicialLDLT<EigenSparseMatrix> EigenSparseSolver;
//typedef Eigen::ConjugateGradient<EigenSparseMatrix, Eigen::Lower|Eigen::Upper>
//    EigenSparseSolver; // causes a segfault
typedef Eigen::LDLT<Eigen::MatrixXd> EigenDenseSolver;
