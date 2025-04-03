#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<int> laplacian;
typedef Eigen::SimplicialLDLT<laplacian> LDLTsolver;

//RCPP_MODULE(solver_module) {

    //Rcpp::class_<LDLTsolver>( "SimplicialLDLT" )

    //.constructor()
    
    //.method("compute", &LDLTsolver::compute)
    //.method("solve", &LDLTsolver::solve<laplacian, double>)
    //;
//}

// [[Rcpp::export()]]
Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> LDLTsolver_new() {
    Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
        ptr (new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(), true);
    return ptr;
}

// [[Rcpp::export()]]
void LDLTsolver_compute(Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> ptr,
        const Eigen::SparseMatrix<double>& L) {
    ptr->compute(L);
}

// [[Rcpp::export()]]
Eigen::VectorXd LDLTsolver_solve(Rcpp::XPtr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> ptr,
        const Eigen::VectorXd& b) {
    return ptr->solve<Eigen::VectorXd>(b);
}
