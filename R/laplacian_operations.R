# Defining a class to wrap around the pointer to LDLTsolver
# boilerplate taken from the book "seamless R and C++ integration..."
setClass("EigenSparseSolver", slots=c(pointer="externalptr"))
setClass("EigenDenseSolver", slots=c(pointer="externalptr"))

setMethod("initialize", "EigenSparseSolver", function(.Object, ...) {
    #.Object@pointer <- .Call(LDLTsolver_method("new"), ...)
    .Object@pointer <- EigenSparseSolver_new()
    .Object
})

setMethod("initialize", "EigenDenseSolver", function(.Object, ...) {
    #.Object@pointer <- .Call(LDLTsolver_method("new"), ...)
    .Object@pointer <- EigenDenseSolver_new()
    .Object
})

new_laplacian = function(L, lambda, woodbury = TRUE) {
    l <- list(L=lambda*L, L_orig=L, lambda=lambda)
    class(l) <- "laplacian"
    if(is.sparseMatrix(L)) {
        l$solver = new("EigenSparseSolver")
        class(l) = c("laplacian_sparse", class(l))
    } else if (woodbury) {
      class(l)          = c("laplacian_woodbury", class(l))
      eig               = eigen(l$L_orig)
      eig_values_toKeep = which(abs(eig$values) > .Machine$double.eps)
      l$vectors         = eig$vectors[, eig_values_toKeep]
      l$values          = eig$values[eig_values_toKeep]
    } else {
      l$solver = new("EigenDenseSolver")
      class(l) = c("laplacian_dense", class(l))
    }
    return(l)
}

compute <- function(l, A) {
  UseMethod("compute")
}

#' @export
compute.laplacian_sparse <- function(l, A) {
  EigenSparseSolver_compute(l$solver@pointer, A)
    
}

#' @export
compute.laplacian_woodbury <- function(l, A) {
  mu       = l$mu
  tmp_idty = (1/(2*mu))*diag(unique(dim(A))) 
  tmp_C    = diag(l$values/(1+(l$lambda*l$values)/mu))
  Linv     = tmp_idty - (l$lambda/(2*mu^2))*l$vectors %*% tmp_C %*% t(l$vectors)
  return(Linv)
}

#' @export
compute.laplacian_dense <- function(l, A, mu) {
  EigenDenseSolver_compute(l$solver@pointer, A)
    
}

solve <- function(l, b) {
    UseMethod("solve")
}

#' @export
solve.laplacian_sparse <- function(l, b) {
  EigenSparseSolver_solve(l$solver@pointer, b)
}

#' @export
solve.laplacian_woodbury <- function(l, b) {
  return(l$Linv %*% b)
}

#' @export
solve.laplacian_dense <- function(l, b) {
  EigenDenseSolver_solve(l$solver@pointer, b)
}
