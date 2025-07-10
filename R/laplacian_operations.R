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

new_laplacian = function(L, lambda) {
    l <- list(L=lambda*L, L_orig=L, lambda=lambda)
    class(l) <- "laplacian"
    if(is.sparseMatrix(L)) {
        l$solver = new("EigenSparseSolver")
        class(l) = c("laplacian_sparse", class(l))
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
compute.laplacian_dense <- function(l, A) {
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
solve.laplacian_dense <- function(l, b) {
  EigenDenseSolver_solve(l$solver@pointer, b)
}
