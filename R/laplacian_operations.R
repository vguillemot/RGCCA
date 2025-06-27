# Defining a class to wrap around the pointer to LDLTsolver
# boilerplate taken from the book "seamless R and C++ integration..."
setClass("LDLTsparse", slots=c(pointer="externalptr"))
setClass("LDLTdense", slots=c(pointer="externalptr"))

setMethod("initialize", "LDLTsparse", function(.Object, ...) {
    #.Object@pointer <- .Call(LDLTsolver_method("new"), ...)
    .Object@pointer <- LDLTsparse_new()
    .Object
})

setMethod("initialize", "LDLTdense", function(.Object, ...) {
    #.Object@pointer <- .Call(LDLTsolver_method("new"), ...)
    .Object@pointer <- LDLTdense_new()
    .Object
})

new_laplacian = function(L, lambda) {
    l <- list(L=lambda*L)
    class(l) <- "laplacian"
    if(is.sparseMatrix(L)) {
        l$solver = new("LDLTsparse")
        class(l) = c("laplacian_sparse", class(l))
    } else {
        l$solver = new("LDLTdense")
        class(l) = c("laplacian_dense", class(l))
    }
    return(l)
}

compute <- function(l, A) {
  UseMethod("compute")
}

#' @export
compute.laplacian_sparse <- function(l, A) {
  LDLTsparse_compute(l$solver@pointer, A)
    
}

#' @export
compute.laplacian_dense <- function(l, A) {
  LDLTdense_compute(l$solver@pointer, A)
    
}

solve <- function(l, b) {
    UseMethod("solve")
}

#' @export
solve.laplacian_sparse <- function(l, b) {
  LDLTsparse_solve(l$solver@pointer, b)
}

#' @export
solve.laplacian_dense <- function(l, b) {
  LDLTdense_solve(l$solver@pointer, b)
}
