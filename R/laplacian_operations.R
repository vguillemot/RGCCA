# Defining a class to wrap around the pointer to LDLTsolver
# boilerplate taken from the book "seamless R and C++ integration..."
setClass("LDLTsolver", slots=c(pointer="externalptr"))

LDLTsolver_method <- function(name) {
    paste("LDLTsolver", name, sep="_")
}

setMethod("$", "LDLTsolver", function(x, name) {
    function(...) {
       print(name)
      .Call(LDLTsolver_method(name), x@pointer, ...)
    }
})

setMethod("initialize", "LDLTsolver", function(.Object, ...) {
    #.Object@pointer <- .Call(LDLTsolver_method("new"), ...)
    .Object@pointer <- LDLTsolver_new()
    .Object
})

create_laplacian = function(L, lambda, ranks=NULL) {
    structure(list(L=lambda*L, solver=new("LDLTsolver")),
        class="laplacian")
}
