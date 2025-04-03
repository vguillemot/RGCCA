block_update <- function(x, grad, mu) {
  UseMethod("block_update")
}

#' @export
block_update.block <- function(x, grad, mu=NULL) {
  x$a <- pm(t(x$x), grad, na.rm = x$na.rm)
  return(block_project(x))
}

#' @export
block_update.dual_block <- function(x, grad, mu=NULL) {
  x$alpha <- grad
  return(block_project(x))
}

#' @export
block_update.graphnet_block <- function(x, grad, mu) {
  a_grad = pm(t(x$x), grad, na.rm = x$na.rm)
    if (is.null(x$mu) || mu != x$mu) {
      x$mu = mu
      #x$graph_laplacian$LDLTsolver$compute(2*mu*diag(length(x$a)) + x$graph_laplacian$L)
      #.Call("LDLTsolver_compute", x$graph_laplacian$solver@pointer, 2*mu*diag(length(x$a)) + x$graph_laplacian$L)
      LDLTsolver_compute(x$graph_laplacian$solver@pointer, 2*mu*diag(length(x$a)) + x$graph_laplacian$L)
    }
  #x$a <- x$graph_laplacian$LDLTsolver$solve(a_grad + mu*block_projet(x)$a)
  #.Call("LDLTsolver_solve", x$graph_laplacian$solver@pointer, a_grad + mu*block_project(x))
  x$a <- LDLTsolver_solve(x$graph_laplacian$solver@pointer, a_grad + mu*block_project(x)$a)
  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}
