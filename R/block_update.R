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
  if (x$graph_laplacian$woodburry) {
    if (is.null(x$mu) || mu != x$mu) {
      x$mu <- mu
      x$inv <- lap_inv(x$graph_laplacian, 2*x$mu)
    }
    x$a <- as.matrix(
      x$inv %*% (a_grad + mu*block_project(x)$a)
    ) # if graph_laplacian is a sparse matrix
    # it necessary to cast the result as a dense matrix
  } else {
    A <- 2*mu*diag(length(x$a)) + x$graph_laplacian$L 
    b <- a_grad + mu*block_project(x)$a
    x$a <- rcpp_solve(A, b)
  }
  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}
