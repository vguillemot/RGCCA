block_update <- function(x, grad, mu) {
  UseMethod("block_update")
}

#' @export
block_update.block <- function(x, grad, mu) {
  x$a <- pm(t(x$x), grad, na.rm = x$na.rm)
  return(block_project(x))
}

#' @export
block_update.dual_block <- function(x, grad, mu) {
  x$alpha <- grad
  return(block_project(x))
}

#' @export
block_update.graphnet_block <- function(x, grad, mu) {
  x$a <- as.matrix(
    ginv(2*mu*diag(length(x$a)) + x$lambda * x$graph_laplacian) %*% (grad + mu*block_project(x)$a)
  ) # if graph_laplacian is a sparse matrix
  # it a necessary to cast the result as a dense matrix
  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}