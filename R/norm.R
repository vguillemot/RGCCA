#' L1-norm, L2-norm or LG-norm of a vector of numerics
#'
#' This documentation page groups the functions to compute the following norms:
#' L1-norm, L2-norm, and the grouped L2 (LG) norm.
#'
#' @param vec, vector of numeric value
#' @param grp, vector of numeric value
#'
#' @return the L1-norm, L2-norm or LG-norm of vec
#'
#' @examples
#' x <- c(-0.1, 1, 0.5)
#' g <- c(1, 1, 2)
#' normL1(x) # = 1.6
#' normL2(x) # ~= 1.12
#' normLG(x, g) # ~= 1.5
#' normalizeL2(x)
#' @name norm
NULL
#' @rdname norm
#' @export
normL1 <- function(vec) {
  return(sum(abs(vec)))
}
#' @rdname norm
#' @export
normL2 <- function(vec) {
  return(sqrt(sum(vec**2)))
}
#' @rdname norm
#' @export
normLG <- function(vec, grp) {
  return(sum(tapply(vec, grp, normL2)))
}


