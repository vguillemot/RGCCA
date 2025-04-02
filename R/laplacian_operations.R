lap_decomp <- function(L, lambda, k=NULL, woodburry=TRUE) {
  if (is.null(L)) return(list(NULL))
  if (is.null(k)) k = ncol(L) - 1

  L = lambda*L
  if (woodburry || k != ncol(L) - 1) {
      #Leig = eigen(L)
      Leig = rcpp_sparse_eigen(lambda*L, k)
      # ordering eigenvalues by magnitude
      #order_value = order(abs(Leig$values), decreasing=TRUE)
      #V = Leig$vectors[,order_value]
      #Delta = Leig$values[order_value]
      # Keeping the rank-k approximation
      #V = V[,1:k]
      #Delta = Delta[1:k]
      Delta = c(Leig$s)
      Delta_inv = 1/Delta
      V = Leig$v
      L = V %*% diag(Delta) %*% t(V)
  } else {
      V = NULL
      Delta_inv = NULL
  }

  return(list(list(V=V, Delta_inv=Delta_inv, L=L, woodburry=woodburry)))
}

lap_inv <- function(L, mu) {
  if (L$woodburry) {
    inv = (diag(nrow(L$V)) - L$V %*% diag(1/(mu*L$Delta_inv + 1)) %*% t(L$V))/mu
  } else {
    inv = solve(mu*diag(ncol(L$L)) + L$L)
  }
  return(inv)
}
