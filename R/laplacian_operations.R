lap_decomp <- function(L, lambda, k=NULL, woodburry=TRUE) {
  if (is.null(L)) return(list(NULL))
  if (is.null(k)) k = ncol(L) - 1
  Leig = eigen(lambda*L)
  
  # ordering eigenvalues by magnitude
  order_value = order(abs(Leig$values), decreasing=TRUE)
  V = Leig$vectors[,order_value]
  Delta = Leig$values[order_value]

  # Keeping the rank-k approximation
  V = V[,1:k]
  Delta = Delta[1:k]
  #Lk = V %*% diag(Delta) %*% t(V)

  return(list(list(V=V, Delta_inv=1/Delta, L=L, woodburry=woodburry)))
}

lap_inv <- function(L, mu) {
  if (L$woodburry) {
    (diag(nrow(L$V)) - L$V %*% diag(1/(mu*L$Delta_inv + 1)) %*% t(L$V))/mu
  } else {
    solve(mu*diag(ncol(L$L)) + L$L)
  }
}