laplacian_filtering = function(L, threshold){
  eigenL                      = eigen(L, symmetric = T)
  dim_L                       = unique(dim(L))
  idx_keep                    = which(eigenL$values <= threshold)
  values                      = eigenL$values[idx_keep]
  vectors                     = eigenL$vectors[, idx_keep]
  L                           = as.matrix(vectors) %*% as.matrix(diag(values)) %*% t(as.matrix(vectors))
  L                           = (L + t(L))/2
  return(list(L = L, vectors = vectors, values = values, filtered = T))
}