rgcca_inner_loop_Laplacian <- function(A, C, g, dg, tau = rep(1, length(A)),
                                       sparsity = rep(1, length(A)),
                                       lambda = rep(0, length(A)),
                                       graph_laplacians = NULL,
                                       verbose = FALSE, init = "svd", bias = TRUE,
                                       mu_init=1,
                                       tol_inner = 1e-04, tol_outer = 1e-08, na.rm = TRUE, n_iter_max = 1000) {
  if (!is.numeric(tau)) {
    # From Schafer and Strimmer, 2005
    tau <- vapply(A, tau.estimate, na.rm = na.rm, FUN.VALUE = 1.0)
  }
  
  # TODO: change this behaviour
  if (any(sparsity == 0)) {
    tau[which(sparsity == 0)] <- 0
    sparsity[which(sparsity == 0)] <- 1
  }
  
  ### Initialization
  block_objects <- lapply(seq_along(A), function(j) {
    create_block(A[[j]], j, bias, na.rm, tau[j], sparsity[j], lambda[j], 
                 graph_laplacians[[j]], tol_inner)
  })
  
  block_objects <- lapply(block_objects, block_init, init = init)
  Y <- do.call(cbind, lapply(block_objects, "[[", "Y"))
  N <- block_objects[[1]]$N
  
  iter_outer <- 1
  crit <- NULL
  crit_old <- sum(C * g(crossprod(Y) / N))
  a_old_inner <- a_old_outer <- lapply(block_objects, "[[", "a")
  mu <- mu_init
  
  repeat{  
    iter_inner <- 1
    crit <- NULL
    crit_old <- sum(C * g(crossprod(Y) / N))
    repeat {
      for (j in seq_along(A)) {
        # Compute grad
        grad <- Y %*% (C[j, ] * dg(crossprod(Y, Y[, j]) / N))
        block_objects[[j]] <- block_update(block_objects[[j]], grad, mu)
        Y[, j] <- block_objects[[j]]$Y
      }
      
      # Print out intermediate fit
      crit <- c(crit, sum(C * g(crossprod(Y) / N)))
      
      if (verbose) {
        cat(
          " iter_inner: ", formatC(iter_inner, width = 3, format = "d"),
          " Fit: ", formatC(crit[iter_inner], digits = 8, width = 10, format = "f"),
          " Dif: ", formatC(crit[iter_inner] - crit_old,
                            digits = 8, width = 10, format = "f"
          ), "\n"
        )
      }
      
      a <- lapply(block_objects, "[[", "a")
      stopping_criteria_inner <- crossprod(unlist(a, FALSE, FALSE) - unlist(a_old_inner, FALSE, FALSE))
      
      if (verbose) print(paste("stop inner:", stopping_criteria_inner))
      if (any(stopping_criteria_inner < tol_inner) || (iter_inner > n_iter_max)) {
        break
      }
      
      crit_old <- crit[iter_inner]
      a_old_inner <- a
      iter_inner <- iter_inner + 1
    }
    
    stopping_criteria_outer <- crossprod(unlist(a, FALSE, FALSE) - unlist(a_old_outer, FALSE, FALSE))

    if (any(stopping_criteria_outer < tol_outer) || (iter_outer > n_iter_max)) {
      #if (flag1) print("tolerance")
      #if (flag2) print("iterations")
      break
    }
    
    crit_old <- crit[iter_inner]
    a_old_outer <- a_old_inner <- a
    iter_outer <- iter_outer + 1
    mu <- 2*mu + 1
    if (verbose) print(paste("MU UPDATE:", mu))
  }
  
  if (iter_inner > n_iter_max) {
    warning(
      "The distance majorization algorithm inner loop for RGCCA did not converge after ", n_iter_max,
      " iterations."
    )
  }
  
  if (iter_outer > n_iter_max) {
    warning(
      "The distance majorization algorithm outer loop for RGCCA did not converge after ", n_iter_max,
      " iterations."
    )
  }
  
  if (verbose) {
    if (iter_outer <= n_iter_max) {
      message(
        "The RGCCA algorithm converged to a stationary point after ",
        iter_outer - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
  }
  
  # Post-process the resulting block-weight and block-component vectors
  ctrl <- all(g(-5:5) == g(5:-5))
  block_objects <- lapply(block_objects, block_postprocess, ctrl)
  a <- lapply(block_objects, "[[", "a")
  Y <- do.call(cbind, lapply(block_objects, "[[", "Y"))
  
  return(list(Y = Y, a = a, crit = crit, tau = tau))
}
