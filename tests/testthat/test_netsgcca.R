library(dplyr)
library(Matrix)
library(ggplot2)
library(RGCCA)

set.seed(565359)
A <- list(
  B1 = matrix(rnorm(10*5), 10, 5),
  B2 = matrix(rnorm(10*5), 10, 5),
  B3 = matrix(rnorm(10*5), 10, 5)
)

lambdas <- c(100, NA, NA)
LapList <- list(
  B1 = as(diag(5), "dgCMatrix"),
  B2 = NULL,
  B3 = NULL)

res.rgcca <- rgcca(blocks = A)
res.rgcca
res.netsgcca <- rgcca(
  blocks = A,
  lambda = lambdas,
  graph_laplacians = LapList,
  sparsity = c(1, 1, 1),
  ncomp = c(1, 1, 1),
  mu_init = 1, verbose = TRUE)
res.netsgcca

tibble(
  RGCCA = res.rgcca$Y[[1]][,1],
  netSGCCA = res.netsgcca$Y[[1]][,1]) %>%
  ggplot(aes(RGCCA, netSGCCA)) +
  geom_point() +
  coord_equal()

