#######################################
############ Thai T. Pham #############
####### thaipham@stanford.edu #########
#######################################

# NOTICE: The codes are not optimized for speed #

# FastICA Algorithm
# MaxIter is set at high value, e.g. 20,000
fastICA <- functioN(X, MaxIter)
{
  # Centering X
  X <- scale(X, center = TRUE, scale = FALSE)
  
  # Whitening X
  Sig <- cov(X)
  s <- svd(Sig)
  K <- s$u %*% diag((sqrt(s$d))^(-1)) %*% t(s$v)
  X <- X %*% K
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  
  # Find W so as to recover S = XW
  W <- matrix(runif(p^2), p, p)
  e <- matrix(1, N, 1)
  
  for (iter in 1:MaxIter)
  {
    for (j in 1:p) {
      W[, j] <- 1/N * t(X) %*% tanh(X %*% W[, j]) - 
        1/N * (t(e) %*% (1 - tanh(X %*% W[, j])^2))[1] * W[, j]
    }
    # Orthogonalize W
    svd <- svd(W)
    W <- svd$u %*% t(svd$v)
  }
  S <- X %*% W
  return(S)
}
