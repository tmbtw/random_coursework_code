#######################################
#            Thai T. Pham             #
#        Stanford University          #
#######################################

# NOTICE: The code is not optimized for speed #

# Convex matrix approximations 
# This is the SOFT-IMPUTE algorithm in Mazunder, Hastie, and Tibs (2010).

# approximation for a given lmda
nuclear.norm.approx <- function(lmda, X) 
{ 
  s <- svd(X)
  D_lmda <- 1:length(s$d)
  for (i in 1:length(s$d)) 
  {
    D_lmda[i] <- max(s$d[i] - lmda/2, 0)
  }
  Z <- s$u %*% diag(D_lmda) %*% t(s$v)
  return(Z)
}

# The main soft_impute algorithm
# input: a missing-value matrix X
# lmda_set (which is usually 20:1)
# eps is the threshold of an extremely small value ~ 10^(-18)
# MaxIter is a large number ~ 10^9
#
# output: an approximation for X, one for each value of lmda 
Soft_Impute <- function(X, lmda_set, eps, MaxIter)
{
  O <- !is.na(X) # Observed values of X
  M <- is.na(X) # Missing values of X
  P_O_X <- matrix(0, dim(X)[1], dim(X)[2])
  P_O_X[O] <- X[O]
  
  # Initialize the approximated matrix Z
  Z_old <- matrix(0, dim(X)[1], dim(X)[2]) 
  
  # Convergence process for the given set of lmda's
  optimal_SVD_set <- list()
  
  for (i in 1:length(lmda_set))
  {
    lmda <- lmda_set[i]
    for (iter in 1:MaxIter)
    {
      P_O_perp_Z_old <- matrix(0, dim(Z_old)[1], dim(Z_old)[2])
      P_O_perp_Z_old[M] <- Z_old[M]
      Z_new <- nuclear.norm.approx(lmda, P_O_X + P_O_perp_Z_old) 
      
      if ((norm(Z_new - Z_old, "F")^2) / (norm(Z_old, "F")^2) < eps) {
        break
      }
      Z_old <- Z_new
    }
    optimal_SVD_set[[i]] <- Z_new
  }
  return(optimal_SVD_set)
}
