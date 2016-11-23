#######################################
#            Thai T. Pham             #
#        Stanford University          #
#######################################

# NOTICE: The code is not optimized for speed #

# Self-Organizing Map (SOM) Algorithm

# find the neighbors of M[k, l, ]
# threshold r for defining a neighbor
neighbors <- function(M, k, l, r)
{
  q1 <- dim(M)[1]
  q2 <- dim(M)[2]
  
  result <- matrix(0, q1*q2, 2)
  
  for (i in 1:q1) {
    for (j in 1:q2) {
      if (sqrt((i - k)^2 + (j - l)^2) < r) {
        result[q1*(i - 1) + j, ] <- c(i, j)
      }
    }
  }
  return(result)
}

# SOM Algorithm
# MaxIter is set at high value, e.g. 5000
# q1, q2 are set at low value, e.g. 5 (so that q1 * q2 < dim(X)[1])
# X is the input matrix to be approximated
SOM <- function(X, MaxIter, q1, q2)
{
  alphas <- seq(1, 0, length = MaxIter)
  Rs <- seq(2, 1, length = MaxIter)
  set.seed(69)
  M <- array(rnorm(q1*q2*3), dim = c(q1,q2,3))
  
  for (iter in 1:MaxIter) {
    u <- sample(dim(X)[1])[1]
    x <- X[u, ]
    p <- closest_point(M, x)
    
    k <- p[[1]]
    l <- p[[2]]
    
    neibh_points <- neighbors(M, k, l, Rs[iter])
    
    for (t in 1:dim(neibh_points)[1]) {
      if (sum(neibh_points[t, ]) > 0) {
        i <- neibh_points[t, 1]
        j <- neibh_points[t, 2]
        
        # update
        M[i, j, ] <- M[i, j, ] + alphas[iter] * (x - M[i, j, ])
      }
    }
  }
  return(M)
}
