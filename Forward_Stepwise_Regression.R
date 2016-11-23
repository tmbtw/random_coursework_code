######################################
############ Thai T. Pham ############
######## Stanford University #########
######################################

# NOTICE: The codes are not optimized for speed #


# This function takes as input the set of indices of 
# variables that have been chosen (in order of being chosen); 
# the current set W of orthonormal vectors; the current 
# value of RSS; the matrix X (columns mean centered); 
# the output vector y.
# 
# The function will return the updated set of orthonormal 
# vectors WW, newRSS, newIndexSet, and the coef. of the 
# added variable.
fwd.onestep <- function(y, X, W, RSS, indexSet)
{
  p <- dim(X)[2] 
  N <- dim(X)[1]
  q <- dim(W)[2]
  
  nextW <- NULL
  nextInd <- 0
  nextMax <- -1
  
  for (i in 1:p) {
    if (!(i %in% indexSet)) {
      R <- X[,i]
      for (j in 1:q) {
        R <- R - ((t(X[,i]) %*% W[,j])*W[,j])
      }
      w <- R / sqrt(sum(R^2)) 
      S <- (t(w) %*% y)^2
      
      if (S > nextMax) {
        nextMax <- S
        nextInd <- i
        nextW <- w
      }
    }
  }
  
  WW <- cbind(W, nextW)
  newRSS <- RSS - nextMax
  newIndexSet <- cbind(indexSet,nextInd)
  
  return(list(WW,newRSS,newIndexSet))
}

# This function takes as inputs the given data X, y
# It returns a sequence of indices of variables in 
# the order of being chosen; a sequence of of RSS
# values and coefficients, one at each step. 
fwd.stepwise <- function(y, X) 
{
  # Scale columns of X
  X <- scale(X, scale = FALSE)
  p <- dim(X)[2]
  N <- dim(X)[1]
  # initialize W, RSS, indexSet
  W <- matrix(1,N,1)
  W <- W/sqrt(sum(W^2))
  RSS <- sum((y-mean(y))^2)
  # we do not include the index corresponding with the intercept
  indexSet <- NULL 
  
  RSSSet <- c(RSS)
  coefMat <- matrix(0,p+1,p+1)
  coefMat[1,1] <- mean(y)
  
  while (length(indexSet) < p) {
    update <- fwd.onestep(y, X, W, RSS, indexSet)
    W <- update[[1]]
    RSS <- update[[2]]
    indexSet <- update[[3]]
    RSSSet <- cbind(RSSSet,RSS)
    
    currentX <- cbind(matrix(1,N,1),X[,indexSet])
    R <- t(W)%*%currentX
    currentCoefs <- solve(R)%*%t(W)%*%y
    
    l <- length(indexSet)
    coefMat[1:(l+1),l+1] <- currentCoefs
  }
  return(list(indexSet,RSSSet,coefMat))
}
