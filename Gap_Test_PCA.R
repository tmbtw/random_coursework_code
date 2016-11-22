#######################################
############ Thai T. Pham #############
####### thaipham@stanford.edu #########
#######################################

# NOTICE: The codes are not optimized for speed #

# Gap-style Test to determine the optimal number of
# components of PCA to use
gap_test_pca <- function(Data, B)
{
  # centering and scaling the data to feed in PCA
  Data <- scale(Data, center = TRUE, scale = TRUE)
  p <- dim(Data)[2]
  
  # obtain principal components
  s <- svd(Data)
  T <- s$u %*% diag(s$d)
  
  # Calculate the variation explained by the principal components
  cov <- var(T)
  d_cov <- diag(cov)
  
  Var_PC <- matrix(0, 1, p)
  for (k in 1:p) {
    Var_PC[k] <- log(sum(d_cov[1:k]) / sum(d_cov)) 
  }
  
  # Calculate the MC - variation explained by the principal components
  # together with their Bootstrapped standard error
  Var_MC_PC <- matrix(0, B, p)
  for (j in 1:B) {
    cur_dat <- gen.MC.data(10*j, Concent_Data)
    cur_dat <- apply(cur_dat, 2, function(x){x - mean(x)})
    cur_s <- svd(cur_dat)
    cur_U <- cur_s$u
    cur_D <- diag(cur_s$d)
    cur_V <- cur_s$v
    
    cur_T <- cur_U %*% cur_D
    
    cur_cov <- var(cur_T)
    cur_d_cov <- diag(cur_cov)
    for (k in 1:p) {
      Var_MC_PC[j, k] <- log(sum(cur_d_cov[1:k]) / sum(cur_d_cov))
    }
  }
  MC_PC <- colMeans(Var_MC_PC)
  
  # Gap statistic
  Gap <- MC_PC - Var_PC
  
  ### Calculate the standard error for the estimate
  s_vals <- matrix(0, 1, 8)
  for (k in 1:8) {
    s_vals[k] <- sqrt(1 + 1/B) * sd(Var_MC_PC[, k]) 
  }
  
  # Determine the optimal number of components
  opt_k <- 0
  for (k in 1:7) {
    if (Gap[k] <= Gap[k + 1] - s_vals[k + 1]) {
      opt_k <- k
      break
    }
  }
  return(opt_k)
}