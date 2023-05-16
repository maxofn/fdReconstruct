# References:
#
# Yao, F., Müller, H. G., & Wang, J. L. (2005). Functional data analysis for sparse
# longitudinal data. Journal of the American statistical association, 100(470), 577-590.
#
# Kneip, A., & Liebl, D. (2020). On the optimal reconstruction of partially observed
# functional data. The Annals of Statistics, 48(3), 1692-1717.

directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

library(doParallel)
library(foreach)
library(FDFM)
library(ReconstPoFD)

set.seed(1)

# Run Simulation A ---------------------------------------------------------


SimReconstructionA(T =  50, N = 51, type.miss = 'A', n.rep = 100)
SimReconstructionA(T = 100, N = 51, type.miss = 'A', n.rep = 100)
SimReconstructionA(T = 200, N = 51, type.miss = 'A', n.rep = 100)

# Run Simulation B ---------------------------------------------------------

SimReconstructionB(T =  50, N = 51, type.miss = 'B', n.rep = 100)
SimReconstructionB(T = 100, N = 51, type.miss = 'B', n.rep = 100)
SimReconstructionB(T = 200, N = 51, type.miss = 'B', n.rep = 100)

# Functions ---------------------------------------------------------------

SimReconstructionA <- function(T, N, r.true = 10, r.max = 15, type.miss, n.rep = 1) {
  # Initialize variables

  MAE <- numeric(6)
  names(MAE) <- c("AFM", "ANo", "ANoCE", "AYes", "AYesCE", "PACE")
  
  cluster <- makeCluster(4)
  registerDoParallel(cluster)

  res <- foreach(rep = 1:n.rep, .combine = 'rbind',
                 .packages = c('FDFM', 'ReconstPoFD')) %dopar% {
  
    data.test <- GenObs(T = 100, N = N, r = r.true, complete = FALSE, type.miss = 'A')
    data.train <- GenObs(T = T, N = N, r = r.true, complete = TRUE, type.miss = type.miss)
    data <- list()
    data$X <- rbind(data.train$X, data.test$X)
    data$Y <- rbind(data.train$Y, data.test$Y)
    data$Y.obs <- rbind(data.train$Y.obs, data.test$Y.obs)
    
    incompletely.obs <- which(rowSums(is.na(data$Y.obs)) > 0)
    
    # Reconstruction via AFM --------------------------------------------------
    reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = "Kraus", r.max = r.max)
    r.hat <- reconst$r
    
    MAE["AFM"] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))
    
    # Reconstruction following Kneip & Liebl ----------------------------------
    grid <- seq(0, 1, length.out = N)
    Lu <- vector(mode = "list", length = T + 100)
    Ly <- vector(mode = "list", length = T + 100)
    for(t in 1:(T + 100)) {
     O.t <- which(!is.na(data$Y.obs[t,]))
     Lu[[t]] <- grid[O.t]
     Ly[[t]] <- data$Y.obs[t,O.t]
    }
    
    r.unique <- unique(r.hat)
    
    # ANo ---------------------------------------------------------------------
    X.hat <- matrix(NA, nrow = 100, ncol = N)
    
    for(r in r.unique) {
     idx <- which(r.hat == r)
     reconst <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'Error>=0_AlignNO',
                                      K = r,
                                      nRegGrid = N)
     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
    }
    
    MAE["ANo"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
    
    # ANoCE -------------------------------------------------------------------
    X.hat <- matrix(NA, nrow = 100, ncol = N)
    
    for(r in r.unique) {
     idx <- which(r.hat == r)
     reconst <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'Error>0_AlignNO_CEscores',
                                      K = r,
                                      nRegGrid = N)
     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
    }
    
    MAE["ANoCE"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
    
    # AYes --------------------------------------------------------------------
    X.hat <- matrix(NA, nrow = 100, ncol = N)
    
    for(r in r.unique) {
     idx <- which(r.hat == r)
     reconst <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'Error>0_AlignYES',
                                      K = r,
                                      nRegGrid = N)
     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
    }
    
    MAE["AYes"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
    
    # AYesCE ------------------------------------------------------------------
    X.hat <- matrix(NA, nrow = 100, ncol = N)
    
    for(r in r.unique) {
     idx <- which(r.hat == r)
     reconst <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'Error>0_AlignYES_CEscores',
                                      K = r,
                                      nRegGrid = N)
     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
    }
    
    MAE["AYesCE"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
    
    # Reconstruction following Yao, Mueller, and Wang -------------------------
    X.hat <- matrix(NA, nrow = 100, ncol = N)
    for(r in r.unique) {
     idx <- which(r.hat == r)
     reconst <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'PACE',
                                      K = r,
                                      nRegGrid = N)
     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
    }
    MAE["PACE"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
    MAE
  }
  
  stopCluster(cluster)
  
  filename <- paste0("T_", T, "_", type.miss)
  
  MAE <- apply(res, 2, mean)
  SD <- apply(res, 2, sd) * sqrt(99/100)
  
  simres <- rbind(MAE, SD)
  colnames(simres) <- names(MAE)
  rownames(simres) <- c("MAE", "SD")
  
  write.csv(round(simres, 4), file = paste0(getwd(), "/Res/Reconstruction/", filename, ".csv"))
  return()
}

SimReconstructionB <- function(T, N, r.true = 10, r.max = 15, type.miss, n.rep = 1) {
  # Initialize variables
  
  MAE <- numeric(2)
  names(MAE) <- c("AFM", "PACE")
  
  cluster <- makeCluster(4)
  registerDoParallel(cluster)
  
  res <- foreach(rep = 1:n.rep, .combine = 'rbind',
                 .packages = c('FDFM', 'ReconstPoFD')) %dopar% {
                   
                   data.test <- GenObs(T = 100, N = N, r = r.true, complete = FALSE, type.miss = 'B')
                   data.train <- GenObs(T = T, N = N, r = r.true, complete = TRUE, type.miss = type.miss)
                   data <- list()
                   data$X <- rbind(data.train$X, data.test$X)
                   data$Y <- rbind(data.train$Y, data.test$Y)
                   data$Y.obs <- rbind(data.train$Y.obs, data.test$Y.obs)
                   
                   incompletely.obs <- which(rowSums(is.na(data$Y.obs)) > 0)
                   
                   # Reconstruction via AFM --------------------------------------------------
                   reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = "Kraus", r.max = r.max)
                   r.hat <- reconst$r
                   
                   MAE["AFM"] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))
                   
                   # Reconstruction following Kneip & Liebl ----------------------------------
                   grid <- seq(0, 1, length.out = N)
                   Lu <- vector(mode = "list", length = T + 100)
                   Ly <- vector(mode = "list", length = T + 100)
                   for(t in 1:(T + 100)) {
                     O.t <- which(!is.na(data$Y.obs[t,]))
                     Lu[[t]] <- grid[O.t]
                     Ly[[t]] <- data$Y.obs[t,O.t]
                   }
                   
                   r.unique <- unique(r.hat)
                   
                   # Reconstruction following Yao, Mueller, and Wang -------------------------
                   X.hat <- matrix(NA, nrow = 100, ncol = N)
                   for(r in r.unique) {
                     idx <- which(r.hat == r)
                     reconst <- reconstructKneipLiebl(Ly, Lu,
                                                      method = 'PACE',
                                                      K = r,
                                                      nRegGrid = N)
                     X.hat[idx,] <- t(matrix(unlist(reconst[['Y_reconst_list']]), ncol = T + 100))[idx + T,]
                   }
                   MAE["PACE"] <- mean(apply(abs(X.hat - data$X[incompletely.obs,]),1, max))
                   MAE
                 }
  
  stopCluster(cluster)
  
  filename <- paste0("T_", T, "_", type.miss)
  
  MAE <- apply(res, 2, mean)
  SD <- apply(res, 2, sd) * sqrt(99/100)
  
  simres <- rbind(MAE, SD)
  colnames(simres) <- names(MAE)
  rownames(simres) <- c("MAE", "SD")
  
  write.csv(round(simres, 4), file = paste0(getwd(), "/Res/Reconstruction/", filename, ".csv"))
  return()
}
