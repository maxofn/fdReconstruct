directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

library(doParallel)
library(foreach)
library(fdReconstruct)
library(fdapace)
library(ReconstPoFD)

set.seed(1)

# Run Simulation A ---------------------------------------------------------

predres1A_exp_01   <- SimReconstruction(T =  50, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2A_exp_01   <- SimReconstruction(T = 100, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1A_exp_005  <- SimReconstruction(T =  50, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2A_exp_005  <- SimReconstruction(T = 100, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1A_poly_01  <- SimReconstruction(T =  50, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2A_poly_01  <- SimReconstruction(T = 100, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1A_poly_005 <- SimReconstruction(T =  50, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2A_poly_005 <- SimReconstruction(T = 100, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres <- rbind(predres1A_exp_01, predres2A_exp_01,
                 predres1A_exp_005, predres2A_exp_005,
                 predres1A_poly_01, predres2A_poly_01,
                 predres1A_poly_005, predres2A_poly_005)

write.csv(predres, file = paste0(getwd(), "/Results/Reconstruction/SettingA.csv"))

# Run Simulation B ---------------------------------------------------------

predres1B_exp_01   <- SimReconstruction(T =  50, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2B_exp_01   <- SimReconstruction(T = 100, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1B_exp_005  <- SimReconstruction(T =  50, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2B_exp_005  <- SimReconstruction(T = 100, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1B_poly_01  <- SimReconstruction(T =  50, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2B_poly_01  <- SimReconstruction(T = 100, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1B_poly_005 <- SimReconstruction(T =  50, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2B_poly_005 <- SimReconstruction(T = 100, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres <- rbind(predres1B_exp_01, predres2B_exp_01,
                 predres1B_exp_005, predres2B_exp_005,
                 predres1B_poly_01, predres2B_poly_01,
                 predres1B_poly_005, predres2B_poly_005)

write.csv(predres, file = paste0(getwd(), "/Results/Reconstruction/SettingB.csv"))

# Functions ---------------------------------------------------------------

SimReconstruction <- function(T, type.miss, ev, eps.sd, n.rep = 1) {
  # Initialize variables

  MAE <- numeric(5)
  names(MAE) <- c("Uni", "Mult", "PACE", "KL", "FLM")

  cluster <- makeCluster(4)
  registerDoParallel(cluster)

  res <- foreach(rep = 1:n.rep, .combine = 'rbind',
                 .packages = c('fdReconstruct', 'fdapace', 'ReconstPoFD')) %dopar% {

    data_test <-  GenObs(T = 50, type.miss = type.miss, ev = ev,
                         eps.sd = eps.sd, complete = FALSE)
    data_train <- GenObs(T = T, type.miss = type.miss, ev = ev,
                         eps.sd = eps.sd, complete = TRUE)
    
    Y0.obs <- rbind(data_test$Y0.obs, data_train$Y0.obs)
    Y1.obs <- rbind(data_test$Y1.obs, data_train$Y1.obs)
    
    Y.obs <- cbind(Y0.obs, Y1.obs)
    incompletely.obs <- which(rowSums(is.na(Y.obs)) > 0)
    completely.obs <- which(rowSums(is.na(Y.obs)) == 0)
    Tc <- length(completely.obs)
    Tm <- length(incompletely.obs)

    # Reconstruction using univariate factor model ----------------------------
    reconst_uni  <- fdReconstruct(Y0.obs, T.set = incompletely.obs)

    MAE["Uni"]  <- mean(apply(abs(reconst_uni$X.hat[incompletely.obs,] -
                                    data_test$X0[incompletely.obs,]),1, max))

    # Reconstruction using multivariate factor model --------------------------

    reconst_mult <- fdReconstruct(Y0.obs, Y1.obs, T.set = incompletely.obs)

    MAE["Mult"] <- mean(apply(abs(reconst_mult$X.hat[incompletely.obs,] -
                                    data_test$X0[incompletely.obs,]),1, max))

    # Reconstruction following Yao, Müller & Wang (2005a) ---------------------
    N <- dim(Y0.obs)[2]
    grid <- seq(0, 1, length.out = N)
    Lu <- vector(mode = "list", length = T + 50)
    Ly <- vector(mode = "list", length = T + 50)
    for(t in 1:(T + 50)) {
      O.t <- which(!is.na(Y0.obs[t,]))
      Lu[[t]] <- grid[O.t]
      Ly[[t]] <- Y0.obs[t,O.t]
    }
    
    fpca <- FPCA(Ly, Lu)
    X.hat.PACE <- fitted(fpca)[incompletely.obs,]

    MAE["PACE"] <- mean(apply(abs(X.hat.PACE -
                                    data_test$X0[incompletely.obs,]),1, max))


    # Reconstruction following Kneip & Liebl ----------------------------------
    reconst_KL <- reconstructKneipLiebl(Ly, Lu,
                                        method = 'Error>0_AlignYES_CEscores',
                                        nRegGrid = N)
    X.hat.KL <- t(matrix(unlist(reconst_KL[['Y_reconst_list']]),
                         ncol = T + 50))[incompletely.obs,]

    MAE["KL"] <- mean(apply(abs(X.hat.KL - data_test$X0[incompletely.obs,]),1, max))

    # Reconstruction following Yao, Müller & Wang (2005b) ---------------------
    Lt <- list()
    Ly <- list()
    Lx <- list()

    for(t in 1:Tc) {
      Lt[[t]] <- grid
      Lx[[t]] <- Y1.obs[completely.obs[t],]
      Ly[[t]] <- Y0.obs[completely.obs[t],]
    }
    X <- list(X = list(Ly = Lx, Lt = Lt))
    Y <- list(Ly = Ly, Lt = Lt)

    Lt <- list()
    Lx.test <- list()
    for(t in 1:Tm) {
      Lt[[t]] <- grid
      Lx.test[[t]] <- Y1.obs[incompletely.obs[t],]
    }
    X.test <- list(X = list(Ly = Lx.test, Lt = Lt))

    reconst_FLM <- FLM1(Y, X, X.test)
    X.hat.FLM <- matrix(unlist(reconst_FLM[['yPred']]), nrow = Tm)

    MAE["FLM"] <- mean(apply(abs(X.hat.FLM - data_test$X0[incompletely.obs,]),1, max))
    MAE
  }

  stopCluster(cluster)

  MAE <- apply(res, 2, mean)
  SD <- apply(res, 2, sd) * sqrt((n.rep - 1)/n.rep)

  simres <- data.frame(matrix(c(rbind(round(MAE, 2), round(SD, 2))), nrow = 1))
  colnames(simres)[2 * 1:5 - 1] <- names(MAE)
  colnames(simres)[2 * 1:5] <- rep("SD", 5)
  rownames(simres) <- c(paste(type.miss, ev, eps.sd, T))
  
  return(simres)
}
