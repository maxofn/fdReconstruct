directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

library(doParallel)
library(fdReconstruct)
library(foreach)

set.seed(1)

# Run simulation for alpha = 0.05  -----------------------------------------

alpha = 0.05

predres1A_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2A_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1A_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2A_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1A_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2A_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1A_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2A_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres1B_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2B_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1B_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2B_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1B_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2B_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1B_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2B_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres <- rbind(predres1A_exp_01, predres2A_exp_01,
                 predres1A_exp_005, predres2A_exp_005,
                 predres1A_poly_01, predres2A_poly_01,
                 predres1A_poly_005, predres2A_poly_005,
                 predres1B_exp_01, predres2B_exp_01,
                 predres1B_exp_005, predres2B_exp_005,
                 predres1B_poly_01, predres2B_poly_01,
                 predres1B_poly_005, predres2B_poly_005)

write.csv(predres, file = paste0(getwd(), "/Results/Bands/pred_res_095.csv"))

# Run simulation for alpha = 0.10  -----------------------------------------

alpha = 0.10

predres1A_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2A_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1A_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2A_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1A_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2A_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1A_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2A_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres1B_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2B_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1B_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2B_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1B_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2B_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1B_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2B_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres <- rbind(predres1A_exp_01, predres2A_exp_01,
                 predres1A_exp_005, predres2A_exp_005,
                 predres1A_poly_01, predres2A_poly_01,
                 predres1A_poly_005, predres2A_poly_005,
                 predres1B_exp_01, predres2B_exp_01,
                 predres1B_exp_005, predres2B_exp_005,
                 predres1B_poly_01, predres2B_poly_01,
                 predres1B_poly_005, predres2B_poly_005)

write.csv(predres, file = paste0(getwd(), "/Results/Bands/pred_res_090.csv"))

# Run simulation for alpha = 0.25  -----------------------------------------

alpha = 0.25

predres1A_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2A_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1A_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2A_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1A_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2A_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1A_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2A_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'A', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres1B_exp_01   <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres2B_exp_01   <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =   0.1, n.rep = 100)
predres1B_exp_005  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres2B_exp_005  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =   "exp", eps.sd =  0.05, n.rep = 100)
predres1B_poly_01  <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres2B_poly_01  <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =   0.1, n.rep = 100)
predres1B_poly_005 <- SimPredint(T =  50, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)
predres2B_poly_005 <- SimPredint(T = 100, alpha = alpha, type.miss = 'B', ev =  "poly", eps.sd =  0.05, n.rep = 100)

predres <- rbind(predres1A_exp_01, predres2A_exp_01,
                 predres1A_exp_005, predres2A_exp_005,
                 predres1A_poly_01, predres2A_poly_01,
                 predres1A_poly_005, predres2A_poly_005,
                 predres1B_exp_01, predres2B_exp_01,
                 predres1B_exp_005, predres2B_exp_005,
                 predres1B_poly_01, predres2B_poly_01,
                 predres1B_poly_005, predres2B_poly_005)

write.csv(predres, file = paste0(getwd(), "/Results/Bands/pred_res_075.csv"))

# Functions ---------------------------------------------------------------

SimPredint <- function(T, alpha, type.miss, ev, eps.sd, n.rep = 1) {
  # Perform Simulation

  cluster <- makeCluster(detectCores() - 1)
  registerDoParallel(cluster)

  res <- foreach(rep = 1:n.rep, .combine = 'cbind', .packages = 'fdReconstruct') %dopar% {

    # Generate test data
    data_test <-  GenObs(T = 100, type.miss = type.miss, ev = ev,
                         eps.sd = eps.sd, complete = FALSE)
    data_train <- GenObs(T = T, type.miss = type.miss, ev = ev,
                         eps.sd = eps.sd, complete = TRUE)

    Y0.obs <- rbind(data_test$Y0.obs, data_train$Y0.obs)
    Y1.obs <- rbind(data_test$Y1.obs, data_train$Y1.obs)

    incompletely.obs <- which(rowSums(is.na(Y0.obs)) > 0)

    # Reconstruction via factor models
    reconst <- fdReconstruct(Y0.obs, Y1.obs, T.set = incompletely.obs,
                             pred.band = TRUE, p = 1 - alpha)

    coverage <- 0
    for(t in 1:length(incompletely.obs)) {
      O.t <- which(!is.na(Y0.obs[incompletely.obs[t],]))
      if((max(data_test$X0[incompletely.obs[t],-O.t] - reconst$U.hat[t,-O.t]) <= 0) &
         (min(data_test$X0[incompletely.obs[t],-O.t] - reconst$L.hat[t,-O.t]) >= 0)) {
        coverage <- coverage + 1
      }
    }
    coverage/length(incompletely.obs)
  }

  stopCluster(cluster)

  simres <- round(data.frame(MEAN = mean(res), SD = sd(res) * sqrt((n.rep - 1)/n.rep)), 2)
  colnames(simres) <- c("MEAN", "SD")
  rownames(simres) <- c(paste(type.miss, ev, eps.sd, T))

  return(simres)
}

