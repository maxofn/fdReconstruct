directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

library(doParallel)
library(FDFM)
library(foreach)

# Run simulation for alpha = 0.25 -----------------------------------------
set.seed(1)
alpha0 = 0.25
predres1A <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres2A <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres3A <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)

predres1B <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres2B <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres3B <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)

predres <- rbind(predres1A, predres2A, predres3A, predres1B, predres2B, predres3B)
write.csv(predres, file = paste0(getwd(), "/Res/Bands/pred_res_075.csv"))

# Run simulation for alpha = 0.10 -----------------------------------------
set.seed(1)
alpha0 = 0.10
predres1A <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres2A <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres3A <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)

predres1B <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres2B <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres3B <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)

predres <- rbind(predres1A, predres2A, predres3A, predres1B, predres2B, predres3B)
write.csv(predres, file = paste0(getwd(), "/Res/Bands/pred_res_090.csv"))

# Run simulation for alpha = 0.05 -----------------------------------------
set.seed(1)
alpha0 = 0.05
predres1A <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres2A <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)
predres3A <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'A', n.rep = 100)

predres1B <- SimPredint(T =  50, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres2B <- SimPredint(T = 100, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)
predres3B <- SimPredint(T = 200, N = 51, alpha = alpha0, type.miss = 'B', n.rep = 100)

predres <- rbind(predres1A, predres2A, predres3A, predres1B, predres2B, predres3B)
write.csv(predres, file = paste0(getwd(), "/Res/Bands/pred_res_095.csv"))

# Functions ---------------------------------------------------------------

SimPredint <- function(T, N, alpha, r.true = 10, r.max = 15, type.miss, n.rep = 1) {
  # Perform Simulation

  cluster <- makeCluster(4)
  registerDoParallel(cluster)

  res <- foreach(rep = 1:n.rep, .combine = 'cbind', .packages = 'FDFM') %dopar% {

    data.test <- GenObs(T = 100, N = N, r = r.true, complete = FALSE, type.miss = type.miss)
    data.train <- GenObs(T = T, N = N, r = r.true, complete = TRUE, type.miss = type.miss)
    data <- list()
    data$X <- rbind(data.train$X, data.test$X)
    data$Y <- rbind(data.train$Y, data.test$Y)
    data$Y.obs <- rbind(data.train$Y.obs, data.test$Y.obs)

    incompletely.obs <- which(rowSums(is.na(data$Y.obs)) > 0)

    # Reconstruction via factor models
    reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs,
                         method = "Kraus",
                         pred.band = TRUE,
                         p = 1 - alpha,
                         r.max = r.max)

    coverage <- 0
    for(t in 1:100) {
      O.t <- which(!is.na(data$Y.obs[T + t,]))
      if((max(data$X[T + t,-O.t] - reconst$U.hat[t,-O.t]) <= 0) &
         (max(data$X[T + t,-O.t] - reconst$L.hat[t,-O.t]) >= 0)) {
        coverage <- coverage + 1
      }
    }
    coverage/100
  }

  stopCluster(cluster)
  
  simres <- round(data.frame(MEAN = mean(res), SD = sd(res) * sqrt(99/100)), 4)
  colnames(simres) <- c("MEAN", "SD")
  rownames(simres) <- c(paste0(type.miss, T))
  
  return(simres)
}

