directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

library(doParallel)
library(foreach)
library(FDFM)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(latex2exp)

set.seed(1)

# Run Simulation ----------------------------------------------------------

SimTuning(T =  50, N = 51, type.miss = 'A', n.rep = 100)
SimTuning(T = 100, N = 51, type.miss = 'A', n.rep = 100)
SimTuning(T = 200, N = 51, type.miss = 'A', n.rep = 100)

SimTuning(T =  50, N = 51, type.miss = 'B', n.rep = 100)
SimTuning(T = 100, N = 51, type.miss = 'B', n.rep = 100)
SimTuning(T = 200, N = 51, type.miss = 'B', n.rep = 100)

filepath <- paste0(getwd(), "/Res/Estimation r/")
plotres <- CreateFigures(filepath)

pdf(paste0(filepath, "SimTuning.pdf"), height = 18,
    width = 20)
plotres
dev.off()

# Functions ---------------------------------------------------------------

SimTuning <- function(T, N, r.true = 10, r.max = 15, type.miss, n.rep = 10) {
  # Initialize variables
  r.vec <- 1:r.max

  error <- numeric(3 + r.max)
  names(error) <- c("BaiNg", "Onatski", "Kraus", paste0("r.hat", 1:r.max))
  
  # Perform Simulation
  cluster <- makeCluster(4)
  registerDoParallel(cluster)
  
  res <- foreach(rep = 1:n.rep, .combine = 'rbind',
                 .packages = 'FDFM') %dopar% {
    
    data.test <- GenObs(T = 100, N = N, r = r.true, complete = FALSE, type.miss = type.miss)
    data.train <- GenObs(T = T, N = N, r = r.true, complete = TRUE, type.miss = type.miss)
    data <- list()
    data$X <- rbind(data.train$X, data.test$X)
    data$Y <- rbind(data.train$Y, data.test$Y)
    data$Y.obs <- rbind(data.train$Y.obs, data.test$Y.obs)

    incompletely.obs <- which(rowSums(is.na(data$Y.obs)) > 0)

    reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = "BaiNg", r.max = r.max)
    error["BaiNg"] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))
    
    reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = "Kraus", r.max = r.max)
    error["Kraus"] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))
    
    reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = "Onatski", r.max = r.max)
    error["Onatski"] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))

    for(r in 1:r.max) {
      reconst <- ReconstFD(data$Y.obs, T.set = incompletely.obs, method = 'rGiven',
                           r = rep(r,100), r.max = r.max)
      error[3 + r] <- mean(apply(abs(reconst$X.hat - data$X[incompletely.obs,]),1, max))
    }
    
    error
  }

  stopCluster(cluster)
  
  # Save results
  filename <- paste0("T_", T, "_", type.miss)
  
  MAE <- apply(res, 2, mean)
  SD <- apply(res, 2, sd) * sqrt(99/100)
  
  simres <- rbind(MAE, SD)
  colnames(simres) <- names(error)
  rownames(simres) <- c("MAE", "SD")
  
  write.csv(round(simres, 3), file = paste0(getwd(), "/Res/Estimation r/", filename, ".csv"))
  return()
}

CreateFigures <- function(filepath) {
  r.max <- 15
  plot.list <- list()
  labels <- c("Bai & Ng", "Kraus", "Onatski", TeX(r'($fixed\,\hat{r}_o$)'))
  col <- c(brewer.pal(3, name="Paired"), "#000000")
  col <- c(col[2], col[1], col[3], col[4])
  counter <- 1
  for(type.miss in c("A", "B")) {
    for(T in c(50, 100, 200)) {
      
      filename <- paste0("T_", T, "_", type.miss)
      error <- read.csv(paste0(filepath, filename, ".csv"))[1,-1]
      
      main <- bquote(T[c] == .(T)*', DGP'*.(type.miss))
      
      df <- data.frame("method" = rep("BN", r.max), "y" = error$BaiN)
      df <- rbind(df, data.frame("method" = rep("K", r.max), "y" = error$Kraus))
      df <- rbind(df, data.frame("method" = rep("O", r.max), "y" = error$Onatski))
      df2 <- data.frame("method" = rep("fixed", r.max), "y" = t(unname(error[4:18])))
      colnames(df2)[2] <- "y"
      df <- rbind(df, df2)
      df$x <- rep(1:r.max, 4)
      
      plot.tuning <- ggplot2::ggplot(data=df, aes(x = x, y = y, col = method, 
                                                  linetype = method, size = method)) +
        ggtitle(main) +
        geom_line() +
        xlab(TeX(r'($\hat{r}_o$)')) +
        ylab("") +
        scale_x_continuous(breaks = seq(1, r.max, length.out = floor(r.max/2) + 1)) +
        scale_y_continuous(limits = c(0.30, 0.86)) +
        theme_minimal(base_size = 32) +
        scale_linetype_manual(breaks = c("BN", "K", "O", "fixed"),
                              labels = labels,
                              values = c(5, 2, 3, 1)) +
        scale_size_manual(breaks = c("BN", "K", "O", "fixed"),
                          labels = labels,
                          values = c(4, 4, 4, 1)) +
        scale_color_manual(breaks = c("BN", "K", "O", "fixed"),
                           labels = labels,
                           values = col) +
        theme(legend.key.width = unit(3.5,"cm"))
      plot.list[[counter]] <- plot.tuning
      counter <- counter + 1
    }
  }
  
  plot <- ggarrange(plot.list[[1]], plot.list[[2]], plot.list[[3]],
                    plot.list[[4]], plot.list[[5]], plot.list[[6]],
                    ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
  return(plot)
}
