directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)
library(fdReconstruct)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(latex2exp)
require(reshape2)

source("./reconstruct_EJS.R")

# Prepare data ------------------------------------------------------------

data <- temp_graz
data$time <- as.POSIXct(data$time)
N <- 48
T <- dim(data)[1]/48
Y0.obs <- matrix(data$temp_east , nrow = T, ncol = N, byrow = TRUE)
Y1.obs <- matrix(data$temp_west, nrow = T, ncol = N, byrow = TRUE)

# Analysis ----------------------------------------------------------------

# How many days are completely observed?
is.complete   <- which(rowSums(is.na(Y0.obs))==0)
is.incomplete <- which(rowSums(is.na(Y0.obs))>0)
length(is.complete) # 66
length(is.incomplete) # 10

# Mean temperatures
mean(Y0.obs, na.rm = TRUE) # 21.6
mean(Y1.obs, na.rm = TRUE) # 20.6
cor(data$temp_east, data$temp_west, use = "complete") # 0.90

# Plot data ---------------------------------------------------------------

plot(data$temp_east, type = "l")
lines(data$temp_south, lty = 2)
matplot(t(Y1.obs), type = "l")

# ggplot of data ----------------------------------------------------------

time <- seq(
  from=as.POSIXct("2022-07-01 0:00", tz="UTC"),
  to=as.POSIXct("2022-07-01 23:30", tz="UTC"),
  by="30 min"
)

colours <- brewer.pal(10, name="Paired")

# Temperature east

Y.miss.df <- melt(t(Y0.obs[is.incomplete,]))
Y.miss.df$Var1 <- rep(time, length(is.incomplete))
Y.comp.df <- melt(t(Y0.obs[is.complete,]))
Y.comp.df$Var1 <- rep(time, length(is.complete))

plot.sample0 <- ggplot2::ggplot() +
  geom_line(data = Y.comp.df, aes(x = Var1, y = value, group = Var2),
            col = "gray70", linewidth = 0.3) +
  geom_line(data = Y.miss.df, aes(x = Var1, y = value, group = Var2,
                                  col = as.factor(Var2)), linewidth = 1) +
  scale_color_manual(values = colours) +
  scale_x_datetime(date_breaks = "360 min",
                   date_labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  xlab("time") +
  ylab("temperature (°C)") +
  theme_classic(base_size = 11) +
  theme(legend.position="none") +
  ggtitle("East")
plot.sample0

# Temperature west

Y.comp.df <- melt(t(Y1.obs))
Y.comp.df$Var1 <- rep(time, T)

plot.sample1 <- ggplot2::ggplot() +
  geom_line(data = Y.comp.df, aes(x = Var1, y = value, group = Var2),
            col = "gray70", linewidth = 0.3) +
  scale_x_datetime(date_breaks = "360 min",
                   date_labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  xlab("time") +
  ylab("temperature (°C)") +
  theme_classic(base_size = 11) +
  theme(legend.position="none") +
  ggtitle("West")
plot.sample1

pdf(paste0(getwd(), "/Results/Real data/", "sample.pdf"), height = 3,
    width = 8.27)
do.call("ggarrange", c(list(plot.sample0, plot.sample1), nrow = 1, ncol = 2))
dev.off()


# Reconstruct data --------------------------------------------------------

reconst_mult <- fdReconstruct(Y0.obs = Y0.obs, Y1.obs = Y1.obs, T.set = is.incomplete,
                              pred.band = TRUE, p = 0.95, r.max = 30)
reconst_uni  <- fdReconstruct(Y0.obs = Y0.obs, T.set = is.incomplete,
                              pred.band = TRUE, p = 0.95, r.max = 30)

reconst_mult$r # [1] 18  8 18 18 15 24 19 25 18 16

# Plot reconstruction of partially observed curves ------------------------

time <- seq(
  from=as.POSIXct("2022-07-01 0:00", tz="UTC"),
  to=as.POSIXct("2022-07-01 23:30", tz="UTC"),
  by="30 min"
)

colours <- brewer.pal(10, name = "Paired")

t.count <- 1
list <- list()
for(t.idx in 1:10) {
  X.hat.mult <- data.frame(value = reconst_mult$X.hat[t.idx,], Var1 = time)
  X.hat.uni  <- data.frame(value =  reconst_uni$X.hat[t.idx,], Var1 = time)
  Y.obs.df<- data.frame(value = Y0.obs[is.incomplete[t.idx],], Var1 = time)
  pred.df <- data.frame(u.df = as.vector(reconst_mult$U.hat[t.idx,]),
                        l.df = as.vector(reconst_mult$L.hat[t.idx,]),
                        Var1 = time)

  plot.new <- ggplot2::ggplot() +
    geom_line(data = X.hat.uni , aes(x = Var1, y = value), linewidth = 0.5,
              linetype = "dotted") +
    geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
                fill = "gray70", alpha = 0.3) +
    geom_point(data = Y.obs.df, aes(x = Var1, y = value),
               size = 3, alpha = 0.5, col = colours[t.count]) +
    geom_line(data = X.hat.mult, aes(x = Var1, y = value),
              col = colours[t.count], linewidth = 1) +
    scale_x_datetime(date_breaks = "360 min",
       date_labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
    xlab("time") +
    ylab("temperature (°C)") +
    theme_classic(base_size = 11) +
    theme(legend.position="none") +
    ggtitle(format(data$time[(is.incomplete[t.idx] - 1) * N + 1], "%Y-%m-%d"))
  list[[t.count]] <- plot.new
  t.count <- t.count + 1
}

pdf(paste0(getwd(), "/Results/Real data/", "weather_reconstructions.pdf"),
    height = 11.69, width = 8.27)
do.call("ggarrange", c(list, nrow = 5, ncol = 2))
dev.off()

pdf(paste0(getwd(), "/Results/Real data/", "weather_examples.pdf"), height = 3,
    width = 8.27)
do.call("ggarrange", c(list(list[[5]], list[[6]]), nrow = 1, ncol = 2))
dev.off()


# Real data simulation study ----------------------------------------------

set.seed(123)

library(fdapace)
library(ReconstPoFD)
library(doParallel)
library(foreach)

cluster <- makeCluster(4)
registerDoParallel(cluster)

res <- foreach(rep = 1:100, .combine = 'rbind',
               .packages = c('fdReconstruct', 'fdapace', 'ReconstPoFD',
                             'fdaPOIFD')) %dopar% {

  T <- 60
  idx.sample <- sample(is.complete, size = T)
  
  Y0.sim <- Y0.obs[idx.sample,]
  Y1.sim <- Y1.obs[idx.sample,]
  X0.sim <- Y0.sim
  X1.sim <- Y1.sim
  
  sim.complete <- 1:50
  sim.incomplete <- 51:60
  
  for(t in sim.incomplete) {
    D <- sample(floor(1 * N/2):floor(3 * N/4), size = 1)
    Y0.sim[t, D:N] <- NA
  }
  
  MAE <- numeric(6)
  names(MAE) <- c("UNI", "MULT", "PACE", "AYESCE", "DEPTH", "FLM")
  
  # Reconstruction using univariate factor model ----------------------------
  reconst_uni  <- fdReconstruct(Y0.sim, T.set = sim.incomplete)
  
  MAE["UNI"]  <- mean(apply(abs(reconst_uni$X.hat - X0.sim[sim.incomplete,]),1, max))
  
  # Reconstruction using multivariate factor model --------------------------
  
  reconst_mult <- fdReconstruct(Y0.sim, Y1.sim, T.set = sim.incomplete)
  
  MAE["MULT"] <- mean(apply(abs(reconst_mult$X.hat - X0.sim[sim.incomplete,]),1, max))
  
  # Reconstruction following Yao, Müller & Wang (2005a) ---------------------
  N <- dim(Y0.sim)[2]
  grid <- seq(0, 1, length.out = N)
  Lu <- vector(mode = "list", length = T)
  Ly <- vector(mode = "list", length = T)
  for(t in 1:T) {
    O.t <- which(!is.na(Y0.sim[t,]))
    Lu[[t]] <- grid[O.t]
    Ly[[t]] <- Y0.sim[t,O.t]
  }
  
  fpca <- FPCA(Ly, Lu)
  X.hat.PACE <- fitted(fpca)[sim.incomplete,]
  
  MAE["PACE"] <- mean(apply(abs(X.hat.PACE -
                                  X0.sim[sim.incomplete,]),1, max))
  
  
  # Reconstruction following Kneip & Liebl ----------------------------------
  reconst_KL <- reconstructKneipLiebl(Ly, Lu,
                                      method = 'Error>0_AlignYES_CEscores',
                                      nRegGrid = N)
  X.hat.KL <- t(matrix(unlist(reconst_KL[['Y_reconst_list']]),
                       ncol = T))[sim.incomplete,]
  
  MAE["AYESCE"] <- mean(apply(abs(X.hat.KL - X0.sim[sim.incomplete,]),1, max))
  
  # Reconstruction following Elias, Jimenez & Wang (2023) -------------------

  X.hat.EJS <- reconstruct_EJS(Y0.sim, T.set = sim.incomplete)
  
  MAE["DEPTH"] <- mean(apply(abs(X.hat.EJS - X0.sim[sim.incomplete,]),1, max))
  
  # Reconstruction following Yao, Müller & Wang (2005b) ---------------------
  Lt <- list()
  Ly <- list()
  Lx <- list()
  
  for(t in sim.complete) {
    Lt[[t]] <- grid
    Lx[[t]] <- Y1.sim[t,]
    Ly[[t]] <- Y0.sim[t,]
  }
  X <- list(X = list(Ly = Lx, Lt = Lt))
  Y <- list(Ly = Ly, Lt = Lt)
  
  Lt <- list()
  Lx.test <- list()
  for(t in 1:length(sim.incomplete)) {
    Lt[[t]] <- grid
    Lx.test[[t]] <- Y1.sim[sim.incomplete[t],]
  }
  X.test <- list(X = list(Ly = Lx.test, Lt = Lt))
  
  reconst_FLM <- FLM1(Y, X, X.test)
  X.hat.FLM <- matrix(unlist(reconst_FLM[['yPred']]), nrow = length(sim.incomplete))
  
  MAE["FLM"] <- mean(apply(abs(X.hat.FLM - X0.sim[sim.incomplete,]),1, max))
  MAE
}

data_long <- melt(as.data.frame(res), variable.name = "method", value.name = "MAE")

plot_comparison <- ggplot(data_long, aes(x = method, y = MAE, fill = method)) +
  geom_boxplot() +
  labs(x = "", y = "MAE") +
  scale_y_continuous(limits = c(0, 6)) +
  theme_minimal()

pdf(paste0(getwd(), "/Results/Real data/", "weather_comparison.pdf"), height = 3,
    width = 6.1)
plot_comparison
dev.off()
