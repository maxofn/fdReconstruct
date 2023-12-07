directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)
library(FDFM)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(latex2exp)
require(reshape2)

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
length(is.incomplete) # 76

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
            col = "gray80", linewidth = 0.3) +
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
            col = "gray80", linewidth = 0.3) +
  scale_x_datetime(date_breaks = "360 min",
                   date_labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  xlab("time") +
  ylab("temperature (°C)") +
  theme_classic(base_size = 11) +
  theme(legend.position="none") +
  ggtitle("West")
plot.sample1

pdf(paste0(getwd(), "/Res/Real data/", "sample.pdf"), height = 3,
    width = 8.27)
do.call("ggarrange", c(list(plot.sample0, plot.sample1), nrow = 1, ncol = 2))
dev.off()


# Reconstruct data --------------------------------------------------------

reconst_mult <- ReconstFD(Y0.obs = Y0.obs, Y1.obs = Y1.obs, T.set = is.incomplete,
                           pred.band = TRUE, p = 0.95, r.max = 30, w = c(10,1))
reconst_uni  <- ReconstFD(Y0.obs = Y0.obs, T.set = is.incomplete,
                             pred.band = TRUE, p = 0.95, r.max = 15)

reconst_mult$r # [1] 17 13 18 16 17 11 18 12 18 15

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

pdf(paste0(getwd(), "/Res/Real data/", "weather_reconstructions.pdf"),
    height = 11.69, width = 8.27)
do.call("ggarrange", c(list, nrow = 5, ncol = 2))
dev.off()

pdf(paste0(getwd(), "/Res/Real data/", "weather_examples.pdf"), height = 3,
    width = 8.27)
do.call("ggarrange", c(list(list[[6]], list[[8]]), nrow = 1, ncol = 2))
dev.off()
