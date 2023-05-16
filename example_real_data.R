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
Y.obs <- matrix(data$temp_east, nrow = T, ncol = N, byrow = TRUE)
Y.aux <- matrix(data$temp_south, nrow = T, ncol = N, byrow = TRUE)

# Analysis ----------------------------------------------------------------

# How many days are completely observed?
is.complete <- which(rowSums(is.na(Y.obs))==0)
is.incomplete <- which(rowSums(is.na(Y.obs))>0)
length(is.complete) # 54

# Decay of eigenvalues
svd <- svd(cov(Y.obs[is.complete, ]/N))
ev.global <- svd$d

O <- which(!is.na(Y.obs[is.incomplete[2],]))
svd <- svd(cov(Y.obs[is.complete, O]/length(O)))
ev.local <- svd$d
plot(ev.local)
points(ev.global, pch = 2)

# Fraction of variance explained by the first eigenvalues
cumsum(ev.local[1:10])/(sum(ev.local))
cumsum(ev.global[1:10])/(sum(ev.global))

# Plot data ---------------------------------------------------------------

plot(data$temp_east, type = "l")
lines(data$temp_south, lty = 2)
matplot(t(Y.obs), type = "l")

# Reconstruct data --------------------------------------------------------

reconst <- ReconstFD(Y.obs = Y.obs, method = "Kraus", pred.band = TRUE,
                     p = 0.95, r.max = 15)
r.hat <- reconst$r

matplot(t(reconst$X.hat), type = "l")

# Plot reconstruction of partially observed curves ------------------------

time <- seq(
  from=as.POSIXct("2022-07-01 0:00", tz="UTC"),
  to=as.POSIXct("2022-07-01 23:30", tz="UTC"),
  by="30 min"
)

colours <- brewer.pal(8, name="Paired")

t.count <- 1
list <- list()
for(t.idx in is.incomplete) {
  X.hat.df<- data.frame(value = reconst$X.hat[t.idx,], Var1 = time)
  Y.obs.df<- data.frame(value = Y.obs[t.idx,], Var1 = time)
  Y.aux.df<- data.frame(value = Y.aux[t.idx,], Var1 = time)
  pred.df <- data.frame(u.df = as.vector(reconst$U.hat[t.idx,]),
                        l.df = as.vector(reconst$L.hat[t.idx,]),
                        Var1 = time)

  plot.new <- ggplot2::ggplot() +
    geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
                fill = "gray70", alpha = 0.3) +
    geom_point(data = Y.obs.df, aes(x = Var1, y = value),
               size = 3, alpha = 0.5, col = colours[t.count]) +
    geom_line(data = Y.aux.df, aes(x = Var1, y = value),
              col = "black", linewidth = 0.5, linetype = "dotted") +
    geom_line(data = X.hat.df, aes(x = Var1, y = value),
              col = colours[t.count], linewidth = 1) +
    scale_x_datetime(date_breaks = "360 min",
       date_labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
    xlab("time") +
    ylab("temperature (°C)") +
    theme_classic(base_size = 11) +
    theme(legend.position="none") +
    ggtitle(format(data$time[(t.idx - 1) * N + 1], "%Y-%m-%d"))
  list[[t.count]] <- plot.new
  t.count <- t.count + 1
}

pdf(paste0(getwd(), "/Res/Real data/", "reconstruction.pdf"), height = 11.69,
    width = 8.27)
do.call("ggarrange", c(list, nrow = 4, ncol = 2))
dev.off()
