directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)
library(FDFM)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(latex2exp)
require(reshape2)

BLUE <- brewer.pal(3, name="Paired")[2]
set.seed(1)

N = 51
grid <- seq(0, 1, length.out = N)
r = 10

data.test <- GenObs(T = 1, N = N, r = r, complete = FALSE, type.miss = 'B')

data   <- GenObs(T = 200, N = N, r = r, complete = TRUE , type.miss = 'B')
data$X <- rbind(data$X, data.test$X)
data$Y <- rbind(data$Y, data.test$Y)
data$Y.obs <- rbind(data$Y.obs, data.test$Y.obs)

reconst <- ReconstFD(data$Y.obs, T.set = c(201), method = "rGiven", r = c(7),
                     pred.band = TRUE, p = 0.95, r.max = 15)
reconst$r

x.hat.df <- data.frame(Var1 = grid, value = as.vector(reconst$X.hat))
x.df <- data.frame(Var1 = grid, value = data$X[201,])
y.df <- data.frame(Var1 = grid, value = data$Y.obs[201,])
pred.df <- data.frame(Var1 = grid, u.df = as.vector(reconst$U.hat),
                      l.df = as.vector(reconst$L.hat))

plot.example <- ggplot2::ggplot() +
  geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
              fill = "gray70", alpha = 0.3) +
  geom_point(data = y.df, aes(x = Var1, y = value),
             size = 5, alpha = 0.5, col = BLUE) +
  geom_line(data = x.df, aes(x = Var1, y = value),
            linewidth = 1.5, linetype = 3) +
  geom_line(data = x.hat.df, aes(x = Var1, y = value),
            linewidth = 1.5, col = BLUE) +
  xlab("") +
  ylab("") +
  theme_minimal(base_size = 24) +
  theme(legend.position="bottom") +
  ggtitle("DGPB")
plot.example

pdf(paste0(getwd(), "/Res/Bands/", "ExampleB.pdf"), height = 6,
    width = 8.5)
plot.example
dev.off()

data.test <- GenObs(T = 1, N = N, r = r, complete = FALSE, type.miss = 'A')

data   <- GenObs(T = 200, N = N, r = r, complete = TRUE , type.miss = 'A')
data$X <- rbind(data$X, data.test$X)
data$Y <- rbind(data$Y, data.test$Y)
data$Y.obs <- rbind(data$Y.obs, data.test$Y.obs)

reconst <- ReconstFD(data$Y.obs, T.set = c(201), method = "rGiven", r = c(7),
                     pred.band = TRUE, p = 0.95, r.max = 15)
reconst$r

x.hat.df <- data.frame(Var1 = grid, value = as.vector(reconst$X.hat))
x.df <- data.frame(Var1 = grid, value = data$X[201,])
y.df <- data.frame(Var1 = grid, value = data$Y.obs[201,])
pred.df <- data.frame(Var1 = grid, u.df = as.vector(reconst$U.hat),
                      l.df = as.vector(reconst$L.hat))

plot.example <- ggplot2::ggplot() +
  geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
              fill = "gray70", alpha = 0.3) +
  geom_point(data = y.df, aes(x = Var1, y = value),
             size = 5, alpha = 0.5, col = BLUE) +
  geom_line(data = x.df, aes(x = Var1, y = value),
            linewidth = 1.5, linetype = 3) +
  geom_line(data = x.hat.df, aes(x = Var1, y = value),
            linewidth = 1.5, col = BLUE) +
  xlab("") +
  ylab("") +
  theme_minimal(base_size = 24) +
  theme(legend.position="bottom") +
  ggtitle("DGPA")
plot.example

pdf(paste0(getwd(), "/Res/Bands/", "ExampleA.pdf"), height = 6,
    width = 8.5)
plot.example
dev.off()

