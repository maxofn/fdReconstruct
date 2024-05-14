directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)
library(fdReconstruct)
library(fdapace)
library(ReconstPoFD)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(latex2exp)
require(reshape2)

BLUE <- brewer.pal(3, name="Paired")[2]

# Exponentially decaying eigenvalues --------------------------------------
set.seed(1)

T <- 100
data_test  <- GenObs(T = 1, N = 51, r = 50, type.miss = "A", ev = "exp",
                     eps.sd = 0.05, complete = FALSE)
data_train <- GenObs(T = T, N = 51, r = 50, type.miss = "A", ev = "exp",
                     eps.sd = 0.05, complete = TRUE)

Y0.obs <- rbind(data_test$Y0.obs, data_train$Y0.obs)
Y1.obs <- rbind(data_test$Y1.obs, data_train$Y1.obs)

reconst_mult <- fdReconstruct(Y0.obs, Y1.obs, T.set = c(1), pred.band = TRUE)

grid <- seq(0, 1, length.out = 51)
x.hat.df <- data.frame(Var1 = grid, value = as.vector(reconst_mult$X.hat))
x.df <- data.frame(Var1 = grid, value = data_test$X0[1,])
y.df <- data.frame(Var1 = grid, value = Y0.obs[1,])
pred.df <- data.frame(Var1 = grid, u.df = as.vector(reconst_mult$U.hat),
                      l.df = as.vector(reconst_mult$L.hat))

plot.exp <- ggplot2::ggplot() +
  geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
              fill = "gray70", alpha = 0.3) +
  geom_point(data = y.df, aes(x = Var1, y = value),
             size = 3, alpha = 0.5, col = BLUE) +
  geom_line(data = x.df, aes(x = Var1, y = value),
            linewidth = 0.5, linetype = 3) +
  geom_line(data = x.hat.df, aes(x = Var1, y = value),
            linewidth = 1, col = BLUE) +
  xlab("") +
  ylab("") +
  theme_minimal(base_size = 11) +
  theme(legend.position="bottom")
plot.exp

# Polynomially decaying eigenvalues --------------------------------------
set.seed(1)

T <- 100
data_test  <- GenObs(T = 1, N = 51, r = 50, type.miss = "A", ev = "poly",
                     eps.sd = 0.05, complete = FALSE)
data_train <- GenObs(T = T, N = 51, r = 50, type.miss = "A", ev = "poly",
                     eps.sd = 0.05, complete = TRUE)


Y0.obs <- rbind(data_test$Y0.obs, data_train$Y0.obs)
Y1.obs <- rbind(data_test$Y1.obs, data_train$Y1.obs)

reconst_mult <- fdReconstruct(Y0.obs, Y1.obs, T.set = c(1), pred.band = TRUE)

grid <- seq(0, 1, length.out = 51)
x.hat.df <- data.frame(Var1 = grid, value = as.vector(reconst_mult$X.hat))
x.df <- data.frame(Var1 = grid, value = data_test$X0[1,])
y.df <- data.frame(Var1 = grid, value = Y0.obs[1,])
pred.df <- data.frame(Var1 = grid, u.df = as.vector(reconst_mult$U.hat),
                      l.df = as.vector(reconst_mult$L.hat))

plot.poly <- ggplot2::ggplot() +
  geom_ribbon(data = pred.df, aes(x = Var1, ymin = l.df, ymax = u.df),
              fill = "gray70", alpha = 0.3) +
  geom_point(data = y.df, aes(x = Var1, y = value),
             size = 3, alpha = 0.5, col = BLUE) +
  geom_line(data = x.df, aes(x = Var1, y = value),
            linewidth = 0.5, linetype = 3) +
  geom_line(data = x.hat.df, aes(x = Var1, y = value),
            linewidth = 1, col = BLUE) +
  xlab("") +
  ylab("") +
  theme_minimal(base_size = 11) +
  theme(legend.position="bottom")
plot.poly

pdf(paste0(getwd(), "/Results/Reconstruction/", "simulation_examples.pdf"),
    height = 3, width = 8.27)
do.call("ggarrange", c(list(plot.exp, plot.poly), nrow = 1, ncol = 2))
dev.off()

