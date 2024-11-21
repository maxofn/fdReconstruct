#' Generate functional data
#'
#' @param T         Sample size (T = 100 by default).
#' @param N         Number of equispaced grid points (N = 51 by default).
#' @param r         Number of factors (r = 50 by default); must be an even number.
#' @param type.miss Type of missingness mechanism.
#' @param ev        Decay of eigenvalues; either "exp" (default) or polynomial else.
#' @param eps.sd    Standard deviation of measurement errors.
#' @param b         Coefficient vector for the covariate.
#' @param complete  If complete = TRUE (default), functions are completely
#'                  observable.
#'
#' @returns
#' A list containing the true data matrices 'X0', 'X1', the noisy data
#' matrices 'Y0', 'Y1', and the imperfectly observed data matrices
#' 'Y0.obs', 'Y1.obs' with NA's.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- GenObs(T = 4, N = 51)
#' par(mfrow = c(1,2))
#' matplot(t(data$Y0.obs), type = "l")
#' matplot(t(data$Y1.obs), type = "l")

GenObs <- function (T = 100, N = 51, r = 50, type.miss = "A", ev = "exp",
                    eps.sd = 0.1, b = c(1.1, 0.7, 0.5, 0.3), complete = FALSE) {
  if (r%%2 == 1) {
    stop("The parameter r must be even.")
  }
  grid <- seq(0, 1, length.out = N)
  mu <- sin(pi * grid)
  mu.mat <- matrix(rep(mu, T), byrow = TRUE, nrow = T, ncol = N)
  b.grid <- as.matrix(1:(r/2), nrow = r/2) %*% t(as.matrix(grid, nrow = N))
  sin.grid <- sqrt(2) * sin(2 * pi * b.grid)
  cos.grid <- sqrt(2) * cos(2 * pi * b.grid)
  r.vec <- 1:(r/2)

  if(ev == "exp") {
    xi1.sd <- sqrt(exp(- 2 * r.vec + 1))
    xi2.sd <- sqrt(exp(- 2 * r.vec    ))
  } else {
    # polynomially decaying eigenvalues
    xi1.sd <- sqrt((2 * r.vec - 1)^(-3)/2)
    xi2.sd <- sqrt((2 * r.vec    )^(-3)/2)
  }

  xi1 <- matrix(rep(xi1.sd, T), nrow = T, ncol = r/2, byrow = TRUE) *
    matrix(stats::rnorm(T * r/2), nrow = T, ncol = r/2)
  xi2 <- matrix(rep(xi2.sd, T), nrow = T, ncol = r/2, byrow = TRUE) *
    matrix(stats::rnorm(T * r/2), nrow = T, ncol = r/2)
  X0 <- xi1 %*% sin.grid + xi2 %*% cos.grid + mu.mat

  beta <- b[1] * t(t(sin.grid[1,])) %*%  sin.grid[1,] +
          b[2] * t(t(sin.grid[1,])) %*%  cos.grid[1,] +
          b[3] * t(t(cos.grid[1,])) %*%  sin.grid[1,] +
          b[4] * t(t(cos.grid[1,])) %*%  cos.grid[1,]

  X1 <- t(beta %*% t(X0) / N)

  Y0 <- X0 + stats::rnorm(length(X0), sd = eps.sd)
  Y1 <- X1 + stats::rnorm(length(X1), sd = eps.sd)

  Y0.obs <- Y0
  Y1.obs <- Y1

  if(complete == FALSE) {
    if (type.miss == "A") {
      for (t in 1:T) {
        D <- sample(floor(1 * N/2):floor(3 * N/4), size = 1)
        Y0.obs[t, D:N] <- NA
      }
    }
    if (type.miss == "B") {
      for (t in 1:T) {
        D <- sample(floor(1 * N/4):floor(3 * N/4), size = 1)
        Y0.obs[t, D:N] <- NA
      }
    }
    if (type.miss == "C") { # Not considered in paper.
      for (t in 1:T) {
        U <- sample(1:floor(N/2), size = 1)
        M <- c(floor(N/4):(floor(N/4) + U), (N - U):N)
        Y0.obs[t, M] <- NA
      }
    }
  }

  return(list(X0 = X0, X1 = X1, Y0 = Y0, Y1 = Y1, Y0.obs = Y0.obs,
              Y1.obs = Y1.obs))
}
