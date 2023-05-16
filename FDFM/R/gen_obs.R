#' Generate functional data
#'
#' @param T Sample size (T = 100 by default).
#' @param N Number of equispaced grid points (N = 51 by default).
#' @param r Number of global factors (r = 10 by default); must be an even number.
#' @param complete Either TRUE (default) indicating that functions are
#' observed on the whole domain, or FALSE otherwise.
#' @param type.miss Type of missingness mechanism; see paper.
#'
#' @returns
#' A list containing the true data matrix 'X', the noisy data matrix 'Y' and
#' the imperfectly observed data matrix 'Y.obs' with NA's.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- GenObs(T = 100, N = 51, r = 10, complete = FALSE, type.miss = "A")
#' matplot(t(data$Y.obs), type = "l")

GenObs <- function(T = 100, N = 51, r = 10, complete = TRUE, type.miss = 'A') {
  # Check if r is even
  if(r %% 2 == 1) {
    stop("The parameter r must be even.")
  }

  # Construct basis
  grid <- seq(0, 1, length.out = N)
  mu <- sin(pi * grid)
  mu.mat <- matrix(rep(mu, T), byrow = TRUE, nrow = T, ncol = N)
  b.grid <- as.matrix(1:(r/2), nrow = r/2) %*% t(as.matrix(grid, nrow = N))
  sin.grid <- sqrt(2) * sin(2 * pi * b.grid)
  cos.grid <- sqrt(2) * cos(2 * pi * b.grid)
  r.vec <- 1:(r/2)
  xi1.sd <- sqrt((2 * r.vec - 1)^(-3))
  xi2.sd <- sqrt((2 * r.vec    )^(-3))

  # Simulate curves
  xi1 <- matrix(rep(xi1.sd, T), nrow = T, ncol = r/2, byrow = TRUE) *
    matrix(stats::rnorm(T*r/2), nrow = T, ncol = r/2)
  xi2 <- matrix(rep(xi2.sd, T), nrow = T, ncol = r/2, byrow = TRUE) *
    matrix(stats::rnorm(T*r/2), nrow = T, ncol = r/2)

  X <- xi1 %*% sin.grid + xi2 %*% cos.grid + mu.mat
  Y <- X + stats::rnorm(T * N, sd = sqrt(0.01))

  # Miss observations
  Y.obs <- Y

  if(!complete) {
    if(type.miss == 'A') {
      for(t in 1:T) {
        d <- stats::runif(1, min = 0.5, max = 2/3)
        M.t <- which(grid > d)
        Y.obs[t, M.t] <- NA
      }
    }
    if(type.miss == 'B') {
      for(t in 1:T) {
        a <- stats::runif(1, min = 0, max = 1/2)
        b <- a + 1/2
        M.t <- which((grid > a) & (grid < b))
        Y.obs[t, M.t] <- NA
      }
    }
  }

  return(list(X = X, Y = Y, Y.obs = Y.obs))
}
