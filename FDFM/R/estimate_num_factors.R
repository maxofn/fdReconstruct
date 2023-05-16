#' Estimate number of factors according to Bai and Ng (2002).
#'
#' @param Y Data matrix.
#' @param r.max Maximum number of factors (r.max = 10 by default).
#' @param center Either TRUE (default) indicating that data matrix must be
#' centered first, or FALSE.
#'
#' @returns
#' The estimated number of factors (between 1 and r.max).
#'
#' @export
#'
#' @references
#' Bai, Jushan and Serena Ng (2002). “Determining the number of factors in
#' approximate factor models”. In: Econometrica 70.1, pp. 191–221.
#'
#' @examples

EstNumBaiNg <- function(Y, r.max = 10, center = TRUE) {
  T <- dim(Y)[1]
  N <- dim(Y)[2]

  if(center) {
    mu.hat <- colMeans(Y)
    Y <- Y - matrix(rep(mu.hat, T), nrow = T, ncol = N, byrow = TRUE)
  }

  IC <- numeric(r.max)
  F.hat <- sqrt(T) * svd(Y, nu = r.max, nv = r.max)$u

  for(r in 1:r.max) {
    gTN <- (N + T)/(N*T) * log((N*T)/(N + T))
    lognorm <- log(mean((Y - 1/T * F.hat[,1:r] %*% t(F.hat[,1:r]) %*% Y)^2))
    IC[r] <- lognorm + r * gTN
  }
  r <- which.min(IC)
  return(r)
}

#' Estimate number of factors similar to Kraus (2015).
#'
#' @param O.set Set on which functions are observable.
#' @param Y.obs Imperfectly observed data matrix.
#' @param r.max Maximum number of factors (r.max = 10 by default).
#'
#' @return
#' The estimated number of factors (between 1 and r.max).
#'
#' @export
#'
#' @references
#' Kraus, David (2015). “Components and completion of partially observed
#' functional data”. In: J. R. Stat. Soc. Ser. B. Stat. Methodol. 77.4,
#' pp. 777–801.
#'
#' @examples

EstNumKraus <- function(O.set, Y.obs, r.max = 10) {
  # Return r.max if there are no missing values
  N <- dim(Y.obs)[2]
  if(length(O.set) == N) {
    return(r.max)
  }

  # Discard incompletely observed curves
  completely.obs <- which(rowSums(is.na(Y.obs)) == 0)
  Yc <- Y.obs[completely.obs, ]
  Tc <- dim(Yc)[1]

  # Run GCV
  RSS <- numeric(r.max)

  for(t in 1:Tc) {
    Yc.dummy <- Yc
    Yc.dummy[t, -O.set] <- NA

    for(r in 1:r.max) {
      if(r >= length(O.set)) {
        RSS[r:r.max] <- Inf
        break
      } else {
        reconst <- ReconstFD(Yc.dummy, T.set = c(t), method = 'rGiven',
                             r = c(r), r.max = r.max)
        X.hat.t <- reconst$X.hat
        RSS[r] <- RSS[r] +
          sum((Yc[t,-O.set] - X.hat.t[1,-O.set])^2)/(N - length(O.set))
      }
    }
  }

  GCV <- RSS/(1 - (1:r.max)/Tc)^2
  r <- which.min(GCV)
  return(r)
}

#' Estimate number of factors according to Onatski (2010).
#'
#' @param Y Data matrix.
#' @param r.max Maximum number of factors (r.max = 10 by default).
#' @param center Either TRUE (default) indicating that data matrix must be
#' centered first, or FALSE.
#' @param niter Maximum number of iterations (niter = 10 by default).
#'
#' @returns
#' The estimated number of factors (between 1 and r.max).
#'
#' @export
#'
#' @references
#' Onatski, Alexei (2010). “Determining the number of factors from empirical
#' distribution of eigenvalues”. In: The Review of Economics and Statistics
#' 92.4, pp. 1004–1016.
#'
#' @examples

EstNumOnatski <- function(Y, r.max = 10, center = TRUE, niter = 10) {
  T <- dim(Y)[1]
  N <- dim(Y)[2]

  if(center) {
    mu.hat <- colMeans(Y)
    Y <- Y - matrix(rep(mu.hat, T), nrow = T, ncol = N, byrow = TRUE)
  }

  lambda <- svd(Y, nu = 1, nv = 1)$d^2/T
  j <- r.max + 1

  for(iter in 1:niter) {
    y <- lambda[j:(j+4)]
    x <- ((j-1):(j+3))^(2/3)
    beta.hat <- stats::lm(y ~ x)$coefficients[2]
    delta <- 2 * abs(beta.hat)
    gap <- lambda[1:r.max] - lambda[2:(r.max+1)]
    idx <- which(gap > delta)
    if(length(idx) == 0) {
      r.hat <- 0
    } else {
      r.hat <- max(idx)
    }
    j <- r.hat + 1
  }
  return(max(1,r.hat))
}

