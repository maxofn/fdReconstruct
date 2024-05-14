#' Reconstruct imperfectly observed functional data
#'
#' @param Y0.obs    Data matrix for target.
#' @param Y1.obs    Data matrix for covariate (NULL for univariate case).
#' @param T.set     Vector indicating functions to be reconstructed. By default,
#'                  all functions are reconstructed.
#' @param method    Method for determining O-specific number of factors. Must be
#'                  either 'CV' for 5-fold cross validation or 'manual' for
#'                  manual choice of number of factors.
#' @param r         Manual choice of O-specific number of factors (a vector).
#' @param r.max     Maximum number of factors.
#' @param w         Manual choice of weights (NULL uses integrated variances for
#'                  scaling).
#' @param pred.band Either TRUE or FALSE (default) indicating whether prediction
#'                  bands should be returned.
#' @param p         Coverage probability for prediction bands (p = 0.95 by
#'                  default.)
#'
#' @return          A list with the reconstructed data matrix 'X.hat', the
#'                  prediction bands 'U.hat' and 'L.hat' (only if
#'                  pred.band = TRUE) as well as vectors containing estimated
#'                  numbers of factors 'r.hat' and weights 'w.hat'.
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' data_comp <- GenObs(T = 50, complete = TRUE)
#' data_inco <- GenObs(T = 50, complete = FALSE)
#'
#' Y0.obs <- rbind(data_inco$Y0.obs, data_comp$Y0.obs)
#' Y1.obs <- rbind(data_inco$Y1.obs, data_comp$Y1.obs)
#'
#' reconst <- fdReconstruct(Y0.obs, Y1.obs, T.set = c(1), pred.band = TRUE)
#' plot(Y0.obs[1,], ylim = c(-0.5, 2.5))
#' lines(data_inco$X0[1,], lty = 2)
#' lines(reconst$X.hat[1,], col = "blue", lwd = 2)
#' lines(reconst$U.hat[1,], lty = 3)
#' lines(reconst$L.hat[1,], lty = 3)

fdReconstruct <- function (Y0.obs, Y1.obs = NULL, T.set = 1:dim(Y0.obs)[1],
                           method = "CV", r = NULL, r.max = 20, w = NULL,
                           pred.band = FALSE, p = 0.95) {
  T <- dim(Y0.obs)[1]
  N <- dim(Y0.obs)[2]

  if(is.null(Y1.obs)) {
    # Univariate case

    X.hat <- matrix(NA, nrow = length(T.set), ncol = N)
    if (pred.band) {
      U.hat <- matrix(NA, nrow = length(T.set), ncol = N)
      L.hat <- matrix(NA, nrow = length(T.set), ncol = N)
    }
    Yc <- Y0.obs[which(rowSums(is.na(Y0.obs)) == 0), ]
    Tc <- dim(Yc)[1]
    mu.hat <- colMeans(Yc)
    Yc <- Yc - matrix(rep(mu.hat, Tc), byrow = TRUE, nrow = Tc,
                      ncol = N)
    Y0.obs <- Y0.obs - matrix(rep(mu.hat, T), byrow = TRUE, nrow = T,
                              ncol = N)

    r.hat <- numeric(length(T.set))
    w.hat <- c(1)

    t.idx <- 1
    for (t in T.set) {
      O.t <- which(!is.na(Y0.obs[t, ]))
      M.t <- which(is.na(Y0.obs[t,]))

      if(is.null(r)) {
        r.hat[t.idx] <- RunCV(M.t, Yc, r.max = r.max)
      } else {
        r.hat[t.idx] <- r[t.idx]
      }

      svd <- svd(Yc[, O.t]/sqrt(Tc * N), nu = r.hat[t.idx], nv = r.hat[t.idx])
      F.hat.t <- 1/sqrt(N) * Y0.obs[t, O.t] %*% svd$v %*%
        diag(1/svd$d[1:r.hat[t.idx]], nrow = r.hat[t.idx], ncol = r.hat[t.idx])
      X.hat[t.idx, ] <- F.hat.t %*% t(svd$u) %*% Yc/sqrt(Tc) + mu.hat
      if (pred.band & length(O.t) < N) {
        band <- PredBand(Yc, O.set = O.t, p = p, r.max = r.max,
                         ro = r.hat[t.idx])
        U.hat[t.idx, -O.t] <- X.hat[t.idx, -O.t] + band$q.hat *
          band$Z.var.sqrt[-O.t]
        L.hat[t.idx, -O.t] <- X.hat[t.idx, -O.t] - band$q.hat *
          band$Z.var.sqrt[-O.t]
      }
      t.idx <- t.idx + 1
    }
  } else {
    # Multivariate case

    N0 <- dim(Y0.obs)[2]
    N1 <- dim(Y1.obs)[2]

    X.hat <- matrix(NA, nrow = length(T.set), ncol = N0)
    if (pred.band) {
      U.hat <- matrix(NA, nrow = length(T.set), ncol = N0)
      L.hat <- matrix(NA, nrow = length(T.set), ncol = N0)
    }

    Y0.c <- Y0.obs[which(rowSums(is.na(Y0.obs)) == 0), ]
    Y1.c <- Y1.obs[which(rowSums(is.na(Y0.obs)) == 0), ]

    Tc <- dim(Y0.c)[1]
    mu0.hat <- colMeans(Y0.c)
    mu1.hat <- colMeans(Y1.c)
    mu.hat  <- c(mu0.hat, mu1.hat)
    Y0.c <- Y0.c - matrix(rep(mu0.hat, Tc), byrow = TRUE, nrow = Tc, ncol = N0)
    Y1.c <- Y1.c - matrix(rep(mu1.hat, Tc), byrow = TRUE, nrow = Tc, ncol = N1)
    Y0.obs <- Y0.obs - matrix(rep(mu0.hat, T), byrow = TRUE, nrow = T,
                              ncol = N0)
    Y1.obs <- Y1.obs - matrix(rep(mu1.hat, T), byrow = TRUE, nrow = T,
                              ncol = N1)
    Y.obs <- cbind(Y0.obs, Y1.obs)

    r.hat <- numeric(length(T.set))

    w.hat <- c(NA, NA)
    if(is.null(w)) {
      w.hat[1] <- 1/mean(diag(stats::cov(Y0.c)))
      w.hat[2] <- 1/mean(diag(stats::cov(Y1.c)))
    } else {
      w.hat <- w
    }

    Yc <- cbind(sqrt(w.hat[1]) * Y0.c, sqrt(w.hat[2]) * Y1.c)

    t.idx <- 1

    for (t in T.set) {
      O.t <- which(!is.na(Y.obs[t, ]))
      M.t <- which(is.na(Y.obs[t,]))

      if(is.null(r)) {
        r.hat[t.idx] <- RunCV(M.t, Yc, r.max = r.max)
      } else {
        r.hat[t.idx] <- r[t.idx]
      }

      Y.obs.w <- cbind(sqrt(w.hat[1]) * Y0.obs[t,], sqrt(w.hat[2]) * Y1.obs[t,])

      svd <- svd(Yc[, O.t]/sqrt(Tc * N0), nu = r.hat[t.idx], nv = r.hat[t.idx])
      F.hat.t <- 1/sqrt(N0) * Y.obs.w[O.t] %*% svd$v %*%
        diag(1/svd$d[1:r.hat[t.idx]], nrow = r.hat[t.idx], ncol = r.hat[t.idx])
      X.hat[t.idx, ] <- F.hat.t %*% t(svd$u) %*% Yc[,1:N0]/sqrt(Tc)
      X.hat[t.idx, ] <- X.hat[t.idx,] / sqrt(w.hat[1]) + mu0.hat
      if (pred.band & length(O.t) < N0 + N1) {
        band <- PredBand(Y0.c, Y1.c, O.set = O.t, p = p, r.max = r.max,
                         ro = r.hat[t.idx])
        U.hat[t.idx, -O.t] <- X.hat[t.idx, -O.t] + band$q.hat *
          band$Z.var.sqrt[-O.t]
        L.hat[t.idx, -O.t] <- X.hat[t.idx, -O.t] - band$q.hat *
          band$Z.var.sqrt[-O.t]
      }
      t.idx <- t.idx + 1
    }
  }

  if (pred.band) {
    ls <- list(X.hat = X.hat, r = r.hat, w = w.hat, U.hat = U.hat,
               L.hat = L.hat)
  }
  else {
    ls <- list(X.hat = X.hat, r = r.hat, w = w.hat)
  }
  return(ls)
}

#' Construct simultaneous prediction bands
#'
#' @param Y0.c  Data matrix for target.
#' @param Y1.c  Data matrix for covariate (NULL for univariate case).
#' @param O.set Set on which functions are observable.
#' @param p     Coverage probability for prediction interval (p = 0.95 by
#'              default).
#' @param r.max Maximum number of factors (r.max = 20 by default).
#' @param ro    Number of factors on O.
#'
#' @returns
#' A list with the quantile 'q.hat' and shape 'Z.var.sqrt' for constructing
#' a simultaneous prediction band.

PredBand <- function(Y0.c, Y1.c = NULL, O.set, p = 0.95, r.max = 20, ro) {
  # Initialize variables
  Tc <- dim(Y0.c)[1]
  N  <- dim(Y0.c)[2]
  Z  <- matrix(NA, nrow = Tc, ncol = N)

  grid <- seq(0, 1, length.out = N)
  X0.hat <- matrix(NA, nrow = Tc, ncol = N)
  for(t in 1:Tc) {
    X0.hat[t,] <- stats::smooth.spline(grid, Y0.c[t,])$y
  }

  for(t in 1:Tc) {
    Y0.obs <- Y0.c
    Y0.obs[t, -O.set] <- NA
    reconst <- fdReconstruct(Y0.obs, Y1.c, T.set = c(t), method = "manual",
                             pred.band = FALSE, r = ro, r.max = r.max)
    Z[t,] <- X0.hat[t,] - as.vector(reconst$X.hat)
  }

  Z.var.sqrt <- sqrt(diag(stats::cov(Z)))
  Z.var.sqrt[O.set] <- NA

  vec <- numeric(Tc)
  for(t in 1:Tc) {
    vec[t] = max(abs(Z[t,-O.set])/Z.var.sqrt[-O.set])
  }
  q.hat <- stats::quantile(vec, probs = p)
  return(list(q.hat = q.hat, Z.var.sqrt = Z.var.sqrt))
}
