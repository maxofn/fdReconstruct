#' Reconstruct imperfectly observed functional data
#'
#' @param Y.obs Data matrix.
#' @param T.set Vector indicating functions to be reconstructed. By default,
#' all functions are reconstructed.
#' @param method Method for determining O-specific number of factors. Must be
#' either 'BaiNg', 'Kraus' (default), 'Onatski' or 'rGiven' (for manual choice).
#' @param pred.band Either TRUE or FALSE (default) indicating whether
#' prediction bands should be returned.
#' @param p Coverage probability for prediction intervals (p = 0.95 by
#' default).
#' @param r Manual choice of O-specific number of factors (only if method =
#' "rGiven").
#' @param r.max Maximum number of factors.
#'
#' @returns
#' A list with the reconstructed data matrix 'X.hat', the prediction bands
#' 'U.hat' and 'L.hat' (only if pred.band = TRUE) as well as a vector
#' containing estimated numbers of factors 'r.hat'.
#'
#' @export
#'
#' @references
#' Bai, Jushan and Serena Ng (2002). “Determining the number of factors in
#' approximate factor models”. In: Econometrica 70.1, pp. 191–221.
#'
#' Kraus, David (2015). “Components and completion of partially observed
#' functional data”. In: J. R. Stat. Soc. Ser. B. Stat. Methodol. 77.4,
#' pp. 777–801.
#'
#' Onatski, Alexei (2010). “Determining the number of factors from empirical
#' distribution of eigenvalues”. In: The Review of Economics and Statistics
#' 92.4, pp. 1004–1016.
#'
#' @examples
#' set.seed(1)
#' data_com <- GenObs(T = 100, N = 51, r = 10, complete = TRUE, type.miss = "A")
#' data_inc <- GenObs(T = 1, N = 51, r = 10, complete = FALSE, type.miss = "A")
#' data <- list(X = rbind(data_com$X, data_inc$X),
#'              Y = rbind(data_com$Y, data_inc$Y),
#'              Y.obs = rbind(data_com$Y.obs, data_inc$Y.obs))
#'
#' # Reconstruct functions
#' reconst <- ReconstFD(data$Y.obs, pred.band = TRUE, method = "Onatski")
#'
#' # Plot reconstructions
#' t.idx <- which(rowSums(is.na(data$Y.obs)) > 0)[1]
#' lim <- range(data$X[t.idx], data$Y.obs[t.idx,], reconst$U.hat[t.idx,],
#'              reconst$L.hat[t.idx,], na.rm = TRUE)
#' plot(data$X[t.idx,], type = "l", lty = 3, ylim = lim)
#' points(data$Y.obs[t.idx,])
#' lines(reconst$X.hat[t.idx,], lwd = 2)
#' lines(reconst$U.hat[t.idx,], lty = 2)
#' lines(reconst$L.hat[t.idx,], lty = 2)

ReconstFD <- function(Y.obs, T.set = 1:dim(Y.obs)[1], method = 'Kraus',
                      pred.band = FALSE, p = 0.95, r = NULL, r.max = 10) {
  # Choose method for tuning r
  method <- switch(method,
                   'BaiNg'   = 1,
                   'Kraus'   = 2,
                   'Onatski' = 3,
                   'rGiven'  = 4,
  )
  if(method == 4 & length(r) != length(T.set)) {
    stop("Length of r does not match number of functions to be reconstructed.")
  }

  # Initialize variables
  T <- dim(Y.obs)[1]
  N <- dim(Y.obs)[2]
  X.hat <- matrix(NA, nrow = length(T.set), ncol = N)

  if(pred.band) {
    U.hat <- matrix(NA, nrow = length(T.set), ncol = N)
    L.hat <- matrix(NA, nrow = length(T.set), ncol = N)
  }

  # Complete observations
  Yc <- Y.obs[which(rowSums(is.na(Y.obs)) == 0),]
  Tc <- dim(Yc)[1]

  # Center data
  mu.hat <- colMeans(Yc)
  Yc <- Yc - matrix(rep(mu.hat, Tc), byrow = TRUE, nrow = Tc, ncol = N)
  Y.obs <- Y.obs - matrix(rep(mu.hat, T), byrow = TRUE, nrow = T, ncol = N)

  # Impute missing values
  r.hat <- numeric(length(T.set))
  t.idx <- 1
  for(t in T.set) {
    O.t <- which(!is.na(Y.obs[t,]))
    No <- length(O.t)
    if(method == 1) {# BaiNg
      r.hat[t.idx] = EstNumBaiNg(Yc[,O.t], r.max = r.max, center = FALSE)
    }
    if(method == 2) {# Kraus
      r.hat[t.idx] = EstNumKraus(O.t, Yc, r.max = r.max)
    }
    if(method == 3) {# Onatski
      r.hat[t.idx] = EstNumOnatski(Yc[,O.t], r.max = r.max, center = FALSE)
    }
    if(method == 4) {# rGiven
      r.hat[t.idx] = r[t.idx]
    }

    # Reconstruct function
    svd <- svd(Yc[,O.t]/sqrt(Tc * No), nu = r.hat[t.idx], nv = r.hat[t.idx])
    F.hat.t <- 1/sqrt(No) * Y.obs[t,O.t] %*% svd$v %*%
      diag(1/svd$d[1:r.hat[t.idx]], nrow = r.hat[t.idx], ncol = r.hat[t.idx])
    X.hat[t.idx,] <- F.hat.t %*% t(svd$u) %*% Yc / sqrt(Tc) + mu.hat

    # Compute prediction bands
    if(pred.band & No < N) {
      band <- PredBand(Yc, O.set = O.t, p = p, r.max = r.max,
                       ro = r.hat[t.idx])
      U.hat[t.idx,-O.t] <- X.hat[t.idx,-O.t] + band$q.hat*band$Z.var.sqrt[-O.t]
      L.hat[t.idx,-O.t] <- X.hat[t.idx,-O.t] - band$q.hat*band$Z.var.sqrt[-O.t]
    }
    t.idx <- t.idx + 1
  }

  if(pred.band) {
    ls <- list(X.hat = X.hat, r = r.hat, U.hat = U.hat, L.hat = L.hat)
  } else {
    ls <- list(X.hat = X.hat, r = r.hat)
  }
  return(ls)
}

#' Construct simultaneous prediction bands
#'
#' @param Yc Data matrix.
#' @param O.set Set on which functions are observable.
#' @param p Coverage probability for prediction interval (p = 0.95 by
#' default).
#' @param r.max Maximum number of factors (r.max = 10 by default).
#' @param ro Number of factors on O.
#'
#' @returns
#' A list with the quantile 'q.hat' and shape 'Z.var.sqrt' for constructing
#' a simultaneous prediction band.

PredBand <- function(Yc, O.set, p = 0.95, r.max = 10, ro) {
  # Initialize variables
  Tc <- dim(Yc)[1]
  N <- dim(Yc)[2]
  Z <- matrix(NA, nrow = Tc, ncol = N)

  r <- max(EstNumOnatski(Yc, r.max = r.max, center = FALSE), ro)
  reconst <- ReconstFD(Yc, method = "rGiven", pred.band = FALSE,
                       r = rep(r, Tc), r.max = r.max)
  X0.hat <- reconst$X.hat

  for(t in 1:Tc) {
    Y.obs <- Yc
    Y.obs[t, -O.set] <- NA
    reconst <- ReconstFD(Y.obs, T.set = c(t), method = "rGiven",
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
