#' Run 5-fold cross validation for estimating the unknown rank.
#'
#' @param M.set   Set on which observations are missing. If M.set = NULL (i.e.
#'                all observations are available) the rank is chosen such that
#'                the fraction of variance explained exceeds 0.99.
#' @param Y.obs   Data matrix.
#' @param r.max   Maximum number of factors (r.max = 20 by default).
#'
#' @return Estimated number of factors.
#' @export
#'
#' @examples
#' set.seed(123)
#' data_comp <- GenObs(T = 50, N = 51, r = 4, complete = TRUE)
#' data_inco <- GenObs(T = 50, N = 51, r = 4, complete = FALSE)
#' Y0.obs <- rbind(data_inco$Y0.obs, data_comp$Y0.obs)
#' Y1.obs <- rbind(data_inco$Y1.obs, data_comp$Y1.obs)
#' RunCV(M.set = 1:10, cbind(Y0.obs, Y1.obs))

RunCV <- function(M.set = NULL, Y.obs, r.max = 20) {

  N <- dim(Y.obs)[2]

  completely.obs <- which(rowSums(is.na(Y.obs)) == 0)
  Yc <- Y.obs[completely.obs, ]
  Tc <- dim(Yc)[1]

  r.max <- min(r.max, min(floor(Tc * 4/5), N - length(M.set)) - 1)

  if(is.null(M.set)) {
    ev <- eigen(stats::cov(Yc))$values
    r.hat <- min(min(which(cumsum(ev)/sum(ev) > 0.99)), r.max)
  } else {
    sse <- numeric(r.max)
    list <- split(1:Tc, factor(sort(rep(1:5, length.out = Tc))))

    for(k in 1:5) {
      # Run 5-fold cross validation
      test.obs <- list[[k]]
      Y.train <- Yc[-test.obs,]
      Y.test  <- Yc[ test.obs,]
      Y.dummy <- Yc[ test.obs,]
      Y.dummy[, M.set] <- NA

      T.train <- dim(Y.train)[1]
      T.test  <- dim(Y.test )[1]

      for (r in 1:r.max) {
        reconst <- fdReconstruct(rbind(Y.dummy, Y.train), T.set = 1:T.test,
                                 method = 'manual', r = rep(r, T.test),
                                 r.max = r.max)
        X.hat <- reconst$X.hat
        sse[r] <- sse[r] + sum((Y.test[, M.set] - X.hat[, M.set])^2)
      }
    }
    r.hat   <- which.min(sse)
  }

  return(r.hat)
}
