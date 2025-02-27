#' The twostep algorithm
#'
#' This function implement the OSMAC method for the input covariate matrix @param X, response vector Y, first step sample size r1, the second step sample size r2, and the method to use.
#' It returns a list with three elements: par, the weighted MLE; se, the standard errors; msg, the fitting message; iter, the number of itterations used; method, the method used.

#' @param X the input covariate matrix
#' @param Y the input response vector
#' @param r1 the first step sample size
#' @param r2 the second step sample size
#' @param method the method to use
#' @keywords getMLE twostep
#' @export
#' @examples
#' library(OSMAC)
#' dat <- adult.train
#' X <- as.matrix(dat[,c(1,3,5,12:13)])
#' X <- t(t(X) / apply(X, 2, sd))
#' X <- cbind(1, X)
#' Y <- as.numeric(dat[,15]) - 1
#' set.seed(0)
#' twostep(X, Y, 200, 800, "mmse")

twostep <- function(X, Y, r1, r2, method=c("mvc", "mmse", "uni")) {
  call <- match.call()
  method <- match.arg(method)
  n <- length(Y)
  if (method == "uni") {
    idx.simp <- sample(1:n, r1, T)
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    if (fit.simp$message == "Successful convergence") {
      p.simp  <- 1 - 1 / (1 + exp(c(x.simp %*% beta.simp)))
      w.simp <- p.simp * (1 - p.simp)
      W.simp <- solve(t(x.simp) %*% (x.simp * (w.simp * n))) * r1 * n
      Vc.simp <- t(x.simp) %*% (x.simp * (y.simp-p.simp)^2 * n^2) / r1^2 / n^2
      V.simp <- W.simp %*% Vc.simp %*% W.simp
      se.simp <- sqrt(diag(V.simp))
    }
    else {
      se.simp = NA
    }
    result <- list(par=beta.simp, se=se.simp, msg=msg, method="uni")
  }
  else {
    n1 <- sum(Y)
    n0 <- n - n1
    PI.prop <- rep(1/(2*n0), n)
    PI.prop[Y==1] <- 1/(2*n1)
    idx.prop <- sample(1:n, r1, T, PI.prop)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- n
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if (anyNA(beta.prop)) {
      result <- list(opt=NA, msg="first stage not converge")
    }
    else if (method == "mmse") {
      P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
      p.prop <- P.prop[idx.prop]
      w.prop <- p.prop * (1 - p.prop)
      W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
      PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
      PI.mMSE <- PI.mMSE / sum(PI.mMSE)
      idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
      
      ru <- length(y.mMSE)
      beta.mMSE <- fit.mMSE$par
      p.mMSE  <- 1 - 1 / (1 + exp(c(x.mMSE %*% beta.mMSE)))
      w.mMSE <- p.mMSE * (1 - p.mMSE)
      W.mMSE <- solve(t(x.mMSE) %*% (x.mMSE * (w.mMSE * pinv.mMSE))) * ru * n
      Vc.mMSE <- t(x.mMSE) %*% (x.mMSE * (y.mMSE-p.mMSE)^2 * pinv.mMSE^2) / ru^2 / n^2
      V.mMSE <- W.mMSE %*% Vc.mMSE %*% W.mMSE
      
      se <- sqrt(diag(V.mMSE))
      msg <- c(fit.prop$message, fit.mMSE$message)
      result <- list(par=beta.mMSE, se=se, msg=msg, method="mmse")
    }
    else if (method == "mvc") {
      P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
      PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
      PI.mVc <- PI.mVc / sum(PI.mVc)
      idx.mVc <- sample(1:n, r2, T, PI.mVc)
      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
      
      ru <- length(y.mVc)
      beta.mVc <- fit.mVc$par
      p.mVc  <- 1 - 1 / (1 + exp(c(x.mVc %*% beta.mVc)))
      w.mVc <- p.mVc * (1 - p.mVc)
      W.mVc <- solve(t(x.mVc) %*% (x.mVc * (w.mVc * pinv.mVc))) * ru * n
      Vc.mVc <- t(x.mVc) %*% (x.mVc * (y.mVc-p.mVc)^2 * pinv.mVc^2) / ru^2 / n^2
      V.mVc <- W.mVc %*% Vc.mVc %*% W.mVc
      ## ## LCC
      ## PI.optU <- abs(Y - P.prop)
      ## PI.optU <- PI.optU / sum(PI.optU)
      ## idx.optU <- sample(1:n, r2, T, PI.optU)
      ## x.optU <- X[c(idx.optU),]
      ## y.optU <- Y[c(idx.optU)]
      ## pinv.optU <- 1
      ## fit.optU <- getMLE(x=x.optU, y=y.optU, w=pinv.optU)
      
      se <- sqrt(diag(V.mVc))
      msg <- c(fit.prop$message, fit.mVc$message)
      result <- list(par=beta.mVc, se=se, msg=msg, method="mvc")
    }
  }
  return(result)
}