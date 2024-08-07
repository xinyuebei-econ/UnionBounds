#' CI_SDRM Function
#'
#' This function calculates confidence intervals for the treatment effect under second differences relative magnitudes relaxation.
#'
#' @param betahat T by 1, regression coefficients, e.g. in Rambachan and Roth(2023) Equation 2.
#' @param betaSigma T by T, covariance matrix of betahat.
#' @param prePeriodIndices T_pre by 1, e.g. (1,2,3,4)
#' @param ell_post T_post by 1, t(ell)*tau_post is the parameter of interest
#' @param M 1 by 1, level of relaxation.
#' @param alpha 1 by 1, 1-alpha is the confidence level
#' @param alphac 1 by 1, a tuning parameter.
#' @param eta 1 by 1, a tuning parameter, see Bei(2024) Equation 40.
#' @param B The number of bootstrap draws used to caluclate c^t
#' @param Blarge The number of bootstrap draws used to caluclate c^t.
#' @param tol Tolerance level
#' @param tol_r Tolerance level, suggested value 10^-3.
#' @return  A list with ConfidenceInterval, delta_opt, c_t.
#' @export

CI_SDRM <- function(betahat, betaSigma, prePeriodIndices, ell_post, M = 1, alpha = 0.05, alphac = 0.8 * alpha, eta = 1e-3, B = 4000, Blarge = 40000, tol = 1e-3, tol_r = 1e-3) {
  # Calculate the relative magnitudes CI in Rambachan and Roth (2023, REStud) Section 2.4.1
  # betahat     T*1       regression coefficients, e.g. in RR23 (2),
  # deltaSigma  T*T       covariance matrix of betahat
  # ell         T_post*1  ell'*tau_post is the parameter of interest

  T_pre <- length(prePeriodIndices)
  T_post <- length(betahat) - T_pre
  ell <- c(rep(0, T_pre), ell_post)

  aux1 <- diag(T_pre)
  aux1 <- aux1[-nrow(aux1), ]
  aux2 <- -2 * diag(T_pre)
  aux2 <- aux2[-1, ]
  aux3 <- diag(T_pre + 1)
  aux3 <- aux3[, -ncol(aux3)]
  aux3 <- aux3[-c(1, 2), ]
  Lpre <- aux1 + aux2 + aux3
  Lpost <- rep(0, T_pre + T_post)
  H <- 1:T_post
  L <- numeric(T_post)
  for (i in 1:T_post) {
    L[i] <- i * (i + 1) / 2
  }
  Lpost[T_pre] <- sum(ell_post * H)
  Lpost <- Lpost + ell

  A <- rbind(cbind(Lpre, matrix(0, T_pre - 1, T_post)), Lpost)
  deltahat <- A %*% betahat
  deltaSigma <- A %*% betaSigma %*% t(A)

  Al <- cbind(c(ell_post %*% L * M) * rbind(diag(T_pre - 1), -diag(T_pre - 1)), matrix(rep(1, 2 * (T_pre - 1))))
  Au <- Al

  CI_result <- CI_ModifiedConditional(deltahat, deltaSigma, Al, Au, alpha, alphac, eta, B, Blarge, tol, tol_r)

  return(list(ConfidenceInterval = CI_result$ConfidenceInterval, delta_opt = CI_result$delta_opt, c_t = CI_result$c_t))
}
