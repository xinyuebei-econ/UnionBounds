#' CI_RM
#'
#' Constructs robust confidence intervals for the difference-in-differences estimand under relative magnitudes relaxation, see Rambachan and Roth (2023, REStud) Section 2.4.1.
#'
#' @param betahat T by 1, difference-in-differences coefficients, e.g. in Rambachan and Roth(2023, REStud) Equation 2.
#' @param betaSigma T by T, covariance matrix of betahat.
#' @param prePeriodIndices T_pre by 1, e.g. (1,2,3,4)
#' @param ell_post T_post by 1, ell%*%tau_post is the parameter of interest, see Rambachan and Roth(2023, REStud) Section 2.3.
#' @param M 1 by 1, level of relaxation, see Rambachan and Roth (2023, REStud) Section 2.4.1.
#' @param alpha 1 by 1, 1-alpha is the confidence level
#' @param alphac 1 by 1, a tuning parameter, see Bei(2024) Section 3.2.
#' @param eta 1 by 1, a tuning parameter, see Bei(2024) Equation 40.
#' @param B The number of bootstrap draws used to caluclate c^t
#' @param Blarge The number of bootstrap draws used to caluclate c^t.
#' @param tol Tolerance level to calculate the confidence interval.
#' @param tol_r Tolerance level to check whether |rho| = 1
#' @return  A list with ConfidenceInterval, c_MC_l, c_MC_u, c_t, delta_opt
#' @examples
#' # Load example data
#' # the sample data is from Dustmann et al (2022), see Bei(2024) Section 5
#' data(example_data, package = "UnionBounds")
#' betahat <- example_data$betahat
#' betaSigma <- example_data$betaSigma
#' prePeriodIndices <- example_data$prePeriodIndices
#'
#' # Call the CI_RM function with the example data
#' ell_post <- c(1, 0)  # the parameter of interest is the treatment effect at time 1
#' CI <- CI_RM(betahat, betaSigma, prePeriodIndices, ell_post, M = 1, alpha = 0.05)
#' print(CI$ConfidenceInterval) # This is the RM confidence interval
#'
#' @export

CI_RM <- function(betahat, betaSigma, prePeriodIndices, ell_post, M = 1, alpha = 0.05, alphac = 0.8 * alpha, eta =  1e-3, B = 4000, Blarge = 40000, tol = 1e-3, tol_r = 1e-3) {
  # Calculate the relative magnitudes CI in Rambachan and Roth (2023, REStud) Section 2.4.1

  T_pre <- length(prePeriodIndices)
  T_post <- length(betahat) - T_pre
  ell <- c(rep(0, T_pre), ell_post)

  aux1 <- -diag(T_pre)
  aux2 <- diag(T_pre + 1)[-1, -ncol(diag(T_pre + 1))]
  L_pre <- aux1 + aux2

  A <- rbind(cbind(diag(T_pre), matrix(0, T_pre, T_post)), ell)
  deltahat <- A %*% betahat
  deltaSigma <- A %*% betaSigma %*% t(A)

  gamma <- sum(abs(ell_post))
  Al <- rbind(cbind(gamma * M * rbind(diag(T_pre), -diag(T_pre)), rep(1, 2 * T_pre))) %*% rbind(cbind(L_pre, matrix(0, T_pre, 1)), matrix(c(rep(0, T_pre), 1), 1, T_pre + 1))
  Au <- Al

  CI_result <- CI_ModifiedConditional(deltahat, deltaSigma, Al, Au, alpha, alphac, eta, B, Blarge, tol, tol_r)

  return(list(ConfidenceInterval = CI_result$ConfidenceInterval))
}
