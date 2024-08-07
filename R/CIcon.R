#' CIcon Function
#'
#' This function calculates the conditional confidence interval.
#'
#' @param deltahat A numeric vector.
#' @param deltaSigma A matrix.
#' @param Al A matrix.
#' @param Au A matrix.
#' @param alphac A numeric value.
#' @param CI_p A numeric vector.
#' @param tol A numeric value.
#' @param tol_r A numeric value.
#' @return A numeric vector representing the confidence interval.
#' @export

CIcon <- function(deltahat, deltaSigma, Al, Au, alphac, CI_p, tol, tol_r) {
  kk <- nrow(Al)
  lambdaSigma <- rbind(Al, Au) %*% deltaSigma %*% t(rbind(Al, Au))
  lambda_sigma <- sqrt(diag(lambdaSigma))
  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  sigma_l <- lambda_sigma[1:kk]
  sigma_u <- lambda_sigma[(kk + 1):(2 * kk)]

  lb <- Al %*% deltahat - sigma_l * stats::qnorm(1 - 10^-4)
  ub <- Au %*% deltahat + sigma_u * stats::qnorm(1 - 10^-4)
  lb <- min(lb)
  ub <- max(ub)

  # CI lower bound
  rej <- 1
  theta <- lb
  while (rej == 1 && theta <= CI_p[1]) {
    Tcl <- min((Al %*% deltahat - theta) / sigma_l)
    Tcu <- min((theta - Au %*% deltahat) / sigma_u)
    Tc <- max(Tcl, Tcu)

    bounds <- CIcon_TNbounds(theta, t(deltahat), Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, tol_r)
    t_val <- (1 - alphac) * stats::pnorm(bounds$th_2) + alphac * stats::pnorm(bounds$th_1)
    rej <- (stats::pnorm(Tc) > t_val)
    theta <- theta + tol
  }
  CI_lower <- min(theta, CI_p[1])

  # CI upper bound
  rej <- 1
  theta <- ub
  while (rej == 1 && theta >= CI_p[2]) {
    Tcl <- min((Al %*% deltahat - theta) / sigma_l)
    Tcu <- min((theta - Au %*% deltahat) / sigma_u)
    Tc <- max(Tcl, Tcu)

    bounds <- CIcon_TNbounds(theta, t(deltahat), Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, tol_r)
    t_val <- (1 - alphac) * stats::pnorm(bounds$th_2) + alphac * stats::pnorm(bounds$th_1)
    rej <- (stats::pnorm(Tc) > t_val)
    theta <- theta - tol
  }
  CI_upper <- max(theta, CI_p[2])

  CIcon <- c(CI_lower, CI_upper)
  return(CIcon)
}
