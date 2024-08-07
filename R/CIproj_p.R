#' CIproj_p Function
#'
#' This function calculates the rejection rate given true value delta and lower truncation c
#'
#' @param c A numeric value. The lower truncation.
#' @param delta A numeric vector. The true value of delta.
#' @param alphac A numeric value.
#' @param deltastar_demean A matrix.
#' @param Al A B*p matrix, see Assumption 1.
#' @param Au A B*p matrix, see Assumption 1.
#' @param sigma_l A numeric vector.
#' @param sigma_u A numeric vector.
#' @param corr_m A matrix defined in Bei(2024) Lemma 1.
#' @param corr_l A matrix defined in Bei(2024) Lemma 1.
#' @param corr_u A matrix defined in Bei(2024) Lemma 1.
#' @param eta A numeric value.
#' @param tol_r A numeric value.
#' @return A numeric value.
#' @export

CIproj_p <- function(c, delta, alphac, deltastar_demean, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, eta, tol_r) {
  # Calculate lambda_l and lambda_u
  lambda_l <- Al %*% delta
  lambda_u <- Au %*% delta
  lb <- min(lambda_l)
  ub <- max(lambda_u)
  mb <- (lb + ub) / 2

  deltastar <- matrix(rep(t(delta), each = nrow(deltastar_demean)), nrow = nrow(deltastar_demean), byrow = F)
  deltastar = deltastar + deltastar_demean
  lambdastar_l <- deltastar %*% t(Al)
  lambdastar_u <- deltastar %*% t(Au)

  # Calculate Tstar_l, Tstar_m, and Tstar_u
  Tstar_l <- apply(cbind(apply(sweep(lambdastar_l - lb, 2, sigma_l, FUN = "/"),1,min), apply(sweep(lb - lambdastar_u, 2, sigma_u, FUN = "/"),1,min)),1,max)
  Tstar_m <- apply(cbind(apply(sweep(lambdastar_l - mb, 2, sigma_l, FUN = "/"),1,min), apply(sweep(mb - lambdastar_u, 2, sigma_u, FUN = "/"),1,min)),1,max)
  Tstar_u <- apply(cbind(apply(sweep(lambdastar_l - ub, 2, sigma_l, FUN = "/"),1,min), apply(sweep(ub - lambdastar_u, 2, sigma_u, FUN = "/"),1,min)),1,max)


  ind_proj_l <- (Tstar_l <= c)
  ind_proj_m <- (Tstar_m <= c)
  ind_proj_u <- (Tstar_u <= c)

  # Function to calculate condition proj
  condition_proj <- function(select, lambdastar_l, lambdastar_u, b, sigma_l, sigma_u) {
    if (sum(select) == 0) {
      ind_c <- rep(0, nrow(deltastar))
    } else {
      sub_deltastar <- deltastar[select == 1, , drop = FALSE]
      sub_lambdastar_l <- lambdastar_l[select == 1, , drop = FALSE]
      sub_lambdastar_u <- lambdastar_u[select == 1, , drop = FALSE]

      Tc1 <- apply(sweep(sub_lambdastar_l - b, 2, sigma_l, FUN = "/"), 1, min)
      Tc2 <- apply(sweep(b - sub_lambdastar_u, 2, sigma_u, FUN = "/"), 1, min)

      Tc <- apply(cbind(Tc1, Tc2),1,max)
      bounds <- CIcon_TNbounds(b, sub_deltastar, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, tol_r)
      phi <- (stats::pnorm(Tc) - stats::pnorm(bounds$th_1)) / (stats::pnorm(bounds$th_2) - stats::pnorm(bounds$th_1))

      sub_ind_c <- (phi <= 1 - alphac)
      ind_c <- rep(0, nrow(deltastar))
      ind_c[select == 1] <- sub_ind_c
    }
    return(ind_c)
  }

  ind_c_l <- condition_proj(1 - ind_proj_l, lambdastar_l, lambdastar_u, lb, sigma_l, sigma_u)
  ind_c_m <- condition_proj(1 - ind_proj_m, lambdastar_l, lambdastar_u, mb, sigma_l, sigma_u)
  ind_c_u <- condition_proj(1 - ind_proj_u, lambdastar_l, lambdastar_u, ub, sigma_l, sigma_u)

  ind_l <- pmax(ind_c_l, ind_proj_l)
  ind_m <- pmax(ind_c_m, ind_proj_m)
  ind_u <- pmax(ind_c_u, ind_proj_u)

  p1 <- mean(1 - ind_l * pmax(ind_m, ind_u))
  p2 <- mean(1 - ind_u * pmax(ind_m, ind_l))
  p <- max(p1, p2) + eta

  return(p)
}
