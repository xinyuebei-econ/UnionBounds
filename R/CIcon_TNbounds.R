#' CIcon_TNbounds Function
#'
#' This function calculates the bounds of the truncated Normal.
#'
#' @param theta A numeric value. H0.
#' @param deltahat A B*p matrix.
#' @param Al A matrix.
#' @param Au A matrix.
#' @param sigma_l A numeric vector.
#' @param sigma_u A numeric vector.
#' @param corr_m A matrix.
#' @param corr_l A matrix.
#' @param corr_u A matrix.
#' @param tol_r A numeric value.
#' @return A list with th_1 and th_2.
#' @export

CIcon_TNbounds <- function(theta, deltahat, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, tol_r) {
  B <- nrow(deltahat)
  lambdahat_l <- deltahat %*% t(Al)
  lambdahat_u <- deltahat %*% t(Au)

  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, FUN = "/")
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, FUN = "/")
  Tl <- apply(Tlb, 1, min)
  Tu <- apply(Tub, 1, min)

  bl <- apply(Tlb, 1, which.min)
  bu <- apply(Tub, 1, which.min)

  corr_m_blb <- corr_m[bl, ]
  corr_m_bub <- t(corr_m[, bu])

  taux_l1 <- matrix(10^10, B, nrow(Au))
  taux_l1a <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  taux_l1[1 + corr_m_blb > tol_r] <- taux_l1a[1 + corr_m_blb > tol_r]

  th_1_1 <- apply(taux_l1, 1, min) * (Tl >= Tu)
  th_1_1[th_1_1 == 10^10] <- -10^10
  th_1_1 <- th_1_1 * (Tl >= Tu)

  taux_l2 <- matrix(10^10, B, nrow(Al))
  taux_l2a <- (1 - corr_l[bl, ])^(-1) * (Tlb - corr_l[bl, ] * Tl)
  taux_l2[1 > corr_l[bl, ] + tol_r] <- taux_l2a[1 > corr_l[bl, ] + tol_r]

  th_2_1 <- apply(taux_l2, 1, min) * (Tl >= Tu)

  taux_u1 <- matrix(10^10, B, nrow(Al))
  taux_u1a <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  taux_u1[1 + corr_m_bub > tol_r] <- taux_u1a[1 + corr_m_bub > tol_r]

  th_1_2 <- apply(taux_u1, 1, min) * (Tl < Tu)
  th_1_2[th_1_2 == 10^10] <- -10^10
  th_1_2 <- th_1_2 * (Tl < Tu)

  taux_u2 <- matrix(10^10, B, nrow(Au))
  taux_u2a <- (1 - corr_u[bu, ])^(-1) * (Tub - corr_u[bu, ] * Tu)
  taux_u2[1 > corr_u[bu, ] + tol_r] <- taux_u2a[1 > corr_u[bu, ] + tol_r]

  th_2_2 <- apply(taux_u2, 1, min) * (Tl < Tu)

  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2

  return(list(th_1 = th_1, th_2 = th_2))
}
