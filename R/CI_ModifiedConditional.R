#' CI_ModifiedConditional Function
#'
#' Calculates the confidence interval for \eqn{\theta \in \left[\min_{b \in \mathcal{B}} \lambda_{\ell,b}, \max_{b \in \mathcal{B}} \lambda_{u,b}\right]} using the modified conditional inference proposed by Bei(2024) "Inference on Union Bounds".
#' Assume that the estimators \eqn{\hat{\lambda_{\ell}} = A_{\ell}\hat{\delta}} and \eqn{\hat{\lambda_{u}} = A_{u}\hat{\delta}}. \eqn{\hat{\delta}} is asymptotically normal with covariance matrix \eqn{deltaSigma}.
#'
#' @param deltahat p by 1 vector, see Bei(2024) Assumption 1.
#' @param deltaSigma p by p matrix, see Bei(2024) Assumption 1.
#' @param Al \eqn{|\mathcal{B}|} by p matrix, see Bei(2024) Assumption 1.
#' @param Au \eqn{|\mathcal{B}|} by p matrix, see Bei(2024) Assumption 1.
#' @param alpha The nominal rejection rate.
#' @param alphac A tuning parameter, suggested value is alphac = 0.8*alpha, see Bei(2024) Section 3.2.
#' @param eta A tuning parameter, suggested value is 0.001, see Bei(2024) eq(40).
#' @param B The number of bootstrap draws used to caluclate c^t, suggested value 4000.
#' @param Blarge The number of bootstrap draws used to caluclate c^t, suggested value B*10.
#' @param tol Tolerance level to do the grid search, suggested value 10^-3.
#' @param tol_r Tolerance level to determine |rho|=1, suggested value 10^-3.
#' @return A list with ConfidenceInterval, c_MC_l, c_MC_u, c_t, delta_opt.
#' @examples
#' # Example:  min{delta_1, delta_2} <= theta <= max{delta_3, delta_2}
#' set.seed(0)
#' deltaSigma <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
#' deltahat <- MASS::mvrnorm(n = 1, c(0,0, 0), deltaSigma)
#' deltahat <- matrix(deltahat, nrow = 3, ncol = 1)
#' Al <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
#' Au <- matrix(c(0, 0, 1, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
#' # Call the CI_ModifiedConditional function
#' CI <- CI_ModifiedConditional(deltahat, deltaSigma, Al, Au, alpha = 0.05)
#' print(CI$ConfidenceInterval) # This is the confidence interval
#' print(CI$c_MC_l) # The lower bound of CI is min{lambda_{l,b} - c_MC_l*sigma_l}
#' print(CI$c_MC_u) # The upper bound of CI is min{lambda_{u,b} - c_MC_u*sigma_l}
#' print(CI$c_t) # c_t is the lower bound of c_MC_l and c_MC_u, see Bei(2024) Section 3.3
#' print(CI$delta_opt) # the maximizer in Bei(2024) eq(26)
#' @export

CI_ModifiedConditional<- function(deltahat, deltaSigma, Al, Au, alpha = 0.05, alphac = alpha*0.8, eta = 10^-3, B = 4000, Blarge = B*10, tol = 10^-3, tol_r = 10^-3) {

  index <- c()

  # Loop over columns of Al
  for (i in 1:ncol(Al)) {
    # Check if all elements in the column are the same in both Al and Au
    if (length(unique(c(Al[, i],Au[, i]))) == 1)  {
      index <- c(index, i)
    }
  }

  # Calculate \hat{c}_{\eta}
  set.seed(0)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, length(deltahat)), Sigma = deltaSigma)
  deviation <- apply( sweep(abs(deltastar_demean_large), 2, sqrt(diag(deltaSigma)), FUN = "/"), 1, max)
  c_bd <- stats::quantile(deviation, 1 - eta,type = 1)

  sigma_delta <- sqrt(diag(deltaSigma))
  lambdaSigma <- rbind(Al, Au) %*% deltaSigma %*% t(rbind(Al, Au))
  sigma_lambda <- sqrt(diag(lambdaSigma))
  kk <- nrow(Al)
  sigma_l <- sigma_lambda[1:kk]
  sigma_u <- sigma_lambda[(kk + 1):(2 * kk)]

  corr_all <- solve(diag(sigma_lambda)) %*% lambdaSigma %*% solve(diag(sigma_lambda))
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]

  # Construct \hat\Delta
  lb <- deltahat - sigma_delta * c_bd
  ub <- deltahat + sigma_delta * c_bd
  delta1 <- deltahat
  lb[index] <- 0
  ub[index] <- 0
  delta1[index] <- 0

  cl <- 0
  cu <- stats::qnorm(1 - alpha / 2)
  c <- (cl + cu) / 2
  delta_fea <- list()
  c_fea <- numeric()

  obj_large <- function(delta, c_check) {
    (alpha - CIproj_p(c_check, delta, alphac, deltastar_demean_large, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, eta, tol_r)) * 100
  }

  k <- 1

  while (cu - cl > tol) {
    set.seed(k)
    k <- k + 1

    # Redraw deltastar_demean to reduce the sampling bias
    deltastar_demean <- MASS::mvrnorm(B, mu = rep(0, length(deltahat)), Sigma = deltaSigma)

    obj <- function(delta) {
      (alpha - CIproj_p(c, delta, alphac, deltastar_demean, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, eta, tol_r)) * 100
    }
    p1 <- obj_large(delta1, c)
    if (p1 >= 0) {

      opts <- list( "algorithm" = "NLOPT_GN_ORIG_DIRECT","xtol_rel" = 1.0e-4)
      result <- nloptr::nloptr( x0 = delta1, eval_f = obj, opts = opts, lb = lb, ub = ub)
      delta1 <- result$solution

      if (result$objective >0) {
        opts <- list( "algorithm" = "NLOPT_GN_CRS2_LM","xtol_rel" = 1.0e-4)
        result2 <- nloptr::nloptr( x0 = delta1, eval_f = obj, opts = opts, lb = lb, ub = ub)
        if (result2$objective < result$objective){
          delta1 <- result2$solution
        }
      }


      # Use a relatively large bootstrap sample B*10 to reduce the sampling bias along direction delta1
      p1 <- obj_large(delta1, c)
      if (p1 >= 0) {
        cu <- c
      } else {
        cl <- c
        delta_fea <- append(delta_fea, list(delta1))
        c_fea <- c(c, c_fea)
      }
    } else {
      cl <- c
      delta_fea <- append(delta_fea, list(delta1))
      c_fea <- c(c, c_fea)
    }
    c <- (cl + cu) / 2
  }

  # Double check using the final direction delta1
  cl <- c
  cu <- stats::qnorm(1 - alpha / 2)
  while (cu - cl > tol) {
    c <- (cl + cu) / 2
    p <- obj_large(delta1, c)
    if (p >= 0) {
      cu <- c
    } else {
      cl <- c
    }
  }

  lambdahat_l <- Al %*% deltahat
  lambdahat_u <- Au %*% deltahat

  CI_p <- c(min(lambdahat_l - c * sigma_l), max(lambdahat_u + c * sigma_u))
  # Calculate the conditional confidence interval
  CI_h <- CIcon(deltahat, deltaSigma, Al, Au, alphac, CI_p, tol, tol_r)

  c_MC_l =  min((lambdahat_l - CI_h[1]) / sigma_l)
  c_MC_u =  min((CI_h[2] - lambdahat_u) / sigma_u)
  return(list(ConfidenceInterval = CI_h, c_MC_l = c_MC_l, c_MC_u = c_MC_u, c_t = c, delta_opt = delta1))
}
