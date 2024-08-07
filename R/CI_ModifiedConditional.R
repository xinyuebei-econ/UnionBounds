#' CI_ModifiedConditional Function
#'
#' This function calculates the CI base on Bei(2024) "Inference on Union Bounds"
#'
#' @param deltahat A  p*1 vector, see Assumption 1.
#' @param deltaSigma A p*p matrix, see Assumption 1.
#' @param Al A B*p matrix, see Assumption 1.
#' @param Au A B*p matrix, see Assumption 1.
#' @param alpha The nominal rejection rate.
#' @param alphac A tuning parameter, suggested value is alphac = 0.8*alpha.
#' @param eta A tuning parameter, suggested value is 0.001, see Bei(2024) eq(40).
#' @param B The number of bootstrap draws used to caluclate c^t, suggested value 4000.
#' @param Blarge The number of bootstrap draws used to caluclate c^t, suggested value B*10.
#' @param tol Tolerance level, suggested value 10^-3.
#' @param tol_r Tolerance level, suggested value 10^-3.
#' @return A list with ConfidenceInterval, delta_opt, c_t.
#' @examples
#' # Load example data
#' data(example_data, package = "UnionBounds")
#' betahat <- example_data$betahat
#' betaSigma <- example_data$betaSigma
#' prePeriodIndices <- example_data$prePeriodIndices
#'
#' # Call the CI_SDRM function with the example data
#' ell_post <- c(1, 0)  # the parameter of interest is tau_1, i.e. the treatment effect at time 1
#' CI <- CI_ModifiedConditional(deltahat, deltaSigma, Al, Au, alpha = 0.05)
#' print(CI$ConfidenceInterval) # This is the confidence interval
#' @export

CI_ModifiedConditional<- function(deltahat, deltaSigma, Al, Au, alpha = 0.05, alphac = alpha*0.8, eta = 10^-3, B = 4000, Blarge = B*10, tol = 10^-3, tol_r = 10^-3) {

  index <- c()

  # Loop over columns of Al
  for (i in 1:ncol(Al)) {
    # Check if all elements in the column are the same in both Al and Au
    if (length(unique(Al[, i])) == 1 && length(unique(Au[, i])) == 1) {
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

  return(list(ConfidenceInterval = CI_h, delta_opt = delta1, c_t = c))
}
