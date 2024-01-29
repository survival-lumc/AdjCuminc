#' @title Hypothetical true competing risk curve computation
#
#' @description This function was developed to calculate the hypothetical true curve for a series of time points, based on the cause-specific hazard model. The hypothetical true curve can be controlled by a list of baseline hazard, and treatment covariate effect coefficients per event. Non-linearity conditional on the continuous covariate can be introduced by choosing two different values for `base_haz`.
#'
#' @param time a vector of time points.
#' @param k the number of events.
#' @param theta_coef  vector where each element represents the coefficient of treatment on the hazard of each event.
#' @param beta_coef a list where each listed element is a vector of covariate coefficients associated with event k. The first two elements introduce a discrete covariate effect on the outcome, and the third element a continuous covariate effect on the outcome.
#' @param base_haz a vector of length 2 with a baseline hazard conditional on the continuous covariate. If these two values differ, scenario 3 is induced, resulting in a non-linear relationship between covariates and outcome.
#' @return A `data.frame` with the hypothetical true curve per treatment condition
#' @export
#' @import riskRegression
#' @import survival
#' @importFrom prodlim jackknife prodlim
#' @importFrom stats model.frame terms
#'
#' @examples
#' df_true <- trueCSH(time = times, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1,  -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))

trueCSH <- function(time, k, theta_coef, beta_coef, base_haz){

  # Create combinations discrete covariates
  X_design <- as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1), 0))
  X_design <- X_design[rowSums(X_design[, 2:3]) <= 1, ]
  X_coefs <- cbind(unname(theta_coef), do.call("rbind", beta_coef))
  X_coefs <- unname(split(X_coefs, 1:nrow(X_coefs)))
  X_beta <- lapply(X_coefs,  function(b) X_design %*% b)

  # Construct curve per discrete combination per event k
  curves_Ik <- vector("list", k)
  for(i in 1:k){
    curves_ij <- matrix(NA, nrow = length(time), ncol = 6)
    for(j in 1:6){
      lambda_k <- function(x) {
        exp(ifelse(x > 1, base_haz[1], base_haz[2])) *
          exp(X_beta[[i]][j] + x * X_coefs[[i]][[4]])
      }

      sum_k <- paste(paste0(" exp(ifelse(x > 1, base_haz[1], base_haz[2])) *
                            exp(X_beta[[", 1:k,
                            "]][j] + x * X_coefs[[", 1:k,
                            "]][4])"), collapse = "+")

      lambda_K <- function(x) {
        eval(parse(text = sum_k))
      }

      for(l in 1:length(time)) {
        surv_fun <- function(x) {
          lambda_k(x) / lambda_K(x) * (1 - exp(-time[l] * lambda_K(x))) * dnorm(x)
        }
        curves_ij[l, j] <- integrate(surv_fun, -10, 10)$value
      }
    }
    curves_Ik[[i]] <- c(rowMeans(curves_ij[, c(1, 3, 5)]),
                        rowMeans(curves_ij[, c(2, 4, 6)]))
  }

  # Construct output
  pstate <- do.call("cbind", curves_Ik)
  pstate <- cbind(1 - rowSums(pstate), pstate)
  colnames(pstate) <- c("(s0)", 1:k)
  output <- cbind(data.frame(time = rep(time, times = 2),
                             treatment = rep(c("control", "treat"), each = length(time))),
                  pstate)
  return(output)
}
