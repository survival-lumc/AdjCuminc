#' @title Confounded sampling from simulated competing risk data under the Cause-Specific Hazards Model
#'
#' @description These functions were created to simulate a biased sample from a cause-specific hazards model. The formulation of the outcome model, from which the observations are sampled. Three variants have been introduced to simulate different scenarios of model misspecification. Adjustment of the regression coefficients in beta_coef determines the strength and direction of the association between the covariates and treatment allocation, where beta_coef = 0 results in strict random selection on a covariate, beta_coef > 0 results in positive selection for higher covariate values, and beta_coef > 0 results in positive selection for lower covariate values. The association of these regression coefficients on the control treatment is reversed.
#'
#' The following scenario's can be simulated with this function:
#'
#' Scenario 1: introduces a logistic association between treatment allocation and the discrete and continuous covariates.
#'
#' Scenario 2: incorporates a non-logistic association by sampling individuals under the treatment group under scenario 1, and the control group from scenario 0. This can be achieved by setting control_ref = TRUE.
#'
#' Scenario 3: introduces scenario 1, but incorporates a non-linear effect by adjusting the baseline hazard conditional under the continuous covariate. This can be induced by using two different values for base_haz.
#'
#' @param N the number of observations to be generated.
#' @param k the number of competing outcomes to be generated.
#' @param theta_coef a vector where each element represents the coefficient of treatment on the hazard of each event.
#' @param beta_coef a list where the number of elements is equal to the number of competing outcomes. Aach listed element is a vector of covariate coefficients associated with a competing outcome. The first element in the vector is the effect size of the cohort, the second element the effect of the discrete covariate, and the third element the effect of the continuous covariate.
#' @param sample_coef a vector of length 3 that introduce non-random assignment to a treatment group. If any of these elements are are not 0, scenario 1 is induced, where the first two elements affect probability of assignment conditional on the discrete covariate, and the third element affects the assignment conditional on the continuous covariate.
#' @param control_ref a logical indicating whether control observations should be resampled under scenario 0. If `TRUE`, Scenario 2 is induced where the covariates for `control`-labeled observations are are resampled under scenario 0 (prior to event time simulation).
#' @param base_haz a baseline hazard for each competing outcome.
#'
#' @return A `data.frame` containing simulated observation times, a status indicator, treatment assignments, a discrete covariate (x1), and a continuous covariate (x2)
#' @export
#' @importFrom stats quantile rbinom rexp rnorm runif
#'
#' @examples
#' confCSH(N = 4000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   sample_coef = c(1, -1, 1),
#'   control_ref = FALSE,
#'   base_haz = c(0, 0))

confCSH <- function(N, k, theta_coef, beta_coef, sample_coef, control_ref = FALSE, base_haz) {

  # Scenario 0. Simulate covariates under random allocation
  x1 <- t(rmultinom(N, 1, c(1/3, 1/3, 1/3)))[, -1]
  x2 <- rnorm(N)

  # Scenario 1. Perform assignment conditional on covariates
  prob <- plogis(x1 %*% sample_coef[1:2] + x2 * sample_coef[3])
  Z <- rbinom(N, 1, prob)
  X <- cbind(x1, x2)

  # Scenario 2. Replace scenario 1 controls by scenario 0 controls
  if(control_ref == TRUE){
    control_ind <- which(Z == 0)
    n_control <- length(control_ind)
    X[control_ind, ] <- cbind(t(rmultinom(n_control, 1, c(1/3, 1/3, 1/3)))[, -1],
                              rnorm(n_control))
  }

  # Scenario 3. Induce non-linear baseline difference conditional on x1
  base_hazard <- ifelse(X[, 3] > 1, base_haz[1], base_haz[2])

  # simulate event times per event k
  event_times <- matrix(NA, nrow = N, ncol = k)
  for(i in 1:k){
    haz_fun <- exp(base_hazard) * exp(Z * theta_coef[i] + X %*% beta_coef[[i]])
    event_times[, i] <- rexp(length(haz_fun), rate = haz_fun)
  }

  # simulate competition for first occurrence between k events
  obs_time <- apply(event_times, 1, min)
  status <- apply(event_times, 1, function(x) which(x == min(x)))

  # Add censoring
  censor_time <- runif(N, quantile(obs_time, 0.20), quantile(obs_time, 0.95))
  status <- ifelse(censor_time < obs_time, 0, status)
  obs_time <- pmin(censor_time, obs_time)

  # construct output
  x1_label <- apply(X[, 1:2], 1,
                    function(x) ifelse(sum(x) == 0, 0, which(x == 1)))
  sim_data <- data.frame(obs_time, status = factor(status),
                         treatment = factor(ifelse(Z == 1, "treat", "control")),
                         x1 = factor(paste("stage", x1_label)),
                         x2 = X[, 3])
  return(sim_data)
}
