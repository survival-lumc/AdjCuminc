#' @title Simulating Balanced Competing Risk Data under the Cause-Specific Hazards Model
#'
#' @description The purpose of this function is to simulate a competing risks data set from the cause-specific hazard model, with the following set of predictors: a strata (treatment), a discrete covariate (x1), and a continuous covariate (x2). The observation times are simulated according to the latent competing risk approach, where the first event can censor the event of interest.
#'
#' @param N the number of observations to be generated.
#' @param k the number of competing outcomes to be generated.
#' @param theta_coef a vector where each element represents the coefficient of treatment on the hazard of each event.
#' @param base_haz a baseline hazard for each competing outcome
#' @param beta_coef a list where the number of elements is equal to the number of competing outcomes. Aach listed element is a vector of covariate coefficients associated with a competing outcome. The first element in the vector is the effect size of the cohort, the second element the effect of the discrete covariate, and the third element the effect of the continuous covariate.
#'
#' @return A `data.frame` containing simulated observation times, a status indicator, treatment assignments, a discrete covariate (x1), and a continuous covariate (x2)
#' @export
#' @importFrom stats quantile rbinom rexp rnorm runif
#'
#' @examples
#' simCSH(N = 10000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))

simCSH <- function(N, k, theta_coef, beta_coef, base_haz) {

  # covariate simulation
  x1 <- t(rmultinom(N, 1, c(1/3, 1/3, 1/3)))[, -1]
  x2 <- rnorm(N)
  X <- cbind(x1, x2)
  Z <- rbinom(N, 1, 0.5)

  # Induce non-linear baseline difference conditional on x2
  base_hazard <- ifelse(X[, 3] > 1, base_haz[1], base_haz[2])

  # simulate event times per event k
  event_times <- matrix(NA, nrow = N, ncol = k)
  for(i in 1:k){
    haz_fun <- exp(base_hazard) * exp(Z * theta_coef[i] + X %*% beta_coef[[i]])
    event_times[, i] <- rexp(length(haz_fun), rate = haz_fun)
  }

  # simulate competition for first occurence between k events
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
