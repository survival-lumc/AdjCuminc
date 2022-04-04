#' @title Simulate Balanced Competing Risk Data under the Cause-Specific Hazards Model
#'
#' @description The purpose of this function is to simulate a competing risks data set from the cause-specific hazard model, with the following set of predictors: a strata (cohort), a discrete covariate (mutation), and a continuous covariate (age). The observation times are simulated according to the latent competing risk approach, where the first event can censor the event of interest.
#'
#' @param N The number of observations to be generated.
#' @param k The number of competing outcomes to be generated.
#' @param base_haz a baseline hazard for each competing outcome
#' @param beta_coef A list where the number of elements is equal to the number of competing outcomes. Aach listed element is a vector of covariate coefficients associated with a competing outcome. The first element in the vector is the effect size of the cohort, the second element the effect of the discrete covariate, and the third element the effect of the continuous covariate.
#'
#' @return A `data.frame` containing simulated observation times, a status indicator, cohort assignments, a discrete covariate (mutation), and a continuous covariate (age)
#' @export
#' @importFrom stats quantile rbinom rexp rnorm runif
#'
#' @examples
#' simCSH(N = 1000,
#' k = 2,
#' base_haz = c(1, 1),
#' beta_coef = list(c(-0.5, 2, 0.5), c(0.3, 1, -0.4)))
simCSH <- function(N, k, base_haz, beta_coef) {

  # covariate simulation
  age <- rnorm(N)
  mutant_type <- sample(c(0, 1, 2), N, c(1/3, 1/3, 1/3), replace = T)
  cohort <- rbinom(N, 1, 0.5)
  X <- as.matrix(data.frame(cohort, mutant_type, age))

  # simulate event times per event k
  event_times <- matrix(NA, nrow = N, ncol = k)
  for(i in 1:k){
    haz_fun <- base_haz[i] * exp(X %*% beta_coef[[i]])
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
  sim_data <- data.frame(obs_time, status = factor(status),
                         cohort = factor(ifelse(cohort == 1, "treat", "control")),
                         mutation = factor(paste("stage", mutant_type)),
                         age)
  return(sim_data)
}
