#' @title Generate Inverse Probability Weights
#'
#' @description A simple function to generates inverse probability weights via standard logistic regression
#'
#' @param formula A `formula` object in the form of `y ~ x`, where `x` represents the collection of covariates which should be balanced between the levels of `y`.
#' @param data An optional `data.frame` used to interpret the variables described by `formula` argument.
#'
#' @return A `vector` containing the inverse probability weights for each observation
#' @export
#' @importFrom stats glm predict
#'
#' @examples
#' df <- simCSH(N = 1000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))
#' df$ipw <- adjIPW(cohort ~ age + mutation, data = df)
#' df_IPW <- survfit(Surv(obs_time, status) ~ cohort, data = df, weights = ipw)

adjIPW <- function(formula, data = NULL) {

  # Deconstruct formula object
  resp <- paste(attr(terms(formula, data = data), which = "variables")[2])
  resp_obs <- factor(model.frame(formula, data = data)[, 1])
  pred <- paste(attr(terms(formula, data = data), which = "term.labels"),
                collapse = " + ")

  # Calculate weights
  ipw_ref <- matrix(NA, nrow = length(resp_obs), ncol = nlevels(resp_obs))
  for(i in 1:nlevels(resp_obs)) {
    ipw_formula <- formula(paste0("I(", resp, "==", "'", levels(resp_obs)[i], "'",")",
                                  "~", pred))
    ipw_ref[, i] <- 1 / predict(glm(ipw_formula, data = data, family = "binomial"),
                                type = "response")
  }

  # Match observed treatment with ipw_ref matrix
  ipw <- numeric(nrow(ipw_ref))
  for(j in 1:nrow(ipw_ref)){
    ipw[j] <- ipw_ref[j, which(resp_obs[j] == levels(resp_obs))]
  }
  return(ipw)
}
