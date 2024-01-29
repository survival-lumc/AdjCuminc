#' @title Extract cumulative incidence curves from a survfit object
#'
#' @description This helper function can be used to quickly extract the strata, times, and cause-specific probabilities from a survfit object for each event in a long format, such that it can be easily plotted (eg. by ggplot).
#'
#' @param object  a survfit object from which one would like to extract the cumulative incidence estimator.
#' @param strata  a boolean operator (TRUE / FALSE), which determines whether to return the marginal probabilities per event. If TRUE, which is by default, the cause-specific probabilities are returned as stacked probabilities for subsequent events.
#'
#' @return A `data.frame` with the output from survfit in long format
#' @export
#' @import survival
#' @importFrom tidyr gather
#'
#' @examples
#' df <- simCSH(N = 1000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))
#' df_long <- longCSH(survfit(Surv(obs_time, status) ~ treatment, data = df))

longCSH <- function(object, stack = F) {

  # Exctract components from survfit
  sf_CI <- object$pstate
  time <- object$time
  strata <- rep(names(object$strata), object$strata)
  if(is.null(strata)){
    strata <- rep(1, length(time))
  }
  tmp <- data.frame(time, strata = factor(strata), sf_CI)
  colnames(tmp)[3:ncol(tmp)] <- object$states

  # Return subsequent events as a stacked probability
  if(stack == T){
    tmp2 <- tmp
    for(i in 4:ncol(tmp)){
      tmp2[i] <- rowSums(tmp[3:i])
    }
    tmp <- tmp2
  }

  # Construct long format and return
  tmp_long <- tidyr::gather(tmp, event, prob, 3:ncol(tmp))
  return(tmp_long)
}
