#' @title Doubly Robust Cumulative Incidence Curves
#'
#' @description Generates covariate adjusted cumulative incidence curves via doubly robust.
#'
#' @param formula Object containing a `Surv()` object as the response, and the strata and covariates as the predictors.
#' @param strata A character string used to identify the strata within the `formula` object.
#' @param ref A character input to specify a reference group from `strata`. If this argument is `NULL`, which is by default, the entire sample is used as a reference instead.
#' @param data An optional `data.frame` used to interpret the variables described by `formula` argument.
#' @param times A series of user-defined time points for which the cumulative incidence should be computed per strata. If `NULL`, which is by default, the unique observation times from the data are used instead.
#'
#' @return A `data.frame` containing the covariate adjusted cumulative incidence curves per strata
#' @export
#' @import riskRegression
#' @import survival
#' @importFrom prodlim jackknife prodlim
#' @importFrom stats model.frame terms
#'
#' @examples
#' df <- simCSH(N = 1000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))
#' df_DR <- adjDR(Surv(obs_time, status) ~ treatment + x1 + x2, data = df,
#'   strata = "treatment", times = seq(0, 1.5, 0.05))

adjDR <- function(formula, strata, ref = NULL, data = NULL, times = NULL){

  # Extract components from formula object and convert to model matrix
  tmp <- model.frame(formula, data)
  status_levels <- levels(tmp[, 1])
  tmp <- cbind(data.frame(time = tmp[, 1][, 1],
                          status = factor(tmp[, 1][, 2])), tmp[2:ncol(tmp)])
  names(tmp)[which(names(tmp) == strata)] <- "strata"
  tmp$strata <- factor(paste0(strata, "=", tmp$strata))
  tmp <- tmp[, c(which(colnames(tmp) %in% c("time", "status", "strata")),
                 which(!colnames(tmp) %in% c("time", "status", "strata")))]
  names(tmp) <- gsub(".*\\$", "", names(tmp))

  # Produce csc model
  pred <- paste(names(tmp)[3:ncol(tmp)], collapse = " + ")
  cox_formula <- formula(paste0("Hist(time, status)", "~" , pred))
  df_csc <- CSC(cox_formula, data = tmp)

  # Prediction times
  if(is.null(times)){
    xtime <- sort(unique(c(0, tmp$time)))
  } else {
    xtime <- times
  }

  # Construct dummy data w/ every observation for each strata
  if(is.null(ref)){
    dummy <- do.call("rbind",
                     replicate(nlevels(tmp$strata), tmp[, -3],
                               simplify = F))
    dummy <- cbind(strata = rep(levels(tmp$strata), each = nrow(tmp)),
                   dummy)
  } else {
    dummy <- do.call("rbind",
                     replicate(nlevels(tmp$strata),
                               tmp[tmp$strata == paste0(strata, "=", ref), -3],
                               simplify = F))
    dummy <- cbind(strata = rep(levels(tmp$strata),
                                each = nrow(tmp[tmp$strata == paste0(strata, "=", ref),])),
                   dummy)
  }
  dummy <- dummy[, c(-2:-3)]

  # Construct pseudo-observations
  pseudo_obs <- prodlim(Hist(time, status) ~ 1, data = tmp, type = "risk")

  # Take mean curve per strata for each status
  curves <- vector("list", nlevels(tmp$status) - 1)
  for(i in 1:(nlevels(tmp$status) - 1)){
    csc_pred <- predict(df_csc,
                        newdata = dummy,
                        times = xtime,
                        cause = levels(tmp$status)[-1][i])
    curves_i <- vector("list", nlevels(tmp$strata))

    # Extract pseudo observations
    pseudo_k <- t(jackknife(pseudo_obs,
                            times = xtime,
                            cause = levels(tmp$status)[-1][i]))
    if(!is.null(ref)) {
      pseudo_k <- pseudo_k[, tmp$strata == paste0(strata,"=",ref)]
    }

    for(j in 1:nlevels(tmp$strata)){

      # Extract propensity scores
      pred_ipw <- paste(names(tmp)[4:ncol(tmp)], collapse = " + ")
      ipw_formula <- formula(paste0("I(strata == levels(strata)[j])", "~" ,
                                    pred_ipw))
      ipw_model <- glm(ipw_formula, data = tmp, family = "binomial")
      pscores <- predict(ipw_model, type = "response")
      if(!is.null(ref)) {
        pscores <- pscores[tmp$strata == paste0(strata, "=", ref)]
      }

      # Extract ROS estimator
      index <- 1:(nrow(dummy) / nlevels(tmp$strata)) +
        (j - 1) * (nrow(dummy) / nlevels(tmp$strata))
      I_ROS <- t(csc_pred$absRisk[index, ])

      # Calculate doubly robust estimator
      ind_z <- as.numeric(tmp$strata == levels(tmp$strata)[j])
      if(!is.null(ref)) {ind_z <- ind_z[tmp$strata == paste0(strata, "=", ref)]}
      DR_est <- matrix(NA, nrow = length(xtime), ncol = ncol(pseudo_k))
      for(l in 1:length(xtime)){
        DR_est[l, ] <- ((ind_z * pseudo_k[l, ]) -
                          ((ind_z - pscores) * I_ROS[l, ])) /
          pscores
      }
      curves_i[[j]] <- rowMeans(DR_est)
    }
    curves[[i]] <- do.call("c", curves_i)
  }

  # Combine states and construct output
  curves_output <- do.call("cbind", curves)
  curves_output <- cbind(1 - rowSums(curves_output), curves_output)
  colnames(curves_output) <- c("(s0)", status_levels)
  output <- cbind(data.frame(time = rep(xtime, times = nlevels(tmp$strata)),
                             strata = factor(rep(levels(tmp$strata),
                                                 each = length(xtime)))),
                  curves_output)
  return(output)
}
