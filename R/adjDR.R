#' @title Doubly Robust Cumulative Incidence Curves
#'
#' @description Generates covariate adjusted cumulative incidence curves via doubly robust.
#'
#' @param formula Object containing a `Hist()` object as the response, and the strata and covariates as the predictors. Censored values should be indicated by a 0 in the status variable.
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
#' df_DR <- adjDR(Hist(obs_time, status) ~ treatment + x1 + x2, data = df,
#'   strata = "treatment", times = seq(0, 1.5, 0.05))

adjDR <- function(formula, strata, ref = NULL, data = NULL, times = NULL){

  # Extract components from formula object and convert to model matrix
  tmp <- model.frame(formula, data = data)[, -1]
  tmp_hist <- eval(formula[[2]], envir = data)
  tmp <- cbind(data.frame(time = tmp_hist[, "time"],
                          status = getEvent(tmp_hist)), tmp)
  tmp$status <- as.character(tmp$status)
  tmp$status[tmp_hist[, 2] == "0"] <- "0"
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
  dummy <- do.call("rbind",
                   replicate(nlevels(tmp$strata), tmp[, -3],
                             simplify = F))
  dummy <- cbind(strata = rep(levels(tmp$strata), each = nrow(tmp)),
                 dummy)
  dummy <- dummy[, c(-2:-3)]

  if(!is.null(ref)) {
    dummy_ref <- do.call("rbind",
                         replicate(nlevels(tmp$strata),
                                   tmp[tmp$strata == paste0(strata, "=", ref), -3],
                                   simplify = F))
    dummy_ref <- cbind(strata = rep(levels(tmp$strata),
                                    each = nrow(tmp[tmp$strata == paste0(strata, "=", ref),])),
                       dummy_ref)
    dummy_ref <- dummy_ref[, c(-2:-3)]
  }

  # Construct pseudo-observations
  pseudo_obs <- prodlim(Hist(time, status) ~ 1, data = tmp, type = "risk")

  # Take mean curve per strata for each status
  curves <- vector("list", length(df_csc$causes))
  for(i in 1:length(df_csc$causes)){
    csc_pred <- predict(df_csc,
                        newdata = dummy,
                        times = xtime,
                        cause = df_csc$causes[i])
    curves_i <- vector("list", length(df_csc$causes))

    if(!is.null(ref)) {
      csc_pred_ref <- predict(df_csc,
                              newdata = dummy_ref,
                              times = xtime,
                              cause = df_csc$causes[i])
    }

    # Extract pseudo observations
    pseudo_k <- t(jackknife(pseudo_obs,
                            times = xtime,
                            cause = df_csc$causes[i]))

    for(j in 1:nlevels(tmp$strata)){

      # Extract propensity scores
      pred_ipw <- paste(names(tmp)[4:ncol(tmp)], collapse = " + ")
      ipw_formula <- formula(paste0("I(strata == levels(strata)[j])", "~" ,
                                    pred_ipw))
      ipw_model <- glm(ipw_formula, data = tmp, family = "binomial")
      pscores <- predict(ipw_model, type = "response")

      # Extract ROS estimator
      index <- 1:(nrow(dummy) / nlevels(tmp$strata)) +
        (j - 1) * (nrow(dummy) / nlevels(tmp$strata))
      I_ROS <- t(csc_pred$absRisk[index, ])

      if(!is.null(ref)) {
        index_ref <- 1:(nrow(dummy_ref) / nlevels(tmp$strata)) +
          (j - 1) * (nrow(dummy_ref) / nlevels(tmp$strata))
        I_ROS_ref <- t(csc_pred_ref$absRisk[index_ref, ])
      }

      # Calculate doubly robust estimator
      ind_z <- as.numeric(tmp$strata == levels(tmp$strata)[j])
      IPW_est <- matrix(NA, nrow = length(xtime), ncol = ncol(pseudo_k))
      OR_est <- matrix(NA, nrow = length(xtime), ncol = ncol(pseudo_k))
      for(l in 1:length(xtime)){
        IPW_est[l, ] <- ((ind_z * (pseudo_k[l, ] - I_ROS[l, ])) / pscores)
      }
      if(is.null(ref)){
        curves_i[[j]] <- rowMeans(IPW_est) + rowMeans(I_ROS)
      } else {
        curves_i[[j]] <- rowMeans(IPW_est) + rowMeans(I_ROS_ref)
      }
    }
    curves[[i]] <- do.call("c", curves_i)
  }

  # Combine states and construct output
  curves_output <- do.call("cbind", curves)
  curves_output <- cbind(1 - rowSums(curves_output), curves_output)
  curves_output[curves_output < 0] <- 0
  curves_output[curves_output > 1] <- 1
  colnames(curves_output) <- c("(s0)", df_csc$causes)
  output <- cbind(data.frame(time = rep(xtime, times = nlevels(tmp$strata)),
                             strata = factor(rep(levels(tmp$strata),
                                                 each = length(xtime)))),
                  curves_output)
  return(output)
}
