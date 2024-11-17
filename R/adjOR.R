#' @title Outcome Regression Adjusted Cumulative Incidence Curves
#'
#' @description Generates covariate adjusted cumulative incidence curves via outcome regression.
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
#' @importFrom stats model.frame terms
#'
#' @examples
#' df <- simCSH(N = 1000, k = 2,
#'   theta_coef = c(-1, -0.5),
#'   beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
#'   base_haz = c(2, -2))
#' df_OR <- adjOR(Hist(obs_time, status) ~ treatment + x1 + x2, data = df,
#'   strata = "treatment", times = seq(0, 1.5, 0.05))

adjOR <- function(formula, strata, ref = NULL, data = NULL, times = NULL){
  
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
  
  # Take mean curve per strata for each status
  curves <- vector("list", length(df_csc$causes))
  for(i in 1:length(df_csc$causes)){
    csc_pred <- predict(df_csc, 
                        newdata = dummy, 
                        times = xtime, 
                        cause = df_csc$causes[i])
    curves_i <- vector("list", nlevels(tmp$strata))
    for(j in 1:nlevels(tmp$strata)){
      index <- 1:(nrow(dummy) / nlevels(tmp$strata)) + 
        (j - 1) * (nrow(dummy) / nlevels(tmp$strata))
      curves_i[[j]] <- colMeans(csc_pred$absRisk[index, ])
    }
    curves[[i]] <- do.call("c", curves_i)
  }
  
  # Combine states and construct output
  curves_output <- do.call("cbind", curves)
  curves_output <- cbind(1 - rowSums(curves_output), curves_output)
  colnames(curves_output) <- c("(s0)", df_csc$causes)
  output <- cbind(data.frame(time = rep(xtime, times = nlevels(tmp$strata)),
                             strata = factor(rep(levels(tmp$strata), 
                                                 each = length(xtime)))),
                  curves_output)
  return(output)
}
