% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confCSH.R
\name{confCSH}
\alias{confCSH}
\title{Confounded sampling from simulated competing risk data under the Cause-Specific Hazards Model}
\usage{
confCSH(
  N,
  k,
  theta_coef,
  beta_coef,
  sample_coef,
  control_ref = FALSE,
  base_haz
)
}
\arguments{
\item{N}{the number of observations to be generated.}

\item{k}{the number of competing outcomes to be generated.}

\item{theta_coef}{a vector where each element represents the coefficient of treatment on the hazard of each event.}

\item{beta_coef}{a list where the number of elements is equal to the number of competing outcomes. Aach listed element is a vector of covariate coefficients associated with a competing outcome. The first element in the vector is the effect size of the cohort, the second element the effect of the discrete covariate, and the third element the effect of the continuous covariate.}

\item{sample_coef}{a vector of length 3 that introduce non-random assignment to a treatment group. If any of these elements are are not 0, scenario 1 is induced, where the first two elements affect probability of assignment conditional on the discrete covariate, and the third element affects the assignment conditional on the continuous covariate.}

\item{control_ref}{a logical indicating whether control observations should be resampled under scenario 0. If `TRUE`, Scenario 2 is induced where the covariates for `control`-labeled observations are are resampled under scenario 0 (prior to event time simulation).}

\item{base_haz}{a baseline hazard for each competing outcome.}
}
\value{
A `data.frame` containing simulated observation times, a status indicator, treatment assignments, a discrete covariate (x1), and a continuous covariate (x2)
}
\description{
These functions were created to simulate a biased sample from a cause-specific hazards model. The formulation of the outcome model, from which the observations are sampled. Three variants have been introduced to simulate different scenarios of model misspecification. Adjustment of the regression coefficients in beta_coef determines the strength and direction of the association between the covariates and treatment allocation, where beta_coef = 0 results in strict random selection on a covariate, beta_coef > 0 results in positive selection for higher covariate values, and beta_coef > 0 results in positive selection for lower covariate values. The association of these regression coefficients on the control treatment is reversed.

The following scenario's can be simulated with this function:

Scenario 1: introduces a logistic association between treatment allocation and the discrete and continuous covariates.

Scenario 2: incorporates a non-logistic association by sampling individuals under the treatment group under scenario 1, and the control group from scenario 0. This can be achieved by setting control_ref = TRUE.

Scenario 3: introduces scenario 1, but incorporates a non-linear effect by adjusting the baseline hazard conditional under the continuous covariate. This can be induced by using two different values for base_haz.
}
\examples{
confCSH(N = 4000, k = 2,
  theta_coef = c(-1, -0.5),
  beta_coef = list(c(1, -1, 0.5), c(-1, 1, -0.5)),
  sample_coef = c(1, -1, 1),
  control_ref = FALSE,
  base_haz = c(0, 0))
}
