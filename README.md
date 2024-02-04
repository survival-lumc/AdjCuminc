# Covariate Adjustment on Cumulative Incidence Curves for Competing Risks

The repository contains the following functions related to covariate adjusted for competing risks:


+ [R/adjOR.R](R/adjOR.R) : A function to produce to perform covariate adjustment of cumulative incidence curves via outcome regression with standardization using cause-specific Cox proportional hazard models
+ [R/adjDR.R](R/adjDR.R) : A function to produce to perform covariate adjustment of cumulative incidence curves via doubly robust estimation using cause-specific Cox proportional hazard models
+ [R/adjIPW.R](R/adjIPW.R): A function which produces observation weights by taking the inverse of probability of allocated treatment estimated by standard logistic regression

The following scripts were used to validate these methods in a simulation study:

+ [adj_comprisk_simstudy.Rmd](adj_comprisk_simstudy.Rmd) : The Rmd-script used to validate the correction methods. Produces the tables and figures of the manuscript related to the simulation study.
+	[R/longCSH.R](R/longCSH.R) : A helper function which translates output from survfit to a long format data frame
+	[R/simCSH.R](R/simCSH.R): Simulation of competing risks data from a cause-specific hazard model as a function of random strata allocation, a discrete covariate, and a continuous covariate
+ [R/trueCSH.R](R/trueCSH.R): A function to analytically obtain the hypothetical true cause-specific hazard curve as a function of strata allocation and covariates 
+ [R/confCSH.R](R/confCSH.R): Simulation of competing risk data from cause-specific hazard models with biased strata allocation as a function of a discrete and continuous covariate

## Installation and usage

Git users 

```bash
git clone https://github.com/survival-lumc/AdjCuminc.git
```

Alternatively, a ZIP file can be found at the top-right corner of this Github page, which can be stored and extracted a directory of choice. 
The code can be accesed by opening the `adjCuminc.Rproj` file in Rstudio. 

## Contributions

| Name                                                         | Affiliation                           | Role                                            |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------------------------------- |
| [Patrick van Hage](https://github.com/pvanhage/)  | Institute of Biology Leiden (IBL) | Author .Rmd files/maintainer                    |
| [Hein Putter](https://www.universiteitleiden.nl/en/staffmembers/hein-putter) | Leiden University Medical Center (NL) | Review of both .R and .Rmd scripts              |
| [Nan van Geloven](https://www.universiteitleiden.nl/medewerkers/nan-van-geloven) | Leiden University Medical Center (NL) | Review of both .R and .Rmd scripts              |
| [Saskia le Cessie](https://www.universiteitleiden.nl/medewerkers/saskia-le-cessie) | University Medical Centre Utrecht (NL)  | Review of both .R and .Rmd scripts              |
