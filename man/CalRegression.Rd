% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zipreg.R
\name{CalRegression}
\alias{CalRegression}
\title{Calculates expected value of the outcome for both trt and ctrl seperately}
\usage{
CalRegression(
  design.trt = c(0, 1),
  design.prog.z,
  design.prog.p,
  design.pred.p,
  coef.prog.z,
  coef.prog.p,
  coef.pred.p
)
}
\arguments{
\item{design.trt}{design treatment =c(0,1)}

\item{design.prog.z}{prognostic design matrix of the logistic regression (zero component)}

\item{design.prog.p}{prognostic design matrix of the Poisson regression (Poisson component)}

\item{design.pred.p}{predictive design matrix of Poisson regression (Poisson component)}

\item{coef.prog.z}{coefficients of prognostic design matrix of the logistic regression (zero component)}

\item{coef.prog.p}{coefficients of prognostic covariate of the Poisson regression (Poisson component)}

\item{coef.pred.p}{coefficients of predictive covariate of Poisson regression (Poisson component)}

\item{design.pred.z}{predictive design matrix of logistic regression (zero component)}

\item{coef.pred.z}{coefficients of predictive covariate of logistic regression (zero component)}
}
\value{
expected value of trt and ctrl
}
\description{
Takes in prognostic and predictive covariates and their corresponding coefficients of both Poisson and Zero components
}
