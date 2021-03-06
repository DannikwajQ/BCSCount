% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zipreg.R
\name{credible.bands}
\alias{credible.bands}
\title{credible.bands calculates the credible bands}
\usage{
credible.bands(
  params,
  design = NULL,
  FUN = function(x, params) {     params \%*\% t(x) },
  cred.level = 0.95,
  method = c("asymptotic", "quantile"),
  sides = c("both", "upper", "lower"),
  est.FUN = mean,
  var.FUN = sd,
  point.estimate = NULL,
  track = numeric(0),
  verbose = FALSE
)
}
\description{
credible.bands calculates the credible bands
}
