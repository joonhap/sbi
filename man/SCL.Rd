% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SCL}
\alias{SCL}
\alias{qscl}
\alias{pscl}
\alias{rscl}
\title{The SCL distribution}
\usage{
qscl(
  p,
  M,
  k,
  num_error_size = 0.01,
  lower = TRUE,
  log_p = FALSE,
  force = FALSE
)

pscl(
  q,
  M,
  k,
  num_error_size = 0.01,
  lower = TRUE,
  log_p = FALSE,
  force = FALSE
)

rscl(n, M, k)
}
\arguments{
\item{p}{vector of probabilities}

\item{M}{the first parameter for the SCL distributions}

\item{k}{the second parameter for the SCL distribution}

\item{num_error_size}{The requested size of numerical error for the outputs of qscl and pscl functions, in terms of the estimated standard deviation of the output. For example num_error_size of 0.01 will output values with the standard deviation of approximately equal to 0.01.}

\item{lower}{logical; if TRUE, probabilities are P(X <= x), otherwise, P(X > x).}

\item{log_p}{logical; if TRUE, probabilities p are given as log(p).}

\item{force}{logical; if TRUE, the function will run regardless of how long it will take. If FALSE, the function will ask if you want to continue, stop, or give a new num_error_size value whenever the expected run time is longer than 15 seconds.}

\item{q}{vector of quantiles}

\item{n}{number of draws}
}
\value{
a list consisting of the numeric vector of quantiles and the num_error_size (numeric) used.
}
\description{
Quantile function, distribution function, and random generation for the SCL distribution family. See Park (2025) for information about the SCL distributions.
}
