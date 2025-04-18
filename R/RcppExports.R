# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' The SCL distribution
#'
#' Quantile function, distribution function, and random generation for the SCL distribution family. See Park (2025) for information about the SCL distributions.
#'
#' @name SCL
#' @param p vector of probabilities
#' @param q vector of quantiles
#' @param n number of draws
#' @param M the first parameter for the SCL distributions
#' @param k the second parameter for the SCL distribution
#' @param num_error_size The requested size of numerical error for the outputs of qscl and pscl functions, in terms of the estimated standard deviation of the output. For example num_error_size of 0.01 will output values with the standard deviation of approximately equal to 0.01.
#' @param lower logical; if TRUE, probabilities are P(X <= x), otherwise, P(X > x).
#' @param log_p logical; if TRUE, probabilities p are given as log(p).
#' @param force logical; if TRUE, the function will run regardless of how long it will take. If FALSE, the function will ask if you want to continue, stop, or give a new num_error_size value whenever the expected run time is longer than 15 seconds. 
#' @return a list consisting of the numeric vector of quantiles and the num_error_size (numeric) used.
#' @export
qscl <- function(p, M, k, num_error_size = 0.01, lower = TRUE, log_p = FALSE, force = FALSE) {
    .Call(`_sbim_qscl`, p, M, k, num_error_size, lower, log_p, force)
}

#' @rdname SCL
#' @export
pscl <- function(q, M, k, num_error_size = 0.01, lower = TRUE, log_p = FALSE, force = FALSE) {
    .Call(`_sbim_pscl`, q, M, k, num_error_size, lower, log_p, force)
}

#' @rdname SCL
#' @export
rscl <- function(n, M, k) {
    .Call(`_sbim_rscl`, n, M, k)
}

