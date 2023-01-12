## Class mclle
#' Monte Carlo log likelihood estimator (MCLLE) class
#'
#' @param llest A numeric vector of log likelihood estimates
#' @param param A numeric vector of one-dimensional parameter values (optional).
#' @returns An mclle object
#'
#' Constructor for a class mclle object
new_mclle <- function(llest, param=NULL) {
    stopifnot(is.numeric(llest))

    structure(
        llest,
        param = param,
        class = "mclle"
    )
}


## Internal validator function for a class mclle object
validate_mclle <- function(x) {
    llest <- unclass(x)
    param <- attr(x, "param")

    if (!is.null(param) && !is.numeric(param)) {
        stop(
            "The 'param' attribute should be a numeric vector or a NULL.",
            call. = FALSE
        )
    }

    if (!is.null(param) && length(param) != 1 && length(llest) != length(param)) {
        stop(
            "The length of the 'param' attribute of an mclle object should be equal to the length of the mclle object or 1.",
            call. = FALSE
        )
    }

    x
}


#' Helper function that creates a class mclle object
mclle <- function(llest, param=NULL) {
    validate_mclle(new_mclle(llest, param))
}


#' Hypothesis tests using Monte Carlo log likelihood estimates




#' Construct a confidence interval using Monte Carlo log likelihood estimates
