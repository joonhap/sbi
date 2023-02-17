## Class mclle
#' Monte Carlo log likelihood estimator (MCLLE) class
#'
#' @param llest A numeric vector of log likelihood estimates
#' @param param A numeric vector of one-dimensional parameter values (optional).
#' @returns A class 'mclle' object
#'
#' Constructor for a class 'mclle' object
#' @export
new_mclle <- function(llest, param=NULL) {
    stopifnot(is.numeric(llest))

    structure(
        llest,
        param = param,
        class = "mclle"
    )
}


## Internal validator function for a class 'mclle' object
validate_mclle <- function(x) {
    llest <- c(unclass(x))
    param <- attr(x, "param")

    if (!is.null(param)) {
        if (!is.numeric(param)) {
            stop(
                "The 'param' attribute should be a numeric vector or a NULL.",
                call. = FALSE
            )
        }
        if (!is.null(attr(param, "dim"))) {
            stop(
                "The 'param' attribute should be a numeric vector with no dim (dimension) attribute, or a NULL.",
                call. = FALSE
            )
        }
        if (length(param) != 1 && length(llest) != length(param)) {
            stop(
                "The length of the 'param' attribute of an mclle object should be equal to the length of the mclle object or 1.",
                call. = FALSE
            )
        }
    }

    x
}


#' Helper function that creates a class mclle object
#' @export
mclle <- function(llest, param=NULL) {
    validate_mclle(new_mclle(llest, param))
}
