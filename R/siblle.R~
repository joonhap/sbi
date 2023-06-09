## Class siblle
#' Simulation based log likelihood estimator (SIBLLE) class
#'
#' @param llest A numeric vector of log likelihood estimates
#' @param param A numeric vector of one-dimensional parameter values (optional).
#' @returns A class 'siblle' object
#'
#' Constructor for a class 'siblle' object
#' @export
new_siblle <- function(llest, param=NULL) {
    stopifnot(is.numeric(llest))

    structure(
        llest,
        param = param,
        class = "siblle"
    )
}


## Internal validator function for a class 'siblle' object
validate_siblle <- function(x) {
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
                "The length of the 'param' attribute of an siblle object should be equal to the length of the siblle object or 1.",
                call. = FALSE
            )
        }
    }

    x
}


#' Helper function that creates a class siblle object
#' @export
siblle <- function(llest, param=NULL) {
    validate_siblle(new_siblle(llest, param))
}
