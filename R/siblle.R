## Class siblle
#' Simulation based log likelihood estimator (SIBLLE) class
#'
#' @param llest A numeric vector of log likelihood estimates
#' @param params A numeric vector of one-dimensional parameter values (optional, used if the log likelihood estimates are obtained at more than one parameter values).
#' @returns A class 'siblle' object
#'
#' @export
siblle <- function(llest, params=NULL) {
    validate_siblle(new_siblle(llest, params))
}


#' Constructor for a class 'siblle' object
#'
#' Constructs a new class 'siblle' object. Note that the 'siblle' function uses the 'new_siblle' function to construct an object and use an interval validator function to check the correctness of the object specification.
#' @export
new_siblle <- function(llest, params=NULL) {
    stopifnot(is.numeric(llest))

    structure(
        llest,
        params = params,
        class = "siblle"
    )
}


## Internal validator function for a class 'siblle' object
validate_siblle <- function(x) {
    llest <- c(unclass(x))
    params <- attr(x, "params")

    if (!is.null(params)) {
        if (!is.numeric(params)) {
            stop(
                "The 'params' attribute should be a numeric vector or a NULL.",
                call. = FALSE
            )
        }
        if (!is.null(attr(params, "dim"))) {
            stop(
                "The 'params' attribute should be a numeric vector with no dim (dimension) attribute, or a NULL.",
                call. = FALSE
            )
        }
        if (length(params) != 1 && length(llest) != length(params)) {
            stop(
                "The length of the 'params' attribute of an siblle object should be equal to the length of the siblle object or 1.",
                call. = FALSE
            )
        }
    }
    x
}

