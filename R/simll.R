## Class simll
#' Simulation Log Likelihood class
#'
#' @param ll A matrix vector of simulation log likelihoods. The (m,i)-th entry is given by the log density of the observation y_i given the simulation X(theta_m).
#' @param params A numeric vector of one-dimensional parameter values (optional, can be omitted if simulation log likelihoods are obtained at only one parameter value.)
#' @param weights A numeric vector of weights, inversely proportional to the variance of simulation log likelihoods (optional)
#' @returns A class 'sll' object
#'
#' @export
simll <- function(ll, params=NULL, weights=NULL) {
    validate_simll(new_simll(ll, params))
}


#' Constructor for a class 'simll' object
#'
#' Constructs a new class 'simll' object. Note that the 'simll' function uses the 'new_simll' function to construct an object and use an interval validator function to check the correctness of the object specification.
#' @export
new_simll <- function(ll, params=NULL, weights=NULL) {
    stopifnot(is.numeric(ll))

    structure(
        ll,
        params = params,
        class = "simll"
    )
}


## Internal validator function for a class 'simll' object
validate_simll <- function(x) {
    ll <- unclass(x)
    params <- attr(x, "params")
    weights <- attr(x, "weights")

    if (is.null(dim(ll))) { # if ll is not a matrix
        ll <- matrix(ll, nrow=1) # coerce into a matrix
        message("Simulation log likelihoods were coerced into a matrix form.")
    }

    if (length(dim(ll)) != 2) {
        stop(
            "The simulation log likelihoods should be a matrix",
            call. = FALSE
        )
    }

    if (!is.null(params)) {
        if (!is.numeric(params)) {
            stop(
                "The 'params' attribute should be a numeric vector or a NULL.",
                call. = FALSE
            )
        }
        if (!is.null(attr(params, "dim"))) {
            stop(
                "The 'params' attribute should be a numeric vector with no dim (dimension) attribute (i.e., not a matrix or an array), or a NULL.",
                call. = FALSE
            )
        }
        if (length(params) != 1 && dim(ll)[1] != length(params)) {
            stop(
                "The length of the 'params' attribute of an simll object should be equal to the number of rows in the matrix of simulation log likelihoods, or 1.",
                call. = FALSE
            )
        }
    }

    if (!is.null(weights)) {
        if (!is.numeric(weights)) {
            stop(
                "The 'weights' attribute should be a numeric vector or a NULL.",
                call. = FALSE
            )
        }
        if (!is.null(attr(weights, "dim"))) {
            stop(
                "The 'weights' attribute should be a numeric vector with no dim (dimension) attribute (i.e., not a matrix or an array), or a NULL.",
                call. = FALSE
            )
        }
        if (length(weights) != 1 && dim(ll)[1] != length(weights)) {
            stop(
                "The length of the 'weights' attribute of an simll object should be equal to the number of rows in the matrix of simulation log likelihoods, or 1.",
                call. = FALSE
            )
        }
    }

    x
}

