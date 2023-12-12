## Class simll
#' Simulation Log Likelihood class
#'
#' @param ll A matrix of simulation log likelihoods. The (i,m)-th entry is given by the simulation log likelihood for y_i obtained by simulating X at theta_m (e.g., the log density of y_i given X).
#' @param params A matrix or a vector of parameter values. If a matrix, the m-th row gives the parameter vector theta_m. If theta is one dimensional, 'params' is can be a numeric vector or a matrix with one column. 'params' can be omitted if simulation log likelihoods are obtained at a single one parameter value.
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
        ll <- matrix(ll, ncol=1) # coerce into a matrix
        x <- new_simll(ll)
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
                "The 'params' attribute should be a numeric matrix, a vector or a NULL.",
                call. = FALSE
            )
        }
        if (!is.null(attr(params, "dim"))) {
            if (dim(params)[1] != dim(ll)[2]) {
                stop(
                    "The number of rows in the 'params' matrix should be equal to the number of columns in the 'll' matrix.",
                    call. = FALSE
                )
            }
        } else {
            if (length(params) != dim(ll)[2]) {
                stop(
                    "The length of the 'params' vector should be equal to the number of columns in the 'll' matrix.",
                    call. = FALSE
                )
            }
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
        if (length(weights) != 1 && dim(ll)[2] != length(weights)) {
            stop(
                "The length of the 'weights' attribute of an simll object should be equal to the number of columns in the matrix of simulation log likelihoods, or 1.",
                call. = FALSE
            )
        }
    }

    x
}

