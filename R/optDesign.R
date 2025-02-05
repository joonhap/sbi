#' @export
optDesign <- function(simll, ...) {
    UseMethod("optDesign")
}

#' Find the next optimal design point for simulation-based inference
#'
#' `optDesign` finds the next design point at which simulation should be carried out for approximately best efficiency in a metamodel-based inference. See Park (2025) for more details on this method. It takes a class `simll` object.
#'
#' @name optDesign
#' @param simll A class `simll` object, containing simulation log likelihoods, the parameter values at which simulations are made, and the weights for those simulations for regression (optional). See help(simll).
#' @param n An integer indicating the number of next Design points to be found
#' @param penalty A positive real that determines how heavily deviation from local quadratic approximation should be penalized. Specifically, a design point is assigned a penalization weight of `exp(-penalty*abs(cubic_term)/quadratic_term)`, where the `quadratic_term` and `cubic_term` respectively indicate the estimated second-order and third-order Taylor expansion terms with respect to the estimated maximum of the expected simulated log-likelihood function. The larger this `penalty` coefficient, the closer the next optimal design points are found near the currently estimated maximizer of the expected simulated log-likelihood. The default value is 2.
#' @param ... Other optional arguments, not currently used.
#'
#' @details
#' This is a generic function, taking a class `simll` object as the first argument.
#' Parameter inference for implicitly defined simulation models can be carried out under a metamodel for the distribution of the log-likelihood estimator.
#' See function `ht` for hypothesis testing and `ci` for confidence interval construction for a one-dimensional parameter.
#' This function `optDesign` founds the next points at which simulations are to be carried out such that the variance of the parameter estimate is reduced approximately the most.
#' This function computes penalization weights given by `exp(-penalty*abs(cubic_term)/quadratic_term)` as described in the explanation for the `penalty` parameter to find design points that are within the scope where local quadratic approximation has relative low error.
#' These penalization weights are multiplied to the original `weights` given to the simulation points specified in the `simll` object.
#' Quadratic regression through the simulated log-likelihoods with these product weights is considered for the selection of next design points.
#'
#' @return A matrix of parameter values at which next simulations are to be carried out for approximately best efficiency.
#' Each row gives a parameter vector.
#'
#' @references Park, J. (2025). Scalable simulation-based inference for implicitly defined models using a metamodel for log-likelihood estimator <https://doi.org/10.48550/arxiv.2311.09446>
#' @export
optDesign.simll <- function(simll, n, penalty=2, ...) {
    validate_simll(simll)
    vech <- function(mat) { # half-vectorization
        if (length(mat)==1) {
            mat <- cbind(mat)
        }
        if (dim(mat)[1] != dim(mat)[2]) {
            stop("The argument to vech should be a square matrix.")
        }
        d <- dim(mat)[1]
        l <- 0
        output <- numeric((d^2+d)/2)
        for (k in 1:d) {
            output[(l+1):(l+d+1-k)] <- mat[k:d,k]
            l <- l+d+1-k
        }
        output
    }
    unvech <- function(vec) { # undo vech
        d <- (-1 + sqrt(1+8*length(vec)))/2
        if (abs(d - round(d)) > 1e-5) {
            stop("The length of the given vector is not equal to that of vech of a symmetric matrix")
        }
        l <- 0
        output <- matrix(0, d, d)
        for (k in 1:d) {
            seg <- vec[(l+1):(l+d+1-k)]
            output[k:d,k] <- seg
            output[k,k:d] <- seg
            l <- l+d+1-k
        }
        output
    }
    matricize <- function(theta) {
        d <- length(theta)
        out <- matrix(0, d, (d^2+d)/2)
        n <- 0
        for (i in 1:d) {
            l <- d+1-i
            out[i:d, (n+1):(n+l)] <- theta[i] * diag(l)
            out[i, (n+1):(n+l)] <- theta[i:d]
            n <- n+l
        }
        out
    }
    vec2 <- function(vec) {
        out <- 2*outer(vec,vec)
        diag(out) <- vec^2
        out
    }
    vec012 <- function(vec) {
        c(1, vec, vech(vec2(vec)))
    }
    ## weighted quadratic regression
    if (!is.null(attr(simll, "weights"))) {
        if (!is.numeric(attr(simll, "weights"))) {
            stop("When the `simll` object has `weights` attribute, it has to be a numeric vector.")
        }
        if (dim(simll)[2] != length(attr(simll, "weights"))) {
            stop("When the `simll` object has `weights` attribute, the length of `weights` should be the same as the number of rows in the simulated log likelihood matrix in `simll`.")
        }
        w <- attr(simll, "weights")
    } else {
        w <- rep(1, dim(simll)[2])
    }
    theta <- cbind(attr(simll, "params")) # coerce into a matrix
    theta_mean <- apply(theta, 2, mean)
    theta_sd <- apply(theta, 2, sd)
    trans_n <- function(vec) { (vec-theta_mean)/theta_sd } # normalize by centering and scaling
    trans_b <- function(vec) { vec*theta_sd + theta_mean } # transform back to the original scale
    theta_n <- apply(theta, 1, trans_n) |> rbind() |> t() # apply trans_n rowwise
    llmat <- unclass(simll)
    ll <- apply(llmat, 2, sum)
    M <- length(ll)
    Theta012 <- t(apply(theta_n, 1, vec012))
    dim012 <- 1 + d + (d^2+d)/2
    WTheta012 <- outer(w,rep(1,dim012))*Theta012
    Ahat <- c(solve(t(Theta012)%*%WTheta012, t(Theta012)%*%(w*ll)))
    ahat <- Ahat[1]
    d <- dim(theta)[2]
    bindex <- 2:(d+1) # the positions in A that correspond to b
    bhat <- Ahat[bindex]
    cindex <- (d+2):((d^2+3*d+2)/2) # the positions in A that correspond to vech(c)
    vech_chat <- Ahat[cindex]
    chat <- unvech(vech_chat)
    resids <- ll - c(Theta012%*%Ahat)
    sigsqhat <- c(resids%*%(w*resids)) / M
    MESLEhat <- unname(-solve(chat,bhat)/2)
    ## cubic test
    if (M > (d+1)*(d+2)*(d+3)/6) { # carry out cubic test if this condition is met
        cubic_test <- TRUE
        vec3 <- function(vec) {
            d <- length(vec)
            l <- 0
            out <- numeric((d^3+2*d^2+d)/6)
            for (k1 in 1:d) {
                for (k2 in 1:k1) {
                    out[(l+1):(l+k2)] <- vec[k1]*vec[k2]*vec[1:k2]
                    l <- l+k2
                }
            }
            out
        }
        Theta0123 <- cbind(Theta012, t(rbind(apply(theta_n, 1, vec3)))) # design matrix for cubic regression to test whether the cubic coefficient = 0
        dim0123 <- dim(Theta0123)[2]
        Ahat_cubic <- c(solve(t(Theta0123)%*%(outer(w,rep(1,dim0123))*Theta0123), t(Theta0123)%*%(w*ll)))
        resids_cubic <- ll - c(Theta0123%*%Ahat_cubic)
        sigsqhat_cubic <- c(resids_cubic%*%(w*resids_cubic)) / M
        pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-(d+1)*(d+2)*(d+3)/6)/(d*(d+1)*(d+2)/6), d*(d+1)*(d+2)/6, sum(w>0)-(d+1)*(d+2)*(d+3)/6, lower.tail=FALSE)
    } else {
        cubic_test <- FALSE
    }

}
