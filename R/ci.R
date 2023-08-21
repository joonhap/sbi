#' @export
ci <- function(x, ...) {
    UseMethod("ci")
}

#' Confidence intervals constructed using simulation log likelihoods
#'
#' `ci` constructs confidence intervals using simulation log likelihoods. See Park (2023) for more information.
#'
#' @name ci
#' @param simll A class 'simll' object, containing simulation log likelihoods, the parameter values at which simulations are made (optional), and the weights for those simulations for regression (optional). See help(simll).
#' @param level A numeric vector of confidence levels.
#' @param ci A character string indicating the quantity for which a confidence interval is to be constructed. Either "MESLE" or "parameter". See Details.
#' @param weights An optional argument. The un-normalized weights of the simulation log likelihoods for regression. A numeric vector of length equal to the 'params' attribute of the 'simll' object. See Details below.
#'
#' @details
#' This is a generic function, taking a class 'simll' object as the first argument.
#' Confidence intervals are constructed under a normal, locally quadratic meta model where the simulation log likelihoods given in the 'simll' object are normally distributed.
#'
#' When 'level' has length greater than one, a confidence interval is constructed for each value in the vector.
#'
#' Quadratic regression for the simulation log likelihoods is carried out to construct confidence intervals, where the x-axis values are the 'params' values of the 'simll' object and the y-axis values are the corresponding simulation log likelihoods.
#' In the case where 'ci' = "parameter", inference on the simulation based surrogate will be carried out under the local asymptotic normality for simulation log likelihoods (see Park (2023) for more information.)
#' The default value of 'ci' is "parameter".
#'
#' If 'ci' = "MESLE", confidence intervals are constructed for the maximum expected simulation likelihood estimate given the observed data.
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified. The length of 'weights' should be equal to that of the 'params' attribute of the 'simll', which is equal to the number of rows in the simulation log likelihood matrix in the 'simll' object. It is important to note that the weights are not supposed to be normalized (i.e., sum to one). Multiplying all weights by the same constant changes the estimation outputs. If not supplied, the 'weights' attribute of the 'simll' object is used. If neither is supplied, 'weights' defaults to the vector of all ones.
#'
#' @return A list consisting of the followings are returned.
#' \itemize{
#' \item{meta_model_MLE_for_*: point estimate for the quantity for which confidence intervals are constructed under a normal meta model,}
#' \item{confidence_interval: a data frame of the lower and upper bounds of the confidence intervals and the corresponding confidence levels. Note that in some unfortunate cases (especially if the quadratic coefficient of the estimated quadratic fit of the log likelihood estimates is close to zero or nonnegative), the confidence interval may be inverted, meaning that it is of the form (-infty, bound1) U (bound2, infty). This case can happen if the signal-to-noise ratio in simulation log likelihoods is too small. The inverted confidence interval will be indicated by the additional column "inverted" in the data frame taking values of 0 or 1.}
#' \item{pval_cubic: The p-value of the test about whether the cubic term in the cubic polynomial regression is significant. If so, the constructed confidence interval may be biased.}
#' }
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @export
ci.simll <- function(simll, level, ci=NULL, weights=NULL) {
    validate_simll(simll)
    if (is.null(ci)) {
        ci <- "parameter"
        message("The 'ci' argument is not supplied. Defaults to 'paramter'.")
    }
    if (!is.null(ci)) {
        match.arg(ci, c("MESLE", "parameter"))
    }
    if (is.null(attr(simll, "params"))) {
        stop("The 'simll' object should have 'params' attributes in order to construct confidence intervals.",
            call. = FALSE)
    }
    if (!is.numeric(level)) {
        stop(
            "'level' should be a numeric vector (length >= 1).",
            call. = FALSE
        )
    }
    ## set weights (vector w)
    if (!is.null(weights)) {
        if (!is.numeric(weights)) {
            stop(
                "When the 'weights' argument is given, it has to be a numeric vector.",
                call. = FALSE
            )
        }
        if (length(weights) != dim(simll)[2]) {
            stop(
                "When the 'weights' argument is given, the length of 'weights' should be equal to the number of rows in 'simll'.",
                call. = FALSE
            )
        }
        if (is.numeric(weights)) {
            w <- weights
        }
    }
    if (is.null(weights)) {
        if (!is.null(attr(simll, "weights"))) {
            if (!is.numeric(attr(simll, "weights"))) {
                stop(
                    "When the 'simll' object has 'weights' attribute, it has to be a numeric vector.",
                    call. = FALSE
                )
            }
            if (dim(simll)[2] != length(attr(simll, "weights"))) {
                stop(
                    "When the 'simll' object has 'weights' attribute, the length of 'weights' should be the same as the number of columns in the simulation log likelihood matrix in 'simll'.",
                    call. = FALSE
                )
            }
            w <- attr(simll, "weights")
        } else {
            w <- rep(1, dim(simll)[2])
        }
    }
    ## weighted quadratic regression
    W  <- diag(w)
    theta <- attr(simll, "params")
    llmat <- unclass(simll)
    ll <- apply(llmat, 2, sum)
    M <- length(ll)
    theta012 <- cbind(1, theta, theta^2)
    Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%ll)) # Ahat=(ahat,bhat,chat)
    resids <- ll - c(theta012%*%Ahat)
    sigsqhat <- c(resids%*%W%*%resids) / M
    theta0123 <- cbind(1, theta, theta^2, theta^3) # cubic regression to test whether the cubic coefficient = 0
    Ahat_cubic <- c(solve(t(theta0123)%*%W%*%theta0123, t(theta0123)%*%W%*%ll))
    resids_cubic <- ll - c(theta0123%*%Ahat_cubic)
    sigsqhat_cubic <- c(resids_cubic%*%W%*%resids_cubic) / M
    pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-4), 1, sum(w>0)-4, lower.tail=FALSE)
    ## ci for MESLE
    if (ci=="MESLE") {
        mtheta1 <- sum(w*theta)/sum(w)
        mtheta2 <- sum(w*theta*theta)/sum(w)
        mtheta3 <- sum(w*theta*theta*theta)/sum(w)
        mtheta4 <- sum(w*theta*theta*theta*theta)/sum(w)
        v11 <- sum(w)*(mtheta2 - mtheta1*mtheta1)
        v12 <- sum(w)*(mtheta3 - mtheta1*mtheta2)
        v22 <- sum(w)*(mtheta4 - mtheta2*mtheta2)
        detV <- v11*v22-v12*v12
        lub <- sapply(1:length(level), function(i) {
            lvl <- level[i]
            q <- qf(lvl, 1, M-3)
            coef2 <- 4*(M-3)*Ahat[3]^2*detV - 4*M*sigsqhat*q*v11
            coef1 <- 4*(M-3)*Ahat[2]*Ahat[3]*detV + 4*M*sigsqhat*q*v12
            coef0 <- (M-3)*Ahat[2]^2*detV - M*sigsqhat*q*v22
            Disc <- coef1^2 - 4*coef2*coef0 # Discriminant
            if (Disc >= 0 && coef2 > 0) { # Case 1
                int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
            } else if (Disc >= 0 && coef2 <= 0) { # Case 2
                int <- -1/(2*coef2)*(coef1+c(-1,1)*sqrt(Disc))
                return(c(level=lvl, lb=int[1], ub=int[2], inverted=1))
            } else { # Case 3
                return(c(level=lvl, lb=-Inf, ub=Inf, inverted=0))
            }
        })
        if (any(lub["inverted",]==1)) { # if for any given level the confidence interval is inverted (Case 2)
            warning(paste0("For level(s) ", toString(unlist(level)[which(lub["inverted",]==1)]), ", the constructed confidence is of the form (-Inf, bound_1) U (bound_2, Inf)."), call.=FALSE)
        } else { # otherwise, remove the "inverted" column
            lub <- lub[-4,]
        }
        out <- list(meta_model_MLE_for_MESLE=c(MESLE=unname(-Ahat[2]/(2*Ahat[3]))),
            confidence_interval=t(lub)
        )
        return(out)
    }
    ## ci for the simulation surrogate under LAN
    if (ci=="parameter") {
        Winv <- diag(1/w)
        uniquetheta <- unique(theta)
        uniquetheta012 <- cbind(1, uniquetheta, uniquetheta^2)
        nobs <- dim(llmat)[1] # number of observations
        slopes <- apply(llmat, 1, function(ll_i) {
            # quadratic regression for each observation piece (each row of simll)
            Ahat_i <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%ll_i))
            return(c(uniquetheta012%*%Ahat_i)) # estimated slopes at uniquetheta
        }) # matrix of dimension length(uniquetheta) times nobs
        var_slope_vec <- apply(slopes, 1, var)
        var_slope_est <- median(var_slope_vec) # estimate of the variance of the slope of the fitted quadratic
        errorvars <- apply(llmat, 1, function(ll_i) {
            Ahat_i <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%ll_i))
            resids_i <- ll_i - c(theta012%*%Ahat_i)
            sigsqhat_i <- c(resids_i%*%W%*%resids_i) / (M-3)
            return(sigsqhat_i) # return estimated error variance for the i-th observation
        }) # vector of length nobs
        MESLE_est <- -Ahat[2]/(2*Ahat[3]) # estimate of the MESLE
        E_condVar_slope <- c(mean(errorvars)*c(0,1,2*MESLE_est)%*%solve(t(theta012)%*%W%*%theta012, c(0,1,2*MESLE_est))) # estimate of the expected value of the conditional variance of the estimated slope given Y (expectation with respect to Y)
        K1hat <- var_slope_est - E_condVar_slope
        if (K1hat <= 0) {
            warning("The estimate of K1 is nonpositive. The results should not be reliable.")
        }
        resids_1 <- ll - c(theta012%*%Ahat) # first stage estimates for residuals
        sigsq_1 <- c(resids_1%*%W%*%resids_1) / M # the first stage estimate of sigma^2
        C <- cbind(-1, diag(rep(1,M-1)))
        Ctheta <- c(C%*%theta)
        Cthetasq <- c(C%*%theta^2)
        Q_1 <- solve(C%*%Winv%*%t(C) + K1hat/sigsq_1^2*outer(Ctheta,Ctheta))
        svdQ_1 <- svd(Q_1)
        sqrtQ_1 <- svdQ_1$u %*% diag(sqrt(svdQ_1$d)) %*% t(svdQ_1$v)
        R_1 <- sqrtQ_1%*%cbind(Ctheta, Cthetasq)
        estEq_2 <- c(solve(t(R_1)%*%R_1, t(R_1)%*%(sqrtQ_1%*%C%*%ll))) # estimating equation for K2 and theta_star: Khat(thetastarhat // -1/2) = (R_1^T R_1)^{-1} R_1^T (Q_1^{1/2} C lhat)
        K2hat <- -2*estEq_2[2]/nobs # second stage estimate of K
        thetastarhat <- -estEq_2[1]/estEq_2[2]/2 # maximum meta model likelihood estimate for theta_star (simulation based surrogate)
        sigsqhat <- 1/(M-1)*sum((sqrtQ_1%*%C%*%ll - nobs*K2hat*R_1%*%c(thetastarhat, -1/2))^2)
        RtR <- t(R_1)%*%R_1
        rho11 <- RtR[1,1]
        rho12 <- RtR[1,2]
        rho22 <- RtR[2,2]
        lub <- sapply(1:length(level), function(i) {
            lvl <- level[i]
            q <- qf(lvl, 1, M-3)
            zeta0 <- sum((sqrtQ_1%*%C%*%ll)^2) - (M-1)*sigsqhat*(q/(M-3)+1)
            zeta12 <- t(R_1)%*%sqrtQ_1%*%C%*%ll
            zeta1 <- zeta12[1,1]
            zeta2 <- zeta12[2,1]
            coef2 <- zeta0*rho11 - zeta1^2
            coef1 <- zeta1*zeta2 - zeta0*rho12
            coef0 <- 1/4*(rho22*zeta0 - zeta2^2)
            Disc <- coef1^2 - 4*coef2*coef0 # Discriminant
            if (Disc >= 0 && coef2 > 0) { # Case 1
                int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
            } else if (Disc >= 0 && coef2 <= 0) { # Case 2
                int <- -1/(2*coef2)*(coef1+c(-1,1)*sqrt(Disc))
                return(c(level=lvl, lb=int[1], ub=int[2], inverted=1))
            } else { # Case 3
                return(c(level=lvl, lb=-Inf, ub=Inf, inverted=0))
            }
        })
        if (any(lub["inverted",]==1)) { # if for any given level the confidence interval is inverted (Case 2)
            warning(paste0("For level(s) ", toString(unlist(level)[which(lub["inverted",]==1)]), ", the constructed confidence is of the form (-Inf, bound_1) U (bound_2, Inf)."), call.=FALSE)
        } else { # otherwise, remove the "inverted" column
            lub <- lub[-4,]
        }
        out <- list(meta_model_MLE_for_parameter=c(parameter=thetastarhat, K1=K1hat, K2=K2hat, error_variance=sigsqhat),
            confidence_interval=t(lub)
        )
        return(out)
    }
}
# TODO: for ci="parameter", 'type' should be 'iid' or 'stationary'.
# The stationary case should be implemented.  A test of stationarity can be carried out, and a warning should be given when the test is positive.

