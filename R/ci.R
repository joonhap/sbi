#' @export
ci <- function(x, ...) {
    UseMethod("ci")
}

#' Confidence intervals constructed using simulation log likelihoods
#'
#' `ci` constructs confidence intervals using simulation log likelihoods. See Park (2023) for more information.
#'
#' @name ci
#' @param simll A class 'simll' object, containing simulation log likelihoods, the parameter values at which simulations are made, and the number of simulations made at each parameter value. See help(simll).
#' @param level A numeric vector of confidence levels.
#' @param type A character string indicating what type of situation is considered. The value should be "regression" or "LAN". See Details.
#' @param ci A character string indicating the quantity for which a confidence interval is to be constructed. One of "MESLE", "information", or "parameter". See Details.
#' @param weights An optional argument. The un-normalized weights of the log likelihood estimates for regression. Either a numeric vector of length equal to the 'params' attribute of the 'simll' object. The default weights are equal to the 'nsims' attribute of the 'simll' object. See Details below.
#'
#' @details
#' This is a generic function, taking a class 'siblle' object as the first argument.
#' Confidence intervals are constructed under a normal meta model where the simulation based likelihood estimates given in the 'siblle' object are normally distributed.
#'
#' When 'level' has length greater than one, a confidence interval is constructed for each value in the vector.
#'
#' The 'type' argument should be either "regression" or "LAN".
#' The case 'type' = "regression" means that the 'siblle' object contains simulation based log likelihood estimates evaluated at more than one parameter values, specified by the 'params' attribute of the 'siblle' object. A local quadratic regression for the estimated log likelihood values is used for constructing confidence intervals, where the x-axis values are given by the 'params' values of the 'siblle' object and the y-axis values are the corresponding log likelihood estimates.
#' The case 'type' = "LAN" means that inference on the model parameter will be carried out under the local asymptotic normality (Le Cam and Yang, 2000) condition.
#' If the 'siblle' object has 'params' attribute whose length is equal to the length of the object, then 'type' defaults to "LAN".
#'

#' In this case a confidence interval for the log likelihood given the observed data is constructed.
#'
#' When 'type' = "regression", 'ci' can be "MESLE" or "information".
#' If 'ci' = "MESLE", confidence intervals are constructed for the location of the maximum likelihood estimate given the observed data.
#' If 'ci' = "information", confidence intervals are constructed for the Fisher information, which is assumed to be equal to (-2) times the value of the quadratic coefficient of the regression polynomial describing the mean of the SIBLLE.
#' When 'type' = "regression", 'ci' = "MESLE" is assumed by default.
#'
#' When 'type' = "LAN", 'ci' can be "MESLE", "information", or "parameter".
#' If 'ci' is "MESLE", or "information", the output is the same as in the case where 'type' = "regression".
#' If 'ci' is "parameter", confidence intervals for the value of the model parameter are constructed under the local asymptotic normality assumption.
#' When 'type' = "LAN", 'ci' = "parameter" is assumed by default.
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be supplied.  The weights can either be given as an attribute 'weights' of the 'siblle' object, or as a function argument 'weights', with the latter being used when both are supplied. In either case, 'weights' should be a numeric vector of length equal to that of 'siblle'. If 'weights' is given as a function argument, it can be specified alternatively as a character string "tricube". In this case, the tricube weight (see Cleveland, 1979) is used, and the specified 'fraction' of the points will have nonzero weights. The 'center' argument determines at which parameter value the tricube weight takes the maximum. If weights are not supplied in either location, all weights are taken to be equal to 1.
#' It is important to note that the weights are not supposed to be normalized (i.e., sum to one). Multiplying all weights by the same constant changes the local regression results. Roughly speaking, the variance of the error in the simulation based log likelihood estimate is assumed to be sigma^2/(the weight for the point). See Park (2023) for more information.
#'
#' @return A list consisting of the followings are returned.
#' \itemize{
#' \item{meta_model_MLE_for_*: simulation based maximum likelihood estimate under an appropriate normal meta model,}
#' \item{confidence_interval: a data frame of the lower and upper bounds of the confidence intervals and the corresponding confidence levels. Note that in some unfortunate cases (especially if the quadratic coefficient of the estimated quadratic fit of the log likelihood estimates is close to zero or nonnegative), the confidence interval may be inverted, meaning that it is of the form (-infty, bound1) U (bound2, infty). This case can happen if the signal-to-noise ratio in simulation based log likelihood estimates is too small. The inverted confidence interval will be indicated by the additional column "inverted" in the data frame taking values of 0 or 1.}
#' }
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @references Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American statistical association, 74(368), 829-836.
#' @references Le Cam, L. and Yang, G. L. (2000). Asymptotics in statistics: some basic concepts. Springer-Verlag, New York.
#' @export
ci.simll <- function(siblle, level, type=NULL, ci=NULL, weights=NULL, fraction=NULL, center=NULL) {
    validate_siblle(siblle)
    if (!is.null(type)) {
        match.arg(type, c("regression", "LAN"))
    }
    if (is.null(attr(siblle, "params"))) {
        stop("The 'siblle' object should have 'params' attributes for constructing confidence intervals.",
            call. = FALSE)
    }
    if (is.null(type) && !is.null(attr(siblle, "params")) && length(siblle) == length(attr(siblle, "params"))) {
        type <- "LAN"
    }

    if (!is.null(ci)) {
        match.arg(ci, c("MESLE", "information", "parameter"))
    }
    if (is.null(ci)) {
        if (type=="regression") {
            ci <- "MESLE"
        }
        if (type=="LAN") {
            ci <- "parameter"
        }
    }

    if (!is.numeric(level)) {
        stop(
            "'level' should be a numeric vector (length >= 1).",
            call. = FALSE
        )
    }
    if (type=="regression" && !ci %in% c("MESLE", "information")) {
        stop(
            "When 'type' = 'regression', 'ci' should be one eiher 'MESLE' (default) or 'information'.",
            call. = FALSE
        )
    }
    if (type=="LAN" && !ci %in% c("MESLE", "information", "parameter")) {
        stop(
            "When 'type' = 'LAN', 'ci' should be one of 'MESLE', 'information', or 'parameter' (default).",
            call. = FALSE
        )
    }
    if (type!="LAN" && ci=="parameter") {
        stop(
            "'ci' = 'parameter' is only available when 'type' = 'LAN'.",
            call. = FALSE
        )
    }

    if (type %in% c("regression", "LAN")) {
        ## set weights (vector w)
        if (!is.null(weights)) {
            if (!is.numeric(weights) && weights!="tricube") {
                stop(
                    "When 'type' = 'regression' or 'LAN' and the 'weights' argument is given, 'weights' have to be a numeric vector or a string 'tricube'.",
                    call. = FALSE
                )
            }
            if (length(weights) != length(siblle) && weights!="tricube") {
                stop(
                    "When 'type' = 'regression' or 'LAN' and the 'weights' argument is given, the length of 'weights' should be the same as the length of 'siblle'.",
                    call. = FALSE
                )
            }
            if (weights=="tricube") {
                if (is.null(fraction) || !is.numeric(fraction) || length(fraction) != 1) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and 'weights' = 'tricube', 'fraction' has to be given as a numeric value.",
                        call. = FALSE
                    )
                }
                if (is.null(center) || !is.numeric(center) || length(center) != 1) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and 'weights' = 'tricube', 'center' has to be given as a numeric value.",
                        call. = FALSE
                    )
                }
                distance <- abs(attr(siblle, "params") - center)
                span <- distance[order(distance)[ceiling(fraction*length(siblle))]]
                tricube <- function(x) { pmax((1-abs(x)^3)^3, 0) }
                w <- tricube(distance / span)
            }
            if (is.numeric(weights)) {
                w <- weights
            }
        }
        if (is.null(weights)) {
            if (!is.null(attr(siblle, "weights"))) {
                if (!is.numeric(attr(siblle, "weights"))) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and the 'siblle' object has 'weights' attribute, it has to be a numeric vector.",
                        call. = FALSE
                    )
                }
                if (length(siblle) != length(attr(siblle, "weights"))) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and the 'siblle' object has 'weights' attribute, the length of 'weights' should be the same as the length of 'siblle'.",
                        call. = FALSE
                    )
                }
                w <- attr(siblle, "weights")
            } else {
                w <- rep(1, length(siblle))
            }
        }
        ## weighted quadratic regression
        W  <- diag(w)
        theta <- attr(siblle, "params")
        llest <- c(unclass(siblle))
        M <- length(llest)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest)) # Ahat=(ahat,bhat,chat)
        resids <- llest - c(theta012%*%Ahat)
        sigsqhat <- c(resids%*%W%*%resids) / M
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
        ## ci for Fisher information
        if (ci=="information") {
            U <- t(theta012)%*%W%*%theta012
            u3gv12 <- U[3,3] - c(U[3,1:2]%*%solve(U[1:2,1:2], U[1:2,3])) # u_{3|12}
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[i]
                q <- qf(lvl, 1, M-3)
                int <- -2*Ahat[3] + c(-1,1)*2*sqrt(M*q*sigsqhat/(M-3)/u3gv12)
                return(c(level=lvl, lb=max(0,int[1]), ub=max(0,int[2])))
            })
            out <- list(meta_model_MLE_for_Fisher_information=c(Fisher_information=unname(-2*Ahat[3])),
                confidence_interval=t(lub)
            )
            return(out)
        }
        ## ci for the model parameter under LAN
        if (type=="LAN" && ci=="parameter") {
            Winv <- diag(1/w)
            K_1 <- -2*Ahat[3] # first stage estimate of K
            resids_1 <- llest - c(theta012%*%Ahat) # first stage estimates for residuals
            sigsq_1 <- c(resids_1%*%W%*%resids_1) / M # the first stage estimate of sigma^2
            C <- cbind(-1, diag(rep(1,M-1)))
            Ctheta <- c(C%*%theta)
            Cthetasq <- c(C%*%theta^2)
            Q_1 <- solve(C%*%Winv%*%t(C) + K_1/sigsq_1^2*outer(Ctheta,Ctheta))
            svdQ_1 <- svd(Q_1)
            sqrtQ_1 <- svdQ_1$u %*% diag(sqrt(svdQ_1$d)) %*% t(svdQ_1$v)
            R_1 <- sqrtQ_1%*%cbind(Ctheta, Cthetasq)
            u <- c(sqrtQ_1 %*% C %*% Winv %*% rep(1,M))
            estEq_2 <- c(solve(t(R_1)%*%R_1, t(R_1)%*%(sqrtQ_1%*%C%*%llest + sigsq_1/2*u))) # estimating equation for the second stage estimate Khat and thetahat: Khat(thetahat // -1/2) = (R_1^T R_1)^{-1} R_1^T (Q_1^{1/2} C lhat + 1/2*sigma^2_1*u)
            Khat <- -2*estEq_2[2] # second stage estimate of K
            thetahat <- estEq_2[1]/Khat # maximum meta model likelihood estimate for theta
            sigsqhat <- 1/(M-1)*sum((sqrtQ_1%*%C%*%llest - Khat*R_1%*%c(thetahat, -1/2) + sigsq_1/2*u)^2)
            RtR <- t(R_1)%*%R_1
            rho11 <- RtR[1,1]
            rho12 <- RtR[1,2]
            rho22 <- RtR[2,2]
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[i]
                q <- qf(lvl, 1, M-3)
                zeta0 <- sum((sqrtQ_1%*%C%*%llest + sigsq_1/2*u)^2) - (M-1)*sigsqhat*(q/(M-3)+1)
                zeta12 <- t(R_1)%*%(sqrtQ_1%*%C%*%llest + sigsq_1/2*u)
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
            out <- list(meta_model_MLE_for_parameter=c(parameter=thetahat, information=Khat, error_variance=sigsqhat),
                confidence_interval=t(lub)
            )
            return(out)
        }
    }
}
