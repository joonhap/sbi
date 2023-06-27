#' @export
ci <- function(x, ...) {
    UseMethod("ci")
}

#' Confidence intervals constructed using log likelihood estimates in simulation-based inference
#'
#' `ci` constructs conservative confidence intervals using simulation based log likelihood estimates. See Park and Won (2023) for more information.
#'
#' @name ci
#' @param siblle A class 'siblle' object, containing the simulation based log likelihood estimates and, if those estimates were obtained for different parameter values, the parameter values at which the log likelihood estimates were obtained.
#' @param level A numeric vector of confidence levels.
#' @param type A character string indicating what type of situation is considered. One of "point", "regression", or "LAN". See Details.
#' @param ci A character string indicating the quantity for which a confidence interval is to be constructed. One of "loglik", "MLE", "information", or "parameter". See Details.
#' @param at.param If 'type' = "regression" or "LAN" and 'ci' = "loglik", a confidence interval is constructed for the value of the log likelihood function evaluated at 'at.param'.
#' @param weight.at.param The relative inverse variance for the simulation based log likelihood estimator at parameter value 'at.param'. The weight for regression is proportional to the inverse variance. This argument is used for the cases 'type' = "regression" or "LAN" and 'test' = "loglik". The default value is 1.
#' @param weights An optional argument for the un-normalized weights of the simulation based log likelihood estimates for regression. Either a numeric vector of length equal to the 'siblle' object, or a character string equal to "tricube". The default weights are equal to one for all the points if not specified here or in the siblle object.
#' @param fraction An optional argument indicating the fraction of points with nonzero weights for the case where 'weights' is specified as "tricube".
#' @param center An optional argument indicating the center of the local regression for the case where 'weights' is specified as "tricube".
#'
#' @details
#' This is a generic function, taking a class 'siblle' object as the first argument.
#' Confidence intervals are constructed under a normal meta model where the simulation based likelihood estimates given in the 'siblle' object are normally distributed.
#'
#' When 'level' has length greater than one, a confidence interval is constructed for each value in the vector.
#'
#' The 'type' argument should be one of "point", "regression", or "LAN".
#' The case 'type' = "point" means that the 'siblle' object contains simulation based log likelihood estimates for a single, fixed parameter value.
#' The case 'type' = "regression" means that the 'siblle' object contains simulation based log likelihood estimates evaluated at more than one parameter values, specified by the 'params' attribute of the 'siblle' object. A local quadratic regression for the estimated log likelihood values is used for constructing confidence intervals, where the x-axis values are given by the 'params' values of the 'siblle' object and the y-axis values are the corresponding log likelihood estimates.
#' The case 'type' = "LAN" means that inference on the model parameter will be carried out under the local asymptotic normality (Le Cam and Yang, 2000) condition.
#' If the 'siblle' object has 'params' attribute whose length is equal to the length of the object, then 'type' defaults to "LAN".
#' If the 'siblle' object does not have 'params' attribute, then 'type' defaults to "point".
#'
#' When 'type' = "point", 'ci' can only be "loglik".
#' In this case a confidence interval for the log likelihood given the observed data is constructed.
#'
#' When 'type' = "regression", 'ci' can be "loglik", "MLE", or "information".
#' If 'ci' = "loglik", confidence intervals are constructed for the value of the log likelihood function evaluated at 'at.param' given the observed data.
#' If 'ci' = "MLE", confidence intervals are constructed for the location of the maximum likelihood estimate given the observed data.
#' If 'ci' = "information", confidence intervals are constructed for the Fisher information, which is assumed to be equal to (-2) times the value of the quadratic coefficient of the regression polynomial describing the mean of the SIBLLE.
#' When 'type' = "regression", 'ci' = "MLE" is assumed by default.
#'
#' When 'type' = "LAN", 'ci' can be "loglik", "MLE", "information", or "parameter".
#' If 'ci' is "loglik", "MLE", or "information", the output is the same as in the case where 'type' = "regression".
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
#' \item{quantiles_numerical_error_size: When 'ci'= "loglik", confidence intervals are found using the MLLR_2 distributions, whose quantile values are numerically found using random number generations. Thus these quantile values have some stochastic error. This 'quantiles_numerical_error_size' component gives the size of the numerical error, which is automatically set to not exceed approximately 0.01. Note that when 'ci'="MLE", "information", or "parameter", the (standard) F distribution is used, so this list component is omitted.}
#' }
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @references Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American statistical association, 74(368), 829-836.
#' @references Le Cam, L. and Yang, G. L. (2000). Asymptotics in statistics: some basic concepts. Springer-Verlag, New York.
#' @export
ci.siblle <- function(siblle, level, type=NULL, ci=NULL, at.param=NULL, weight.at.param=NULL, weights=NULL, fraction=NULL, center=NULL) {
    validate_siblle(siblle)
    if (!is.null(type)) {
        match.arg(type, c("point", "regression", "LAN"))
    }
    if (is.null(type) && is.null(attr(siblle, "params"))) {
        type <- "point"
    }
    if (is.null(type) && !is.null(attr(siblle, "params")) && length(siblle) == length(attr(siblle, "params"))) {
        type <- "LAN"
    }

    if (!is.null(ci)) {
        match.arg(ci, c("loglik", "MLE", "information", "parameter"))
    }
    if (is.null(ci)) {
        if (type=="point") {
            ci <- "loglik"
        }
        if (type=="regression") {
            ci <- "MLE"
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
    if (!is.null(at.param)) {
        if (!is.numeric(at.param)) {
            stop(
                "'at.param' should be a numeric value (or NULL).",
                call. = FALSE
            )
        }
        if (length(at.param)!=1) {
            stop(
                "'at.param' should be a single numeric value (or NULL).",
                call. = FALSE
            )
        }
    }
    if (!is.null(weight.at.param)) {
        if (!is.numeric(weight.at.param)) {
            stop(
                "'weight.at.param' should be a numeric value (or NULL).",
                call. = FALSE
            )
        }
        if (length(weight.at.param)!=1) {
            stop(
                "'weight.at.param' should be a single numeric value (or NULL).",
                call. = FALSE
            )
        }
    }
    if (type=="point" && ci != "loglik") {
        stop(
            "When 'type' = 'point', 'ci' should be 'loglik'.",
            call. = FALSE
        )
    }
    if (type=="regression" && !ci %in% c("loglik", "MLE", "information")) {
        stop(
            "When 'type' = 'regression', 'ci' should be one of 'loglik', 'MLE' (default), or 'information'.",
            call. = FALSE
        )
    }
    if (type=="LAN" && !ci %in% c("loglik", "MLE", "information", "parameter")) {
        stop(
            "When 'type' = 'LAN', 'ci' should be one of 'loglik', 'MLE', 'information', or 'parameter' (default).",
            call. = FALSE
        )
    }
    if (type %in% c("regression", "LAN") && ci=="loglik" && is.null(at.param)) {
        stop(
            "If 'type' is 'regression' or 'LAN' and 'ci' = 'loglik', then the point at which the log likelihood is to be estimated should be specified by 'at.param'.",
            call. = FALSE
        )
    }
    if (type!="LAN" && ci=="parameter") {
        stop(
            "'ci' = 'parameter' is only available when 'type' = 'LAN'.",
            call. = FALSE
        )
    }

    if (type=="point") {
        llest <- c(unclass(siblle))
        muhat <- mean(llest)
        Ssq <- var(llest)
        M <- length(llest)
        sighat <- sqrt((M-1)/M*Ssq)
        if (ci=="loglik") {
            num.error.size <- 0.01
            qmllr2out <- qmllr2(1-level, M, 1, 1/M, sighat, num_error_size=num.error.size) # level = 1-alpha
            if (length(qmllr2out)==0) { # execution of qmllr2 stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qmllr2out$quantiles
            num.error.size <- qmllr2out$numerical_error_size
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[i]
                f <- function(gamma) {
                    -gamma + log(gamma)
                } # a function that takes its maximum at 1
                ## Here gamma = (M-1)/M*Ssq/sigma_0^2. The interval (B) defined by {gamma; f(gamma) + 1 - 2/M*MLLR1_{1-alpha}(M,1,1/M,sighat) >= 0} gives the range for which the radicand in the expression for the confidence interval is non-negative.
                fprime <- function(gamma) { # derivative of gamma
                    -1 + 1/gamma
                }
                tol <- .Machine$double.eps^.5 # uniroot's default tolerance level for numerical root finding
                gamma_min <- uniroot(function(g) {f(g)+1-2/M*q[i]-tol}, interval=c(tol, 1))$root ## the lower limit of the interval B is between 0 and M
                gamma_max <- uniroot(function(g) {f(g)+1-2/M*q[i]-tol}, interval=c(1, 2+(2*q[i]/M-1-f(2))/fprime(2)))$root ## the upper limit of interval B is between M and 2+(2*q/M-1-f(2))/f'(2) (think of a tangent line at gamma=2)
                sigma0_sq_min <- (M-1)/M*Ssq/gamma_max # (numerically) smallest sigma0_sq such that the radicand is nonnegative
                sigma0_sq_max <- (M-1)/M*Ssq/gamma_min # (numerically) largest sigma0_sq such that the radicand is nonnegative
                lm <- function(x) { # lower margin (x = sigma_0^2)
                    x/2 - sqrt(x)*sqrt(-(M-1)/M*Ssq/x + log((M-1)/M*Ssq/x) + 1-2*q[i]/M)
                }
                um <- function(x) { # upper margin
                    x/2 + sqrt(x)*sqrt(-(M-1)/M*Ssq/x + log((M-1)/M*Ssq/x) + 1-2*q[i]/M)
                }
                lb <- muhat + optimize(lm, c(sigma0_sq_min, sigma0_sq_max), maximum=FALSE)$objective # lower bound
                ub <- muhat + optimize(um, c(sigma0_sq_min, sigma0_sq_max), maximum=TRUE)$objective # upper bound
                return(c(level=lvl, lb=lb, ub=ub))
            })
            out <- list(meta_model_MLE_for_log_lik=c(log_lik=muhat+(M-1)/(2*M)*Ssq),
                confidence_interval=t(lub),
                quantiles_numerical_error_size=num.error.size
            )
            return(out)
        }
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
        ## ci for log likelihood
        if (ci=="loglik") {
            tau <- c(c(1, at.param, at.param^2) %*% solve(t(theta012)%*%W%*%theta012, c(1, at.param, at.param^2)))
            if (is.null(weight.at.param)) {
                weight.at.param <- 1 # default value for w(at.param)
            }
            num.error.size <- 0.01 # target size of numerical error in evaluation of the quantiles of MLLR2 distribution
            qmllr2out <- qmllr2(1-level, M, 3, weight.at.param^2*tau, sqrt(sigsqhat), num_error_size=num.error.size)
            if (length(qmllr2out)==0) { # execution of qmllr1 stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qmllr2out$quantiles
            num.error.size <- qmllr2out$numerical_error_size
            mu.at <- sum(c(1,at.param,at.param^2)*Ahat) # estimated mean of SIBLLE at at.param
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[i]
                f <- function(gamma) {
                    -gamma + log(gamma)
                } # a function that takes its maximum at 1
                ## gamma = sigsqhat/sigma0^2. The interval (B) defined by {gamma; f(gamma) + 1 - 2/M*MLLR2{1-alpha}(M,3,w(theta)^2*tau,sighat)} gives the range in which the radicand in the expression for the confidence interval is non-negative.
                fprime <- function(gamma) { # derivative of gamma
                    -1 + 1/gamma
                }
                tol <- .Machine$double.eps^.5 # uniroot's default tolerance level for numerical root finding
                gamma_min <- uniroot(function(g) {f(g)+1-2/M*q[i]-tol}, interval=c(0, 1))$root ## the lower limit of the interval B is between 0 and 1
                gamma_max <- uniroot(function(g) {f(g)+1-2/M*q[i]-tol}, interval=c(1, 2+(2*q[i]/M-1-f(2))/fprime(2)))$root ## the upper limit of interval B is between 1 and 2+(2q/M-1-f(2))/f'(2) (think of a tangent line at gamma=2)
                sigma0_sq_min <- sigsqhat/gamma_max # (numerically) smallest sigma0_sq such that the radicand is nonnegative
                sigma0_sq_max <- sigsqhat/gamma_min # (numerically) largest sigma0_sq such that the radicand is nonnegative
                lm <- function(x) { # lower margin
                    x/(2*weight.at.param) - sqrt(x)*sqrt(tau*M*(log(sigsqhat/x) - sigsqhat/x + 1) - 2*tau*q[i])
                }
                um <- function(x) { # upper margin
                    x/(2*weight.at.param) + sqrt(x)*sqrt(tau*M*(log(sigsqhat/x) - sigsqhat/x + 1) - 2*tau*q[i])
                }
                lb <- mu.at + optimize(lm, c(sigma0_sq_min, sigma0_sq_max), maximum=FALSE)$objective # lower bound
                ub <- mu.at + optimize(um, c(sigma0_sq_min, sigma0_sq_max), maximum=TRUE)$objective # upper bound
                return(c(level=lvl, lb=lb, ub=ub))
            })
            out <- list(meta_model_MLE_for_log_lik=c(log_lik=unname(mu.at+sigsqhat/2)),
                confidence_interval=t(lub),
                MLLR2_quantiles_numerical_error_size=num.error.size
            )
            return(out)
        }
        ## ci for MLE
        if (ci=="MLE") {
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
                if (Disc >= 0 && coef2 < 0) { # Case 1
                    int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
                } else if (Disc >= 0 && coef2 >= 0) { # Case 2
                    int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
                } else { # Case 3
                    return(c(level=lvl, lb=-Inf, ub=Inf, inverted=0))
                }
            })
            if (any(lub["inverted",]==1)) { # if for any given level the confidence interval is inverted (Case 2)
                warning(paste0("For level(s) ", toString(unlist(level)[which(lub["inverted",]==1)]), ", the constructed confidence is of the form (-Inf, bound_1) U (bound_2, Inf)."), call.=FALSE)
            } else { # otherwise, remove the "inverted" column
                lub <- lub[-4,]
            }
            out <- list(meta_model_MLE_for_MLE=c(MLE=unname(-Ahat[2]/(2*Ahat[3]))),
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
                if (Disc >= 0 && coef2 < 0) { # Case 1
                    int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
                } else if (Disc >= 0 && coef2 >= 0) { # Case 2
                    int <- -1/(2*coef2)*(coef1+c(1,-1)*sqrt(Disc))
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
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
