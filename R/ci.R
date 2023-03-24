#' @export
ci <- function(x, ...) {
    UseMethod("ci")
}

#' Confidence intervals constructed using Monte Carlo log likelihood estimates in simulation-based inference
#'
#' `ci` constructs conservative confidence intervals using Monte Carlo log likelihood estimates. See Park and Won (2023) for more information.
#'
#' @name ci
#' @param mclle A class 'mclle' object, containing the Monte Carlo log likelihood estimates and (optional) the parameter values for which those estimates were obtained.
#' @param level The conservative confidence level. Either a numeric (vector of length one) or a list of numerics (for constructing multiple confidence intervals).
#' @param type A character string indicating what type of situation is considered. One of "point", "regression", or "LAN". See Details.
#' @param ci A character string indicating the quantity for which a confidence interval is to be constructed. One of "loglik", "MLE", "information", or "parameter". See Details.
#' @param param.at If 'ci' = "loglik", a confidence interval is constructed for the value of the log likelihood function evaluated at 'param.at' (for the cases 'type' = "regression" or "LAN".)
#' @param weights An optional argument for the (non-relative) weights of the Monte Carlo log likelihood estimates for regression. Either a numeric vector of length equal to the 'mclle' object, or a character string equal to "tricube".
#' @param fraction An optional argument indicating the fraction of points with nonzero weights for the case where 'weights' is specified as "tricube".
#' @param center An optional argument indicating the center of the local regression for the case where 'weights' is specified as "tricube".
#'
#' @details
#' This is a generic function, taking a class 'mclle' object as the first argument.
#' Confidence intervals are constructed under the assumption that the Monte Carlo likelihood estimator whose values are given in the 'mclle' object is (approximately) normally distributed.
#' 
#' When 'level' is a list, a conservative confidence interval is constructed for each level specified in the list.
#'
#' The constructed confidence intervals are conservative, in that all points for which there exists a null distribution that is not rejected at the specified confidence level will be included in the constructed interval. This complication arises because a point in the confidence interval does not fully specify a distribution (instead only defines a subspace of the parameter space.) See Park and Won (2023) for more detailed explanation.
#'
#' The 'type' argument should be one of "point", "regression", or "LAN".
#' The case 'type' = "point" means that the 'mclle' object contains Monte Carlo log likelihood estimates for a single, fixed parameter value.
#' The case 'type' = "regression" means that the 'mclle' object contains Monte Carlo log likelihood estimates evaluated at a range of parameter values, specified by the 'param' attribute of the 'mclle' object. A local quadratic regression for the estimated log likelihood values will be used for constructing confidence intervals, where the x-axis values are given by the 'param' values of the 'mclle' object.
#' The case 'type' = "LAN" means that inference on the model parameter will be carried out under the local asymptotic normality (Le Cam and Yang, 2000) condition.
#' If the 'mclle' object has 'param' attribute whose length is equal to the length of the object, then 'type' defaults to "LAN".
#' If the 'mclle' object does not have 'param' attribute, then 'type' defaults to "point".
#'
#' When 'type' = "point", 'ci' can only be "loglik".
#' In this case a conservative confidence interval for the log likelihood given the observed data is constructed.
#' When 'type' = "point", 'test' = "loglik" is assumed by default.
#'
#' When 'type' = "regression", 'ci' can be "loglik", "MLE", or "information".
#' If 'ci' = "loglik", confidence intervals are constructed for the value of the log likelihood function evaluated at 'param.at' given the observed data.
#' If 'ci' = "MLE", confidence intervals are constructed for the location of the maximum likelihood estimate given the observed data.
#' If 'ci' = "information", confidence intervals are constructed for the Fisher information, which is assumed to be equal to (-2) times the value of the quadratic coefficient of the regression polynomial describing the mean of the MCLLE.
#' When 'type' = "regression", 'ci' = "MLE" is assumed by default.
#'
#' When 'type' = "LAN", 'ci' can be "loglik", "MLE", "information", or "parameter".
#' If 'ci' is "loglik", "MLE", or "information", the output is the same as in the case where 'type' = "regression".
#' If 'ci' is "parameter", confidence intervals for the value of the model parameter are constructed under the local asymptotic normality assumption.
#' When 'type' = "LAN", 'ci' = "parameter" is assumed by default.
#'
#' When quadratic regression is carried out, the weights for the Monte Carlo likelihood estimates can be supplied.  The weights can either be given as an attribute 'weights' of the 'mclle' object, or as a function argument 'weights', with the latter being used when both are supplied. In either case, 'weights' should be a numeric vector of length equal to that of 'mclle'. If 'weights' is given as a function argument, it can be specified alternatively as a character string "tricube". In this case, the tricube weight (see Cleveland, 1979) is used, and the specified 'fraction' of the points will have nonzero weights. The 'center' argument determines at which parameter value the tricube weight takes the maximum. If weights are not supplied in either location, all weights are taken to be equal to 1.
#' It is important to note that the weights should NOT be normalized. Multiplying all weights by the same constant changes the local regression results. Roughly speaking, the variance of the error in the Monte Carlo log likelihood estimate is assumed to be sigma^2/(the weight for the point). See Park & Won (2023) for more information.
#'
#' @return A list consisting of the followings are returned.
#' \itemize{
#' \item{Monte Carlo maximum likelihood estimate,}
#' \item{a data frame of the lower and upper bounds of the confidence intervals and the corresponding (conservative) confidence levels,}
#' \item{the precision of the quantile values of the SCL distribution used in the construction of confidence levels (0.01 by default).}
#' }
#' 
#' @references Park, J. and Won, S. (2023). Simulation-based inference for partially observed, implicitly defined models
#' @references Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American statistical association, 74(368), 829-836.
#' @references Le Cam, L. and Yang, G. L. (2000). Asymptotics in statistics: some basic concepts. Springer-Verlag, New York.
#' @export
ci.mclle <- function(mclle, level, type=NULL, ci=NULL, param.at=NULL, weights=NULL, fraction=NULL, center=NULL) {
    validate_mclle(mclle)
    if (!is.null(type)) {
        match.arg(type, c("point", "regression", "LAN"))
    }
    if (is.null(type) && is.null(attr(mclle, "param"))) {
        type <- "point"
    }
    if (is.null(type) && !is.null(attr(mclle, "param")) && length(mclle) == length(attr(mclle, "param"))) {
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
        if (!is.list(level)) {
            stop(
                "'level' should be numeric or a list of numeric values",
                call. = FALSE
            )
        }
        if (is.list(level)) {
            if (!all(sapply(level, is.numeric))) {
                stop(
                    "If 'level' is a list, all of its components should be numeric.",
                    call. = FALSE
                )
            }
            if (!all(sapply(level, length)==1)) {
                stop(
                    "If 'level' is a list, all of its components should be a numeric vector of length 1.",
                    call. = FALSE
                )
            }
        }
    } else if (length(level) != 1) {
            stop(
                "If 'level' is a numeric vector (and not a list), its length should be 1.",
                call. = FALSE
            )
    }
    if (!is.null(param.at)) {
        if (!is.numeric(param.at)) {
            stop(
                "'param.at' should be a numeric value (or NULL).",
                call. = FALSE
            )
        }
        if (length(param.at)!=1) {
            stop(
                "'param.at' should be a single numeric value (or NULL).",
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
    if (type %in% c("regression", "LAN") && ci=="loglik" && is.null(param.at)) {
        stop(
            "If 'type' is 'regression' or 'LAN' and 'ci' = 'loglik', then the point at which the log likelihood is to be estimated should be specified by 'param.at'.",
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
        llest <- c(unclass(mclle))
        muhat <- mean(llest)
        Ssq <- var(llest)
        M <- length(llest)
        if (ci=="loglik") {
            if (!is.list(level)) {
                level <- list(level)
            }
            prec <- 0.01
            qsclout <- qscl(1-unlist(level), M, 1, precision=prec)
            if (length(qsclout)==0) { # execution of qscl stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qsclout$quantiles
            prec <- qsclout$precision
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[[i]]
                f <- function(gamma) {
                    -gamma/2 + M/2*log(gamma/M)
                } # a function that takes its maximum at M
                ## the interval (B) defined by {gamma; f(gamma) >= SCL_{1-alpha}(M,1)} is related to the interval for sigma_0^2 such that the radicand in the expression for the confidence interval is non-negative. The relationship is gamma = (M-1)*Ssq/sigma_0^2.
                fprime <- function(gamma) { # derivative of gamma
                    -1/2 + M/2/gamma
                } 
                tol <- .Machine$double.eps^.25 # uniroot's default tolerance level for numerical root finding
                gamma_min <- uniroot(function(g) {f(g)-q[i]-tol}, interval=c(M*exp(2/M*q[i]), M))$root ## the lower limit of the interval B is between M*exp(2/M*q) and M
                gamma_max <- uniroot(function(g) {f(g)-q[i]-tol}, interval=c(M, 2*M+(q[i]-f(2*M))/fprime(2*M)))$root ## the upper limit of interval B is between M and 2M+(q-f(2M))/f'(2M)
                sigma0_sq_min <- (M-1)*Ssq/gamma_max # (numerically) smallest sigma0_sq such that the radicand is nonnegative
                sigma0_sq_max <- (M-1)*Ssq/gamma_min # (numerically) largest sigma0_sq such that the radicand is nonnegative
                lm <- function(x) { # lower margin
                    x/2 - sqrt(2*x/M*(-(M-1)/2*Ssq/x + M/2*log((M-1)*Ssq/(M*x)) - q[i]))
                }
                um <- function(x) { # upper margin
                    x/2 + sqrt(2*x/M*(-(M-1)/2*Ssq/x + M/2*log((M-1)*Ssq/(M*x)) - q[i]))
                }
                lb <- muhat + optimize(lm, c(sigma0_sq_min, sigma0_sq_max), maximum=FALSE)$objective # lower bound
                ub <- muhat + optimize(um, c(sigma0_sq_min, sigma0_sq_max), maximum=TRUE)$objective # upper bound
                return(c(level=lvl, lb=lb, ub=ub))
            })
            out <- list(Monte_Carlo_MLE=c(log_lik=muhat+(M-1)/(2*M)*Ssq),
                conservative_confidence_interval=t(lub),
                confidence_level_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
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
            if (length(weights) != length(mclle) && weights!="tricube") {
                stop(
                    "When 'type' = 'regression' or 'LAN' and the 'weights' argument is given, the length of 'weights' should be the same as the length of 'mclle'.",
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
                distance <- abs(attr(mclle, "param") - center)
                span <- distance[order(distance)[ceiling(fraction*length(mclle))]]
                tricube <- function(x) { pmax((1-abs(x)^3)^3, 0) }
                w <- tricube(distance / span)
            }
            if (is.numeric(weights)) {
                w <- weights
            }
        }
        if (is.null(weights)) {
            if (!is.null(attr(mclle, "weights"))) {
                if (!is.numeric(attr(mclle, "weights"))) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and the 'mclle' object has 'weights' attribute, it has to be a numeric vector.",
                        call. = FALSE
                    )
                }
                if (length(mclle) != length(attr(mclle, "weights"))) {
                    stop(
                        "When 'type' = 'regression' or 'LAN' and the 'mclle' object has 'weights' attribute, the length of 'weights' should be the same as the length of 'mclle'.",
                        call. = FALSE
                    )
                }
                w <- attr(mclle, "weights")
            } else {
                w <- rep(1, length(mclle))
            }
        }
        ## weighted quadratic regression
        W  <- diag(w)
        theta <- attr(mclle, "param")
        llest <- c(unclass(mclle))
        M <- length(llest)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest)) # Ahat=(ahat,bhat,chat)
        resids <- llest - c(theta012%*%Ahat)
        sig2hat <- c(resids%*%W%*%resids) / M
        ## ci for log likelihood
        if (ci=="loglik") {
            if (!is.list(level)) {
                level <- list(level)
            }
            prec <- 0.01
            qsclout <- qscl(1-unlist(level), M, 3, precision=prec)
            if (length(qsclout)==0) { # execution of qscl stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qsclout$quantiles
            prec <- qsclout$precision
            nu <- c(c(1, param.at, param.at^2)%*%solve(t(theta012)%*%W%*%theta012, c(1, param.at, param.at^2)))
            mu.at <- sum(c(1,param.at,param.at^2)*Ahat) # estimated mean of MCLLE at param.at
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[[i]]
                f <- function(gamma) {
                    log(gamma) - gamma
                } # a function that takes its maximum at 1
                ## the interval (B) defined by {gamma; M*f(gamma) >= 2*SCL_{1-alpha}(M,3)} is related to the interval for sigma0_sq such that the radicand in the expression for the confidence interval is non-negative. The relationship is gamma = sig2hat/sigma0_sq
                fprime <- function(gamma) { # derivative of gamma
                    1/gamma - 1
                } 
                tol <- .Machine$double.eps^.25 # uniroot's default tolerance level for numerical root finding
                gamma_min <- uniroot(function(g) {f(g)-2*q[i]/M-tol}, interval=c(exp(2*q[i]/M), 1))$root ## the lower limit of the interval B is between exp(2q/M) and 1
                gamma_max <- uniroot(function(g) {f(g)-2*q[i]/M-tol}, interval=c(1, 2+(2*q[i]/M-f(2))/fprime(2)))$root ## the upper limit of interval B is between 1 and 2+(2q/M-f(2))/f'(2)
                sigma0_sq_min <- sig2hat/gamma_max # (numerically) smallest sigma0_sq such that the radicand is nonnegative
                sigma0_sq_max <- sig2hat/gamma_min # (numerically) largest sigma0_sq such that the radicand is nonnegative
                lm <- function(x) { # lower margin
                    x/2 - sqrt(x*nu*(M*log(sig2hat/x) - M*sig2hat/x - 2*q[i]))
                }
                um <- function(x) { # upper margin
                    x/2 + sqrt(x*nu*(M*log(sig2hat/x) - M*sig2hat/x - 2*q[i]))
                }
                lb <- mu.at + optimize(lm, c(sigma0_sq_min, sigma0_sq_max), maximum=FALSE)$objective # lower bound
                ub <- mu.at + optimize(um, c(sigma0_sq_min, sigma0_sq_max), maximum=TRUE)$objective # upper bound
                return(c(level=lvl, lb=lb, ub=ub))
            })
            out <- list(Monte_Carlo_MLE=c(log_lik=unname(mu.at+sig2hat/2)),
                conservative_confidence_interval=t(lub),
                quantile_SCL_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## ci for MLE
        if (ci=="MLE") {
            if (!is.list(level)) {
                level <- list(level)
            }
            mtheta1 <- sum(w*theta)/sum(w)
            mtheta2 <- sum(w*theta*theta)/sum(w)
            mtheta3 <- sum(w*theta*theta*theta)/sum(w)
            mtheta4 <- sum(w*theta*theta*theta*theta)/sum(w)
            v11 <- sum(w)*(mtheta2 - mtheta1*mtheta1)
            v12 <- sum(w)*(mtheta3 - mtheta1*mtheta2)
            v22 <- sum(w)*(mtheta4 - mtheta2*mtheta2)
            detV <- v11*v22-v12*v12
            Vl <- sum(w*llest^2) - sum(w*llest)^2/sum(w)
            prec <- 0.01
            qsclout <- qscl(1-unlist(level), M, 3, precision=prec)
            if (length(qsclout)==0) { # execution of qscl stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qsclout$quantiles
            prec <- qsclout$precision
            xi <- M*(exp(-1-2*q/M)-1)
            D <- detV * xi * sig2hat * (Vl - (xi+M)*sig2hat)
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[[i]]
                if (xi[i]*sig2hat < Ahat[3]^2*detV/v11) { # Case 1
                    int <- (-(Ahat[2]*Ahat[3]*detV+v12*xi[i]*sig2hat) + c(-1,1)*sqrt(D[i])) /2/(Ahat[3]^2*detV-v11*xi[i]*sig2hat)
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
                } else if (Ahat[3]^2*detV/v11 <= xi[i]*sig2hat && xi[i]*sig2hat < Vl-M*sig2hat) { # Case 2
                    int <- (-(Ahat[2]*Ahat[3]*detV+v12*xi[i]*sig2hat) + c(1,-1)*sqrt(D[i])) /2/(Ahat[3]^2*detV-v11*xi[i]*sig2hat)
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=1))
                } else { # Case 3
                    return(c(level=lvl, lb=-Inf, ub=Inf, inverted=0))
                }
            })
            if (any(lub["inverted",]==1)) { # if for any given level the confidence interval is inverted (Case 2)
                warning(paste0("For level(s) ", toString(unlist(level)[which(lub["inverted",]==1)]), ", the constructed confidence is of the form (-Inf, lb) U (ub, Inf)."), call.=FALSE)
            } else { # otherwise, remove the "inverted" column
                lub <- lub[-4,]
            }
            out <- list(Monte_Carlo_MLE=c(MLE=unname(-Ahat[2]/(2*Ahat[3]))),
                conservative_confidence_interval=t(lub),
                quantile_SCL_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## ci for Fisher information
        if (ci=="information") {
            if (!is.list(level)) {
                level <- list(level)
            }
            U <- t(theta012)%*%W%*%theta012
            u3gv12 <- U[3,3] - c(U[3,1:2]%*%solve(U[1:2,1:2], U[1:2,3])) # u_{3|12}
            prec <- 0.01
            qsclout <- qscl(1-unlist(level), M, 3, precision=prec)
            if (length(qsclout)==0) { # execution of qscl stopped by user input
                stop("Construction of confidence interval stopped by user input", call. = FALSE)
            }
            q <- qsclout$quantiles
            prec <- qsclout$precision
            xi <- M*(exp(-1-2*q/M)-1)
            lub <- sapply(1:length(level), function(i) {
                lvl <- level[[i]]
                int <- -2*Ahat[3] + c(-1,1)*2*sqrt(xi[i]/u3gv12*sig2hat)
                return(c(level=lvl, lb=max(0,int[1]), ub=max(0,int[2])))
            })
            out <- list(Monte_Carlo_MLE=c(Fisher_information=unname(-2*Ahat[3])),
                conservative_confidence_interval=t(lub),
                quantile_SCL_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
    }
    ## ci for the model parameter under LAN
    if (type=="LAN" && ci=="parameter") {
        warning("For parameter estimation under the LAN assumption, all Monte Carlo log likelihood estimates in the 'mclle' object are used with weights equal to 1. If the 'weights' argument is supplied to the 'ht' function or if the 'mclle' object has the 'weights' attribute, it is ignored.", call.=FALSE)
        theta <- attr(mclle, "param")
        llest <- c(unclass(mclle))
        M <- length(llest)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%theta012, t(theta012)%*%llest))
        resids_fs <- llest - c(theta012%*%Ahat)
        sig2hat_fs <- sum(resids_fs*resids_fs) / M # the first stage estimate of sigma^2
        if (!is.list(level)) {
            level <- list(level)
        }
        theta_chk <- theta - mean(theta)
        thetasq_chk <- theta^2 - mean(theta^2)
        llest_chk <- llest - mean(llest)
        G1 <- diag(rep(1,M)) - outer(theta_chk, theta_chk)/(sum(theta_chk*theta_chk)-sig2hat_fs/(2*Ahat[3])) # G_{(1)} in the paper
        est_ss <- solve(t(cbind(theta_chk, thetasq_chk))%*%G1%*%cbind(theta_chk, thetasq_chk), t(cbind(theta_chk, thetasq_chk))) %*% G1 %*% llest_chk # second-stage estimates vector, to be used in the next two lines
        K_ss <- unname(-2*est_ss[2,1])
        theta_ss <- unname(est_ss[1,1]/K_ss)
        resids_ss <- unname(llest_chk - cbind(theta_chk, thetasq_chk)%*%est_ss)
        sig2hat_ss <- 1/(M-1)*c(t(resids_ss)%*%G1%*%resids_ss)
        prec <- 0.01
        qsclout <- qscl(1-unlist(level), M-1, 2, precision=prec)
        if (length(qsclout)==0) { # execution of qscl stopped by user input
            stop("Construction of confidence interval stopped by user input", call. = FALSE)
        }
        q <- qsclout$quantiles
        prec <- qsclout$precision
        xi <- c(llest_chk%*%G1%*%llest_chk) - (M-1)*sig2hat_ss*exp(-1-2*q/(M-1))
        lub <- sapply(1:length(level), function(i) {
            lvl <- level[[i]]
            qco <- c((llest_chk%*%G1%*%theta_chk)^2 - xi[i]*theta_chk%*%G1%*%theta_chk) # quadratic term coefficient for the quadratic polynomial that determines the Monte Carlo CI
            lco <- -c((llest_chk%*%G1%*%thetasq_chk)*(llest_chk%*%G1%*%theta_chk) - xi[i]*theta_chk%*%G1%*%thetasq_chk) # linear term coefficient
            con <- 1/4*c((llest_chk%*%G1%*%thetasq_chk)^2 - xi[i]*thetasq_chk%*%G1%*%thetasq_chk)
            D <- lco^2 - 4*qco*con
            if (D > 0) {
                int <- 1/(2*qco)*(-lco+c(1,-1)*sqrt(D))
                if (qco < 0) {
                    return(c(level=lvl, lb=int[1], ub=int[2], inverted=0))
                } else {
                    return(c(level=lvl, lb=int[2], ub=int[1], inverted=1))
                }
            } else {
                if (qco > 0) {
                    return(c(level=lvl, lb=-Inf, ub=Inf, inverted=0))
                } else {
                    return(c(level=lvl, lb=NA, ub=NA, inverted=0)) 
                }
            }
        })
        if (any(lub["inverted",]==1)) { # if for any given level the confidence interval is inverted 
            warning(paste0("For level(s) ", toString(unlist(level)[which(lub["inverted",]==1)]), ", the constructed confidence is of the form (-Inf, lb) U (ub, Inf)."), call.=FALSE)
        } else { # otherwise, remove the "inverted" column
            lub <- lub[-4,]
        }
        out <- list(Monte_Carlo_MLE=c(parameter=theta_ss),
            conservative_confidence_interval=t(lub),
            quantile_SCL_precision=prec
        )
        print(out, row.names=FALSE)
        invisible(out)
    }
}


