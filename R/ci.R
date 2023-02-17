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
#' When quadratic regression is carried out, the weights for the Monte Carlo likelihood estimates can be supplied.  The weights can either be given as an attribute 'weights' of the 'mclle' object, or as a function argument 'weights', with the latter being used when both are supplied. In either case, 'weights' should be a numeric vector of length equal to that of 'mclle'. If 'weights' is given as a function argument, it can be specified alternatively as a character string "tricube". In this case, the tricube weight (see Cleveland, 1979) is used, and the specified 'fraction' of the points will have nonzero weights. If weights are not supplied in either locations, all weights are taken to be equal to 1.
#' It is important to note that the weights should NOT be normalized. Multiplying all weights by the same constant changes the local regression results. Roughly speaking, the variance of the error in the Monte Carlo log likelihood estimate is assumed to be sigma^2/(the weight for the point). See Park & Won (2023) for more information.
#'
#' @return A list consisting of the followings are returned.
#' \itemize{
#' \item{Monte Carlo maximum likelihood estimate,}
#' \item{a data frame of the lower and upper bounds of the confidence intervals and the corresponding (conservative) confidence levels,}
#' \item{the precision (0.01 or 0.001) for the (conservative) confidence levels.}
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
            } ## TODO: revise from here.
            sigmaxsq <- sapply(level, function(x) 2*(M*(muhat - x)^2+(M-1)*Ssq)/(M+sqrt(M*(M*(muhat-x)^2+(M-1)*Ssq)+M^2))) # the value of sigma0^2 that maximizes the LLR statistics
            teststats <- sapply(1:length(null.value), function(i) -.5*(muhat-null.value[[i]]+sigmaxsq[i]/2)^2/(sigmaxsq[i]/M) - (M-1)/2*Ssq/sigmaxsq[i] + M/2*log((M-1)*Ssq/(M*sigmaxsq[i])))
            prec <- 0.01
            pval <- pscl(teststats, M, 1, precision=prec)
            if (any(pval < .01)) {
                prec <- 0.001
                pval <- pscl(teststats, M, 1, precision=prec)
            }
            precdigits = ifelse(prec==0.01, 2, 3)
            dfout <- data.frame(
                log_lik_null=unlist(null.value),
                conservative_pvalue=round(pval, digits=precdigits)
            )
            out <- list(Monte_Carlo_MLE=c(log_lik=muhat+(M-1)/(2*M)*Ssq),
                Hypothesis_Tests=dfout,
                pvalue_precision=prec
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
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest))
        resids <- llest - c(theta012%*%Ahat)
        sig2hat <- c(resids%*%W%*%resids) / M
        ## test about moments
        if (ci=="moments") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            if (any(sapply(null.value, function(x) x[4]<=0))) {
                stop("The fourth component of null.value (the variance of Monte Carlo log likelihood estimator) should be positive.",
                    call. = FALSE
                )
            }
            teststats <- sapply(null.value,
                function(x) {
                    err <- llest - c(theta012%*%x[1:3])
                    .5*M*log(sig2hat/x[4]) - .5*c(err%*%W%*%err)/x[4]
                })
            prec <- 0.01
            pval <- pscl(teststats, M, 3, precision=prec)
            if (any(pval < .01)) {
                prec <- 0.001
                pval <- pscl(teststats, M, 3, precision=prec)
            }
            precdigits = ifelse(prec==0.01, 2, 3)
            dfout <- data.frame(
                a_null=sapply(null.value, function(x) x[1]),
                b_null=sapply(null.value, function(x) x[2]),
                c_null=sapply(null.value, function(x) x[3]),
                sigma_sq_null=sapply(null.value, function(x) x[4]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(Monte_Carlo_MLE=c(a=Ahat[1], b=Ahat[2], c=Ahat[3], sigma_sq=sig2hat),
                Hypothesis_Tests=dfout,
                pvalue_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## test about log likelihood
        if (ci=="loglik") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            nu <- c(c(1, param.at, param.at^2)%*%solve(t(theta012)%*%W%*%theta012, c(1, param.at, param.at^2)))
            mu.at <- sum(c(1,param.at,param.at^2)*Ahat) # estimated mean of MCLLE at param.at
            sigmaxsq <- sapply(null.value,
                function(x) {
                    nu1 <- M/2*sig2hat + (x-mu.at)^2/(2*nu)
                    nu2 <- M/2
                    nu3 <- 1/(8*nu)
                    1/(nu2/(2*nu1)+sqrt(nu3/nu1+nu2^2/(4*nu1^2)))
                }) # the value of sigma0^2 that maximizes the LLR statistics
            teststats <- sapply(1:length(null.value), function(i) M/2*log(sig2hat/sigmaxsq[i]) - M*sig2hat/(2*sigmaxsq[i]) - (null.value[[i]]-sigmaxsq[i]/2-mu.at)^2/(2*nu*sigmaxsq[i]))
            prec <- 0.01
            pval <- pscl(teststats, M, 3, precision=prec)
            if (any(pval < .01)) {
                prec <- 0.001
                pval <- pscl(teststats, M, 3, precision=prec)
            }
            precdigits = ifelse(prec==0.01, 2, 3)
            dfout <- data.frame(
                log_lik_null=unlist(null.value),
                conservative_pvalue=round(pval, digits=precdigits)
            )
            out <- list(Monte_Carlo_MLE=c(log_lik=unname(mu.at+sig2hat/2)),
                Hypothesis_Tests=dfout,
                pvalue_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## ci about MLE
        if (ci=="MLE") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            mtheta1 <- sum(w*theta)/sum(w)
            mtheta2 <- sum(w*theta*theta)/sum(w)
            mtheta3 <- sum(w*theta*theta*theta)/sum(w)
            mtheta4 <- sum(w*theta*theta*theta*theta)/sum(w)
            v11 <- sum(w)*(mtheta2 - mtheta1*mtheta1)
            v12 <- sum(w)*(mtheta3 - mtheta1*mtheta2)
            v22 <- sum(w)*(mtheta4 - mtheta2*mtheta2)
            teststats <- sapply(null.value,
                function(x) {
                    -M/2*log(1 + 1/(M*sig2hat)*(Ahat[2]+2*x*Ahat[3])^2/(v22-4*v12*x+4*v11*x*x)*(v11*v22-v12^2)) - M/2
                })
            prec <- 0.01
            pval <- pscl(teststats, M, 3, precision=prec)
            if (any(pval < .01)) {
                prec <- 0.001
                pval <- pscl(teststats, M, 3, precision=prec)
            }
            precdigits = ifelse(prec==0.01, 2, 3)
            dfout <- data.frame(
                MLE_null=unlist(null.value),
                conservative_pvalue=round(pval, digits=precdigits)
            )
            out <- list(Monte_Carlo_MLE=c(MLE=unname(-Ahat[2]/(2*Ahat[3]))),
                Hypothesis_Tests=dfout,
                pvalue_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## ci about Fisher information
        if (ci=="information") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            U <- t(theta012)%*%W%*%theta012
            u3gv12 <- U[3,3] - c(U[3,1:2]%*%solve(U[1:2,1:2], U[1:2,3])) # u_{3|12}
            teststats <- sapply(null.value,
                function(x) {
                    -M/2*log(1+ 1/(M*sig2hat)*(Ahat[3]+x/2)^2*u3gv12) - M/2
                })
            prec <- 0.01
            pval <- pscl(teststats, M, 3, precision=prec)
            if (any(pval < .01)) {
                prec <- 0.001
                pval <- pscl(teststats, M, 3, precision=prec)
            }
            precdigits = ifelse(prec==0.01, 2, 3)
            dfout <- data.frame(
                information_null=unlist(null.value),
                conservative_pvalue=round(pval, digits=precdigits)
            )
            out <- list(Monte_Carlo_MLE=c(Fisher_information=unname(-2*Ahat[3])),
                Hypothesis_Tests=dfout,
                pvalue_precision=prec
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
    }
    ## ci about the model parameter under LAN
    if (type=="LAN" && ci=="parameter") {
        warning("For parameter estimation under the LAN assumption, all Monte Carlo log likelihood estimates in the 'mclle' object are used with weights equal to 1. If the 'weights' argument is supplied to the 'ht' function or if the 'mclle' object has the 'weights' attribute, it is ignored.", call.=FALSE)
        theta <- attr(mclle, "param")
        llest <- c(unclass(mclle))
        M <- length(llest)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%theta012, t(theta012)%*%llest))
        resids_fs <- llest - c(theta012%*%Ahat)
        sig2hat_fs <- sum(resids_fs*resids_fs) / M # the first stage estimate of sigma^2
        if (!is.list(null.value)) {
            null.value <- list(null.value)
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
        cat("sig2_fs", sig2hat_fs, "K_fs", -2*Ahat[3], "sig2_ss", sig2hat_ss, "K_ss", K_ss, "\n")
        teststats <- sapply(null.value,
            function(x) {
                thdsq <- cbind(theta_chk, thetasq_chk)%*%c(x, -.5) # an intermediate step in computation
                iota <- c(llest_chk%*%G1%*%llest_chk - (t(llest_chk)%*%G1%*%thdsq)^2/(t(thdsq)%*%G1%*%thdsq))
                -(M-1)/2 + (M-1)/2*log((M-1)*sig2hat_ss/iota)
            })
        cat(teststats, "\n")
        prec <- 0.01
        pval <- pscl(teststats, M-1, 2, precision=prec)
        if (any(pval < .01)) {
            prec <- 0.001
            pval <- pscl(teststats, M-1, 2, precision=prec)
        }
        precdigits = ifelse(prec==0.01, 2, 3)
        dfout <- data.frame(
            parameter_null=unlist(null.value),
            conservative_pvalue=round(pval, digits=precdigits)
        )
        out <- list(Monte_Carlo_MLE=c(parameter=theta_ss, information=K_ss, error_variance=sig2hat_ss),
            Hypothesis_Tests=dfout,
            pvalue_precision=prec
        )
        print(out, row.names=FALSE)
        invisible(out)
    }
}


