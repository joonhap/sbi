#' @export
ht <- function(x, ...) {
    UseMethod("ht")
}

#' Hypothesis tests using simulation based log likelihood estimates
#'
#' `ht` outputs results of hypothesis tests carried out using simulation based log likelihood estimates. See Park (2023) for more information.
#'
#' @name ht
#' @param siblle A class 'siblle' object, containing the simulation based log likelihood estimates and, if relevant, the parameter values for which those estimates were obtained.
#' @param null.value The null value for the hypothesis test. Either a numeric vector (for running a single test) or a list of numeric vectors (for running multiple tests).
#' @param type A character string indicating what type of situation is considered. One of "point", "regression", or "LAN". See Details.
#' @param test A character string indicating the quantity to be tested about. One of "loglik", "moments", "MLE", "information", or "parameter". See Details.
#' @param param.at For the cases 'type' = "regression" or "LAN" and 'test' = "loglik", the hypothesis test is about the value of the log likelihood function evaluated at 'param.at' 
#' @param weight.param.at The relative inverse variance for the simulation based log likelihood estimator at 'param.at'. Used for the cases 'type' = "regression" or "LAN" and 'test' = "loglik". The default value is 1.
#' @param weights An optional argument. The un-normalized weights of the log likelihood estimates for regression. Either a numeric vector of length equal to the 'siblle' object, or a character string equal to "tricube". The default weights are equal to one for all the points if not specified here or in the siblle object. See Details below.
#' @param fraction An optional argument used when the 'weights' argument is equal to "tricube". This argument specifies the fraction of points with nonzero weights when the tricube function is used for weight assignment.
#' @param center An optional argument indicating the center of the local regression for the case where 'weights' is specified as "tricube".
#'
#' @details
#' This is a generic function, taking a class 'siblle' object as the first argument.
#' Hypothesis tests are carried out under a normal meta model--that is, the log likelihood estimates (whose values are given in the 'siblle' object) are normally distributed.
#'
#' When 'null.value' is a list, a hypothesis test is carried out for each null value specified in the list.
#'
#' Some tests are exact under the normal meta model (e.g., tests on the mean and the variance of the log likelihood estimator) while others are approximate.  See Park (2023) for more information.
#'
#' The 'type' argument should be one of "point", "regression", or "LAN".
#' The case 'type' = "point" means that the 'siblle' object contains simulation based log likelihood estimates for a single, fixed parameter value.
#' The case 'type' = "regression" means that the 'siblle' object contains simulation based log likelihood estimates obtained at more than one parameter values, specified by the 'param' attribute of the 'siblle' object. A local quadratic regression for the estimated log likelihood values will be used for hypothesis tests, where the x-axis values are given by the 'param' values of the 'siblle' object.
#' The case 'type' = "LAN" means that inference on the model parameter will be carried out under the local asymptotic normality (Le Cam and Yang, 2000) condition.
#' If the 'siblle' object has 'param' attribute whose length is equal to the length of the object, then 'type' defaults to "LAN".
#' If the 'siblle' object does not have 'param' attribute, then 'type' defaults to "point".
#'
#' When 'type' = "point", 'test' can only be "loglik" or "moments".
#' In this case 'test' = "loglik" means the hypothesis test \eqn{H_0: l = null.value} versus \eqn{H_1: l != null.value} will be performed.
#' If 'test' = "moments", a test about the mean and the variance of the simulation based log likelihood estimator is conducted.
#' The 'null.value' should be a numeric vector of length two (the first component being the mean and the second being the variance), or a list of numeric vectors of length two.
#' When 'type' = "point", 'test' = "loglik" is assumed by default, unless the 'test' argument is supplied.
#'
#' When 'type' = "regression", 'test' can be "loglik", "moments", "MLE", or "information".
#' If 'test' = "loglik", the test is about the value of the log likelihood function evaluated at 'param.at'.
#' If 'test' = "moments", the test is about the quadruple \eqn{a, b, c, sigma^2} where \eqn{a, b, c} are coefficients of the polynomial describing the mean of the simulation based likelihood estimator (i.e., \eqn{l(\theta) = a + b \theta + c \theta^2}) and \eqn{\sigma^2} is the variance of the SIBLLE.
#' If 'test' = "MLE", the test is about the location of the maximum likelihood estimate.
#' If 'test' = "information", the test is about the Fisher information, which is (-2) times the value of \eqn{c}, the quadratic coefficient of the mean function of the SIBLLE.
#' When 'type' = "regression", 'test' = "MLE" is assumed by default.
#'
#' When 'type' = "LAN", 'test' can be "loglik", "moments", "MLE", "information", or "parameter".
#' If 'test' is "loglik", "moments", "MLE", or "information", the output is the same as in the case where 'type' = "regression".
#' If 'test' is "parameter", a test about the value of the model parameter is conducted under the local asymptotic normality assumption.
#' When 'type' = "LAN", 'test' = "parameter" is assumed by default.
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified.  The weights can either be given as an attribute 'weights' of the 'siblle' object, or as a function argument 'weights', with the latter being used when both are supplied. In either case, 'weights' should be a numeric vector of length equal to that of 'siblle'. If 'weights' is given as an argument to the "ht" function, it can be specified alternatively as a character string "tricube". In this case, the tricube weight (see Cleveland, 1979) is used, and the specified 'fraction' of the points will have nonzero weights. The 'center' argument determines at which parameter value the tricube weight takes the maximum. If weights are not supplied in either location, all weights are taken to be equal to 1.
#' It is important to note that the weights are un-normalized. Multiplying all weights by the same constant changes the local regression results. Roughly speaking, the variance of the simulation based log likelihood estimator is assumed to be sigma^2/(the weight for the point). See Park (2023) for more information.
#'
#' @return A list consisting of the following components are returned.
#' \itemize{
#' \item{meta model maximum likelihood estimate,}
#' \item{a data frame of the null values and the corresponding (approximate) p-values,}
#' \item{approximate size of error in numerical evaluation of p-values (0.01 or 0.001 or 0). When 'test'="moments" or "loglik", the MLLR_1 or MLLR_2 distributions are used to compute the p-values. The cumulative distribution functions for these distributions are numerically evaluated using random number generations, thus having some stochastic error. The size of the numerical error is automatically set to either approximately 0.01 or 0.001 for the output p-values. When 'test'="MLE", "information", or "parameter", the standard F distribution is used, so the numerical error is basically zero.}
#' }
#' When 'test' = "moments", exact p-values are shown (here "exact" means that the formula for the p-value is not based on approximation; this does not mean that size of the numerical evaluation is equal to zero.)
#' In other cases, approximate p-values are shown. See Park, J. (2023) for how approximations are made.
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @references Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American statistical association, 74(368), 829-836.
#' @references Le Cam, L. and Yang, G. L. (2000). Asymptotics in statistics: some basic concepts. Springer-Verlag, New York.
#' @export
ht.siblle <- function(siblle, null.value, type=NULL, test=NULL, param.at=NULL, weight.param.at=NULL, weights=NULL, fraction=NULL, center=NULL) {
    validate_siblle(siblle)
    if (!is.null(type)) {
        match.arg(type, c("point", "regression", "LAN"))
    }
    if (is.null(type) && is.null(attr(siblle, "param"))) {
        type <- "point"
    }
    if (is.null(type) && !is.null(attr(siblle, "param")) && length(siblle) == length(attr(siblle, "param"))) {
        type <- "LAN"
    }

    if (!is.null(test)) {
        match.arg(test, c("loglik", "moments", "MLE", "information", "parameter"))
    }
    if (is.null(test)) {
        if (type=="point") {
            test <- "loglik"
        }
        if (type=="regression") {
            test <- "MLE"
        }
        if (type=="LAN") {
            test <- "parameter"
        }
    }

    if (!is.numeric(null.value)) {
        if (!is.list(null.value)) {
            stop(
                "'null.value' should be numeric or a list of numeric values",
                call. = FALSE
            )
        }
        if (is.list(null.value) && !all(sapply(null.value, is.numeric))) {
            stop(
                "If 'null.value' is a list, all of its components should be numeric.",
                call. = FALSE
            )
        }
    }

    if (is.list(null.value)) {
        if (test=="loglik" && !all(sapply(null.value, length)==1)) {
            stop(
                "If 'null.value' is a list and 'test' is 'loglik', all components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
        if (test=="moments") {
            if (type=="point" && !all(sapply(null.value, length)==2)) {
                stop(
                    "If 'null.value' is a list, 'test' is 'moments', and 'type' is 'point', all components of 'null.value' should be a numeric vector of length 2 (mean and variance of SIBLLE).",
                    call. = FALSE
                )
            }
            if (type %in% c("regression", "LAN") && !all(sapply(null.value, length)==4)) {
                stop(
                    "If 'null.value' is a list, 'test' is 'moments', and 'type' is 'regression' or 'LAN', all components of 'null.value' should be a numeric vector of length 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the SIBLLE).",
                    call. = FALSE
                )
            }
        }
        if (test=="MLE" && !all(sapply(null.value, length)==1)) {
         stop(
                "If 'null.value' is a list and 'test' is 'MLE', all components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
        if (test=="information" && !all(sapply(null.value, length)==1)) {
         stop(
                "If 'null.value' is a list and 'test' is 'information', all components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
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
    if (!is.null(weight.param.at)) {
        if (!is.numeric(weight.param.at)) {
            stop(
                "'weight.param.at' should be a numeric value (or NULL).",
                call. = FALSE
            )
        }
        if (length(param.at)!=1) {
            stop(
                "'weight.param.at' should be a single numeric value (or NULL).",
                call. = FALSE
            )
        }
    }
    if (type=="point" && !test %in% c("loglik", "moments")) {
        stop(
            "When 'type' = 'point', 'test' should be either 'loglik' (default) or 'moments'.",
            call. = FALSE
        )
    }
    if (type=="regression" && !test %in% c("loglik", "moments", "MLE", "information")) {
        stop(
            "When 'type' = 'regression', 'test' should be one of 'loglik', 'moments', 'MLE' (default), or 'information'.",
            call. = FALSE
        )
    }
    if (type=="LAN" && !test %in% c("loglik", "moments", "MLE", "information", "parameter")) {
        stop(
            "When 'type' = 'LAN', 'test' should be one of 'loglik', 'moments', 'MLE', 'information', or 'parameter' (default).",
            call. = FALSE
        )
    }
    if (type %in% c("regression", "LAN") && test=="loglik" && is.null(param.at)) {
        stop(
            "If 'type' is 'regression' or 'LAN' and 'test' = 'loglik', then the point at which the log likelihood is to be estimated should be specified by 'param.at'.",
            call. = FALSE
        )
    }
    if (type!="LAN" && test=="parameter") {
        stop(
            "'test' = 'parameter' is only available when 'type' = 'LAN'.",
            call. = FALSE
        )
    }

    if (type=="point") {
        llest <- c(unclass(siblle))
        muhat <- mean(llest)
        Ssq <- var(llest)
        M <- length(llest)
        if (test=="moments") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            if (any(sapply(null.value, function(x) x[2]<=0))) {
                stop("The second component of null.value (the variance of simulation based log likelihood estimator) should be positive.",
                    call. = FALSE
                )
            }
            teststats <- sapply(null.value, function(x) -.5*M*(muhat - x[1])^2/x[2] - (M-1)/2*Ssq/x[2] + M/2*log((M-1)*Ssq/(M*x[2])) + M/2)
            num.error.size <- 0.01
            pvalout <- pmllr1(teststats, M, 1, precision=num.error.size)
            if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$precision
            if (any(pval < .01)) {
                num.error.size <- 0.001
                pvalout <- pmllr1(teststats, M, 1, precision=num.error.size)
                if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                prec <- pvalout$precision
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                mu_null=sapply(null.value, function(x) x[1]),
                sigma_sq_null=sapply(null.value, function(x) x[2]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE=c(mu=muhat, sigma_sq=(M-1)/M*Ssq),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        if (test=="loglik") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            sigmaxsq <- sapply(null.value, function(x) 2*(sqrt( (muhat-x)^2 + (M-1)/M*Ssq + 1 ) - 1)) # the value of sigma0^2 that maximizes the MLLR statistic
            teststats <- sapply(1:length(null.value), function(i) -.5*M*(muhat-null.value[[i]]+sigmaxsq[i]/2)^2/sigmaxsq[i] - (M-1)/2*Ssq/sigmaxsq[i] + M/2*log((M-1)*Ssq/(M*sigmaxsq[i])) + M/2)
            sighat <- sqrt((M-1)/M*Ssq)
            num.error.size <- 0.01
            pvalout <- pmllr2(teststats, M, 1, 1/M, sighat, precision=num.error.size)
            if (length(pvalout)==0) { # execution of pmllr2 stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$precision
            if (any(pval < .01)) {
                num.error.size <- 0.001
                pvalout <- pmllr2(teststats, M, 1, 1/M, sighat, precision=num.error.size)
                if (length(pvalout)==0) { # execution of pmllr2 stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$precision
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                log_lik_null=unlist(null.value),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE=c(log_lik=muhat+(M-1)/(2*M)*Ssq),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
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
                distance <- abs(attr(siblle, "param") - center)
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
        theta <- attr(siblle, "param")
        llest <- c(unclass(siblle))
        M <- length(llest)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest))
        resids <- llest - c(theta012%*%Ahat)
        sig2hat <- c(resids%*%W%*%resids) / M
        ## test about moments
        if (test=="moments") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            if (any(sapply(null.value, function(x) x[4]<=0))) {
                stop("The fourth component of null.value (the variance of simulation based log likelihood estimator) should be positive.",
                    call. = FALSE
                )
            }
            teststats <- sapply(null.value,
                function(x) {
                    err <- llest - c(theta012%*%x[1:3])
                    .5*M*log(sig2hat/x[4]) - .5*c(err%*%W%*%err)/x[4] + M/2
                })
            num.error.size <- 0.01
            pvalout <- pmllr1(teststats, M, 3, precision=num.error.size)
            if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$precision
            if (any(pval < .01)) {
                prec <- 0.001
                pvalout <- pmllr1(teststats, M, 3, precision=num.error.size)
                if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$precision
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                a_null=sapply(null.value, function(x) x[1]),
                b_null=sapply(null.value, function(x) x[2]),
                c_null=sapply(null.value, function(x) x[3]),
                sigma_sq_null=sapply(null.value, function(x) x[4]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE=c(a=Ahat[1], b=Ahat[2], c=Ahat[3], sigma_sq=sig2hat),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## test about log likelihood
        if (test=="loglik") {
            if (is.null(weight.param.at)) {
                weight.param.at <- 1 # default value for w(param.at)
            }
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            tau <- c(c(1, param.at, param.at^2)%*%solve(t(theta012)%*%W%*%theta012, c(1, param.at, param.at^2)))
            mu.param.at <- sum(c(1,param.at,param.at^2)*Ahat) # estimated mean of SIBLLE at param.at
            sigmaxsq <- sapply(null.value, function(x) 2*weight.param.at^2*tau*M*(sqrt( 1/weight.param.at^2/tau/M^2 * (1/tau*(x-mu.param.at)^2+M*sig2hat) + 1) - 1)) # the value of sigma0^2 that maximizes the LLR statistics
            teststats <- sapply(1:length(null.value), function(i) M/2*log(sig2hat/sigmaxsq[i]) - M*sig2hat/(2*sigmaxsq[i]) - (null.value[[i]]-sigmaxsq[i]/2/weight.param.at-mu.param.at)^2/(2*tau*sigmaxsq[i]) + M/2)
            num.error.size <- 0.01
            pvalout <- pmllr2(teststats, M, 3, weight.param.at^2*tau, sqrt(sig2hat), precision=num.error.size)
            if (length(pvalout)==0) { # execution of pmllr2 stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$precision
            if (any(pval < .01)) {
                num.error.size <- 0.001
                pvalout <- pmllr2(teststats, M, 3, weight.param.at^2*tau, sqrt(sig2hat), precision=num.error.size)
                if (length(pvalout)==0) { # execution of pmllr2 stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$precision
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                log_lik_null=unlist(null.value),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE=c(log_lik=unname(mu.param.at+sig2hat/2)),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## test about MLE
        if (test=="MLE") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            mtheta1 <- sum(w*theta)/sum(w)
            mtheta2 <- sum(w*theta*theta)/sum(w)
            mtheta3 <- sum(w*theta*theta*theta)/sum(w)
            mtheta4 <- sum(w*theta*theta*theta*theta)/sum(w)
            v11 <- sum(w)*(mtheta2 - mtheta1^2)
            v12 <- sum(w)*(mtheta3 - mtheta1*mtheta2)
            v22 <- sum(w)*(mtheta4 - mtheta2^2)
            teststats <- sapply(null.value,
                function(x) {
                    (M-3)*(Ahat[2]+2*x*Ahat[3])^2*(v11*v22-v12^2)/(v22-4*v12*x+4*v11*x*x)/(M*sig2hat)
                })
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                MLE_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE=c(MLE=unname(-Ahat[2]/(2*Ahat[3]))),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=0
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## test about Fisher information
        if (test=="information") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            U <- t(theta012)%*%W%*%theta012
            u3gv12 <- U[3,3] - c(U[3,1:2]%*%solve(U[1:2,1:2], U[1:2,3])) # u_{3|12}
            teststats <- sapply(null.value, function(x) (M-3)*u3gv12*(Ahat[3]+x/2)^2/(M*sig2hat))
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                information_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE=c(Fisher_information=unname(-2*Ahat[3])),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=0
            )
            print(out, row.names=FALSE)
            invisible(out)
        }
        ## test about the model parameter under LAN
        ## TODO: start revising from here
        if (type=="LAN" && test=="parameter") {
            W  <- diag(w)
            theta <- attr(siblle, "param")
            llest <- c(unclass(siblle))
            M <- length(llest)
            theta012 <- cbind(1, theta, theta^2)
            Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest))
            resids_1 <- llest - c(theta012%*%Ahat) # first stage estimates for residuals
            sig2hat_1 <- sum(resids_fs*resids_fs) / M # the first stage estimate of sigma^2
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
            teststats <- sapply(null.value,
                function(x) {
                    thdsq <- cbind(theta_chk, thetasq_chk)%*%c(x, -.5) # an intermediate step in computation
                    iota <- c(llest_chk%*%G1%*%llest_chk - (t(llest_chk)%*%G1%*%thdsq)^2/(t(thdsq)%*%G1%*%thdsq))
                    -(M-1)/2 + (M-1)/2*log((M-1)*sig2hat_ss/iota)
                })
            prec <- 0.01
            pvalout <- pmllr1(teststats, M-1, 2, precision=prec)
            if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            prec <- pvalout$precision
            if (any(pval < .01)) {
                prec <- 0.001
                pvalout <- pmllr1(teststats, M-1, 2, precision=prec)
                if (length(pvalout)==0) { # execution of pmllr1 stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                prec <- pvalout$precision
            }
            precdigits <- max(-floor(log10(prec)), 1)
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
}

