#' @export
ht <- function(x, ...) {
    UseMethod("ht")
}

#' Hypothesis tests using simulation based log likelihood estimates
#'
#' `ht` outputs results of hypothesis tests carried out using simulation based log likelihood estimates. See Park (2023) for more information.
#'
#' @name ht
#' @param siblle A class 'siblle' object, containing the simulation based log likelihood estimates and, if those estimates were obtained for different parameter values, the parameter values at which the log likelihood estimates were obtained. See help(siblle).
#' @param null.value The null value(s) for the hypothesis test. Either a numeric vector (for running a test for a single null value, which can have one or more components) or a list of numeric vectors (for running tests for multiple null values).
#' @param type A character string indicating what type of situation is considered. One of "point", "regression", or "LAN". See Details.
#' @param test A character string indicating the quantity to be tested about. One of "moments", "MLE", "information", or "parameter". See Details.
#' @param weights An optional argument. The un-normalized weights of the log likelihood estimates for regression. Either a numeric vector of length equal to the 'siblle' object, or a character string equal to "tricube". The default weights are equal to one for all the points if not specified here or in the siblle object. See Details below.
#' @param fraction An optional argument used when the 'weights' argument is equal to "tricube". This argument specifies the fraction of points with nonzero weights when the tricube function is used for weight assignment.
#' @param center An optional argument indicating the center of the local regression for the case where 'weights' is specified as "tricube".
#'
#' @details
#' This is a generic function, taking a class 'siblle' object as the first argument.
#' Hypothesis tests are carried out under a normal meta model--that is, the log likelihood estimates (whose values are given in the 'siblle' object) are normally distributed.
#'
#' When 'null.value' is a list, a hypothesis test is carried out for each null value specified in the list. For example, in order to run tests for more than one null values for a single-component parameter, you can let 'null.value=as.list(vector_of_null_values)'.
#'
#' Some tests are exact under the normal meta model (e.g., tests on the mean and the variance of the log likelihood estimator) while others are approximate.  See Park (2023) for more information.
#'
#' The 'type' argument should be one of "point", "regression", or "LAN".
#' The case 'type' = "point" means that the 'siblle' object contains simulation based log likelihood estimates for a single, fixed parameter value.
#' The case 'type' = "regression" means that the 'siblle' object contains simulation based log likelihood estimates obtained at more than one parameter values, specified by the 'params' attribute of the 'siblle' object. A local quadratic regression for the estimated log likelihood values will be used for hypothesis tests, where the x-axis values are given by the 'params' values of the 'siblle' object and the y-axis values are the corresponding log likelihood estimates.
#' The case 'type' = "LAN" means that inference on the model parameter will be carried out under the local asymptotic normality (Le Cam and Yang, 2000) condition.
#' If the 'siblle' object has 'params' attribute whose length is equal to the length of the object, then 'type' defaults to "LAN".
#' If the 'siblle' object does not have 'params' attribute, then 'type' defaults to "point".
#'
#' When 'type' = "point", 'test' can only be "moments".
#' If 'test' = "moments", a test about the mean and the variance of the simulation based log likelihood estimator is conducted.
#' The 'null.value' should be a numeric vector of length two (the first component being the mean and the second being the variance), or a list of numeric vectors of length two.
#' When 'type' = "point", 'test' = "moments" is assumed by default.
#'
#' When 'type' = "regression", 'test' can be "moments", "MLE", or "information".
#' If 'test' = "moments", the test is about the quadruple \eqn{a, b, c, sigma^2} where \eqn{a, b, c} are coefficients of the polynomial describing the mean of the simulation based likelihood estimator (i.e., \eqn{l(\theta) = a + b \theta + c \theta^2}) and \eqn{\sigma^2} is the variance of the SIBLLE.
#' If 'test' = "MLE", the test is about the location of the maximum likelihood estimate.
#' If 'test' = "information", the test is about the Fisher information, which is (-2) times the value of \eqn{c}, the quadratic coefficient of the mean function of the SIBLLE.
#' When 'type' = "regression", 'test' = "MLE" is assumed by default.
#'
#' When 'type' = "LAN", 'test' can be "moments", "MLE", "information", or "parameter".
#' If 'test' is "moments", "MLE", or "information", the output is the same as in the case where 'type' = "regression".
#' If 'test' is "parameter", a test about the value of the model parameter is conducted under the local asymptotic normality assumption.
#' When 'type' = "LAN", 'test' = "parameter" is assumed by default.
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified.  The weights can either be given as an attribute 'weights' of the 'siblle' object, or as a function argument 'weights', with the latter being used when both are supplied. In either case, 'weights' should be a numeric vector of length equal to that of 'siblle'. If 'weights' is given as an argument to the "ht" function, it can be specified alternatively as a character string "tricube". In this case, the tricube weight (see Cleveland, 1979) is used, and the specified 'fraction' of the points will have nonzero weights. The 'center' argument determines at which parameter value the tricube weight takes the maximum. If weights are not supplied in either location, all weights are taken to be equal to 1.
#' It is important to note that the weights are not supposed to be normalized (i.e., sum to one). Multiplying all weights by the same constant changes the local regression results. Roughly speaking, the variance of the simulation based log likelihood estimator is assumed to be sigma^2/(the weight for the point). See Park (2023) for more information.
#'
#' @return A list consisting of the following components are returned.
#' \itemize{
#' \item{meta_model_MLE_for_*: simulation based maximum likelihood estimate under an appropriate meta model,}
#' \item{Hypothesis_Tests: a data frame of the null values and the corresponding (approximate) p-values,}
#' \item{pvalue_numerical_error_size: When 'test'="moments", approximate size of error in numerical evaluation of p-values (automatically set to approximately 0.01 or 0.001). For these case, p-values are found using the SCL distributions, whose cumulative distribution functions are numerically evaluated using random number generations. Thus p-values have some stochastic error. The size of the numerical error is automatically set to approximately 0.01, but if p-value found is less than 0.01 for any of the provided null values, more computations are carried out to reduce the numerical error size to approximately 0.001. Note that when 'test'="MLE", "information", or "parameter", the (standard) F distribution is used, so this list component is omitted.}
#' }
#' When 'test' = "moments", exact p-values are shown (here "exact" means that the formula for the p-value is not based on approximation; this does not mean that size of the numerical evaluation is equal to zero.)
#' In other cases, approximate p-values are shown. See Park, J. (2023) for how approximations are made.
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @references Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American statistical association, 74(368), 829-836.
#' @references Le Cam, L. and Yang, G. L. (2000). Asymptotics in statistics: some basic concepts. Springer-Verlag, New York.
#' @export
ht.siblle <- function(siblle, null.value, type=NULL, test=NULL, weights=NULL, fraction=NULL, center=NULL) {
    # TODO: combine ht and ci
    # TODO: remove a section on observed Fisher information and do the same for ci
    # TODO: remove the tricube weight and the same for ci
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

    if (!is.null(test)) {
        match.arg(test, c("moments", "MLE", "information", "parameter"))
    }
    if (is.null(test)) {
        if (type=="point") {
            test <- "moments"
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

    if (is.numeric(null.value)) {
        if (test=="moments") {
            if (type=="point" && length(null.value)!=2) {
                stop(
                    "If 'test' is 'moments', 'type' is 'point', and 'null.value' is not a list, the length of 'null.value' should be 2 (mean and variance of 'siblle').",
                    call. = FALSE
                )
            }
            if (type %in% c("regression", "LAN") && length(null.value)!=4) {
                stop(
                    "If 'test' is 'moments', and 'type' is 'regression' or 'LAN', and 'null.value' is not a list, the length of 'null.value' should be 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the SIBLLE).",
                    call. = FALSE
                )
            }
        }
        if (test=="MLE" && length(null.value)!=1) {
         stop(
                "If 'test' is 'MLE' and 'null.value' is not a list, the length of 'null.value' should be 1.",
                call. = FALSE
            )
        }
        if (test=="information" && length(null.value)!=1) {
         stop(
                "If 'test' is 'information' and 'null.value' is not a list, the length of 'null.value' should be 1.",
                call. = FALSE
            )
        }
        if (test=="parameter" && length(null.value)!=1) {
            stop(
                "If 'test' is 'parameter' and 'null.value' is not a list, the length of 'null.value' should be 1.",
                call. = FALSE
            )
        }
    }
    if (is.list(null.value)) {
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
        if (test=="parameter" && !all(sapply(null.value, length)==1)) {
            stop(
                "If 'null.value' is a list and 'test' is 'parameter', all components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
    }
    if (type=="point" && !test %in% c("moments")) {
        stop(
            "When 'type' = 'point', 'test' should be 'moments'.",
            call. = FALSE
        )
    }
    if (type=="regression" && !test %in% c("moments", "MLE", "information")) {
        stop(
            "When 'type' = 'regression', 'test' should be one of 'moments', 'MLE' (default), or 'information'.",
            call. = FALSE
        )
    }
    if (type=="LAN" && !test %in% c("moments", "MLE", "information", "parameter")) {
        stop(
            "When 'type' = 'LAN', 'test' should be one of 'moments', 'MLE', 'information', or 'parameter' (default).",
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
            pvalout <- pscl(teststats, M, 1, num_error_size=num.error.size)
            if (length(pvalout)==0) { # execution of pscl stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$numerical_error_size
            if (any(pval < .01)) {
                num.error.size <- 0.001
                pvalout <- pscl(teststats, M, 1, num_error_size=num.error.size)
                if (length(pvalout)==0) { # execution of pscl stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                prec <- pvalout$numerical_error_size
            }
            precdigits <- max(-floor(log10(num.error.size)), 1) + 1
            dfout <- data.frame(
                mu_null=sapply(null.value, function(x) x[1]),
                sigma_sq_null=sapply(null.value, function(x) x[2]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE_for_moments=c(mu=muhat, sigma_sq=(M-1)/M*Ssq),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
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
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%llest))
        resids <- llest - c(theta012%*%Ahat)
        sigsqhat <- c(resids%*%W%*%resids) / M
        theta0123 <- cbind(1, theta, theta^2, theta^3) # cubic regression to test whether the cubic coefficient = 0
        Ahat_cubic <- c(solve(t(theta0123)%*%W%*%theta0123, t(theta0123)%*%W%*%llest))
        resids_cubic <- llest - c(theta0123%*%Ahat_cubic)
        sigsqhat_cubic <- c(resids_cubic%*%W%*%resids_cubic) / M
        pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-4), 1, sum(w>0)-4, lower.tail=FALSE)
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
                    .5*M*log(sigsqhat/x[4]) - .5*c(err%*%W%*%err)/x[4] + M/2
                })
            num.error.size <- 0.01
            pvalout <- pscl(teststats, M, 3, num_error_size=num.error.size)
            if (length(pvalout)==0) { # execution of pscl stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$numerical_error_size
            if (any(pval < .01)) {
                prec <- 0.001
                pvalout <- pscl(teststats, M, 3, num_error_size=num.error.size)
                if (length(pvalout)==0) { # execution of pscl stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$numerical_error_size
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                a_null=sapply(null.value, function(x) x[1]),
                b_null=sapply(null.value, function(x) x[2]),
                c_null=sapply(null.value, function(x) x[3]),
                sigma_sq_null=sapply(null.value, function(x) x[4]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE_for_moments=c(a=Ahat[1], b=Ahat[2], c=Ahat[3], sigma_sq=sigsqhat),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size,
                pval_cubic=pval_cubic
            )
            return(out)
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
                    (M-3)*(Ahat[2]+2*x*Ahat[3])^2*(v11*v22-v12^2)/(v22-4*v12*x+4*v11*x*x)/(M*sigsqhat)
                })
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                MLE_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_MLE=c(MLE=unname(-Ahat[2]/(2*Ahat[3]))),
                Hypothesis_Tests=dfout,
                pval_cubic=pval_cubic
            )
            return(out)
        }
        ## test about Fisher information
        if (test=="information") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            U <- t(theta012)%*%W%*%theta012
            u3gv12 <- U[3,3] - c(U[3,1:2]%*%solve(U[1:2,1:2], U[1:2,3])) # u_{3|12}
            teststats <- sapply(null.value, function(x) (M-3)*u3gv12*(Ahat[3]+x/2)^2/(M*sigsqhat))
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                information_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_observed_Fisher_information=c(observed_Fisher_information=unname(-2*Ahat[3])),
                Hypothesis_Tests=dfout,
                pval_cubic=pval_cubic
            )
            return(out)
        }
        ## test about the model parameter under LAN
        if (type=="LAN" && test=="parameter") {
            Winv <- diag(1/w)
            K_1 <- -2*Ahat[3] # first stage estimate of K
            resids_1 <- llest - c(theta012%*%Ahat) # first stage estimates for residuals
            sigsq_1 <- c(resids_1%*%W%*%resids_1) / M # the first stage estimate of sigma^2
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
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
            teststats <- sapply(null.value,
                function(x) {
                    v <- c(R_1%*%c(x, -1/2))
                    (M-3)*(sum(((diag(rep(1,M-1))-outer(v,v)/sum(v*v))%*%(sqrtQ_1%*%C%*%llest+sigsq_1/2*u))^2)/(M-1)/sigsqhat - 1)
                })
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                parameter_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_parameter=c(parameter=thetahat, information=Khat, error_variance=sigsqhat),
                Hypothesis_Tests=dfout,
                pval_cubic=pval_cubic
            )
            return(out)
        }
    }
}

