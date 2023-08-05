#' @export
ht <- function(x, ...) {
    UseMethod("ht")
}

#' Hypothesis tests based on simulation log likelihoodsbased log likelihood estimates
#'
#' `ht` outputs results of hypothesis tests carried out using simulation log likelihoods. See Park (2023) for more information.
#'
#' @name ht
#' @param simll A class 'simll' object, containing simulation log likelihoods, the parameter values at which simulations are made, and the number of simulations made at each parameter value. See help(simll).
#' @param null.value The null value(s) for the hypothesis test. Either a numeric vector (for running a test for a single null value, which can have one or more components) or a list of numeric vectors (for running tests for multiple null values).
#' @param test A character string indicating the quantity to be tested about. One of "moments", "MESLE", or "parameter". See Details.
#' @param type When 'test' is "moments", the 'type' argument needs to be specified. 'type' = "point" means that the test about the mean and the variance of simulation log likelihoods at a given parameter point is considered. 'type' = "regression" means that the test about the mean function and the variance of simulation log likelihoods at various parameter values is considered. See Details.
#' @param weights An optional argument. The un-normalized weights of the log likelihood estimates for regression. Either a numeric vector of length equal to the 'params' attribute of the 'simll' object. See Details below.
#'
#' @details
#' This is a generic function, taking a class 'simll' object as the first argument.
#' Hypothesis tests are carried out under a normal meta model--that is, the simulation log likelihoods (whose values are given in the 'simll' object) are normally distributed.
#'
#' When 'null.value' is a list, a hypothesis test is carried out for each null value specified in the list. For example, in order to run tests for more than one null values for a single-component parameter, you can let 'null.value=as.list(vector_of_null_values)'.
#'
#' If 'test' = "moments", the 'type' argument needs to be either "point" or "regression".
#' If 'type' = "point", a test about the mean and the variance of the simulation log likelihood at a single parameter value is conducted.
#' The 'null.value' should be a numeric vector of length two (the first component being the mean and the second being the variance), or a list of numeric vectors of length two.
#' If 'type' = "regression", the 'simll' object should contain simulation log likelihoods obtained at more than one parameter values, specified by the 'params' attribute of the 'simll' object. A (weighted) quadratic regression for the simulation log likelihoods will be used for hypothesis tests, where the x-axis values are given by the 'params' values of the 'simll' object and the y-axis values are the corresponding simulation log likelihoods.
#' The test is about the quadruple \eqn{a, b, c, sigma^2} where \eqn{a, b, c} are coefficients of the polynomial describing the mean of the simulation log likelihood (i.e., \eqn{l(\theta) = a + b \theta + c \theta^2}) and \eqn{\sigma^2} is the variance of the simulation log likelihood.
#' If 'test' = "moments" and 'type' is not specified, 'type' defaults to "point" if the 'params' attribute of the 'simll' object is not supplied or has length one, and defaults to "regression" otherwise.
#'
#' When 'test' = "MESLE" or "parameter", the 'simll' object should have the 'params' attribute.
#'
#' If 'test' = "MESLE", the test is about the location of the maximum expected simulation log likelihood estimate.
#'
#' If 'test' = "parameter", inference on the model parameter will be carried out under the local asymptotic normality for simulation log likelihood (see Park (2023) for more information.)
#'
#' The default value for 'test' is "parameter".
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified. The length of 'weights' should be equal to that of the 'params' attribute of the 'simll' and the number of rows in the simulation log likelihood matrix in the 'simll' object. It is important to note that the weights are not supposed to be normalized (i.e., sum to one). Multiplying all weights by the same constant changes the estimation outputs.
#'
#' @return A list consisting of the following components are returned.
#' \itemize{
#' \item{meta_model_MLE_for_*: point estimate for the tested quantity under a normal meta model,}
#' \item{Hypothesis_Tests: a data frame of the null values and the corresponding p-values,}
#' \item{pvalue_numerical_error_size: When 'test'="moments", approximate size of error in numerical evaluation of p-values (automatically set to approximately 0.01 or 0.001). For these case, p-values are found using the SCL distributions, whose cumulative distribution functions are numerically evaluated using random number generations. Thus p-values have some stochastic error. The size of the numerical error is automatically set to approximately 0.01, but if p-value found is less than 0.01 for any of the provided null values, more computations are carried out to reduce the numerical error size to approximately 0.001. Note that when 'test'="MESLE", "information", or "parameter", the (standard) F distribution is used, so this list component is omitted.}
#' }
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @export
ht.simll <- function(simll, null.value, test=NULL, type=NULL, weights=NULL) {
    validate_simll(simll)
    if (is.null(test)) {
        test <- "parameter"
        message("The 'test' argument is not supplied. Defaults to 'parameter'.")
    }
    if (!is.null(test)) {
        match.arg(test, c("moments", "MESLE", "parameter"))
    }
    if (test=="moments") {
        if (is.null(type)) {
            if (is.null(attr(simll, "params")) || length(attr(simll, "params"))==1) {
                type="point"
            } else {
                type="regression"
            }
        }
        match.arg(type, c("point", "regression"))
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
            if (type=="regression" && length(null.value)!=4) {
                stop(
                    "If 'test' is 'moments', 'type' is 'regression', and 'null.value' is not a list, the length of 'null.value' should be 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the SIBLLE).",
                    call. = FALSE
                )
            }
        }
        if (test=="MESLE" && length(null.value)!=1) {
         stop(
                "If 'test' is 'MESLE' and 'null.value' is not a list, the length of 'null.value' should be 1.",
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
            if (type=="regression" && !all(sapply(null.value, length)==4)) {
                stop(
                    "If 'null.value' is a list, 'test' is 'moments', and 'type' is 'regression' or 'LAN', all components of 'null.value' should be a numeric vector of length 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the SIBLLE).",
                    call. = FALSE
                )
            }
        }
        if (test=="MESLE" && !all(sapply(null.value, length)==1)) {
            stop(
                "If 'null.value' is a list and 'test' is 'MESLE', all list components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
        if (test=="parameter" && !all(sapply(null.value, length)==1)) {
            stop(
                "If 'null.value' is a list and 'test' is 'parameter', all list components of 'null.value' should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
    }

    if (test=="moments" && type=="point") {
        llmat <- unclass(simll)
        ll <- apply(llmat, 1, sum) # simulation log likelihood for y_{1:n}
        muhat <- mean(ll)
        Ssq <- var(ll)
        M <- dim(ll)[1]
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
    if ((test=="moments" && type=="regression") || test=="MESLE" || test=="parameter") {
        ## set weights (vector w)
        if (!is.null(weights)) {
            if (!is.numeric(weights)) {
                stop(
                    "When 'type' = 'regression' and the 'weights' argument is given, 'weights' have to be a numeric vector.",
                    call. = FALSE
                )
            }
            if (length(weights) != dim(simll)[1]) {
                stop(
                    "When 'type' = 'regression' and the 'weights' argument is given, the length of 'weights' should be equal to the number of rows in the simulation log likelihood matrix in 'simll'.",
                    call. = FALSE
                )
            }
            w <- weights
        }
        if (is.null(weights)) {
            if (!is.null(attr(simll, "weights"))) {
                if (!is.numeric(attr(simll, "weights"))) {
                    stop(
                        "When 'type' = 'regression' and the 'simll' object has 'weights' attribute, it has to be a numeric vector.",
                        call. = FALSE
                    )
                }
                if (dim(simll)[1] != length(attr(simll, "weights"))) {
                    stop(
                        "When 'type' = 'regression' and the 'simll' object has 'weights' attribute, the length of 'weights' should be the same as the number of rows in the simulation log likelihood matrix in 'simll'.",
                        call. = FALSE
                    )
                }
                w <- attr(simll, "weights")
            } else {
                w <- rep(1, dim(simll)[1])
            }
        }
        ## weighted quadratic regression
        W  <- diag(w)
        theta <- attr(simll, "params")
        llmat <- unclass(simll)
        ll <- apply(llmat, 1, sum)
        M <- length(ll)
        theta012 <- cbind(1, theta, theta^2)
        Ahat <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%ll))
        resids <- ll - c(theta012%*%Ahat)
        sigsqhat <- c(resids%*%W%*%resids) / M
        theta0123 <- cbind(1, theta, theta^2, theta^3) # cubic regression to test whether the cubic coefficient = 0
        Ahat_cubic <- c(solve(t(theta0123)%*%W%*%theta0123, t(theta0123)%*%W%*%ll))
        resids_cubic <- ll - c(theta0123%*%Ahat_cubic)
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
                    err <- ll - c(theta012%*%x[1:3])
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
        ## test about MESLE
        if (test=="MESLE") {
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
                MESLE_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_MESLE=c(MESLE=unname(-Ahat[2]/(2*Ahat[3]))),
                Hypothesis_Tests=dfout,
                pval_cubic=pval_cubic
            )
            return(out)
        }
        ## test about the model parameter under LAN
        if (test=="parameter") {
            Winv <- diag(1/w)
            uniquetheta <- unique(theta)
            uniquetheta012 <- cbind(1, uniquetheta, uniquetheta^2)
            nobs <- dim(llmat)[2] # number of observations
            slopes <- apply(llmat, 2, function(ll_i) {
                # quadratic regression for each observation piece (each column of simll)
                Ahat_i <- c(solve(t(theta012)%*%W%*%theta012, t(theta012)%*%W%*%ll_i))
                return(c(uniquetheta012%*%Ahat_i))
            }) # length(uniquetheta) times nobs
            K1_vec <- apply(slopes, 1, var)
            K1hat <- median(K1_vec)
            resids_1 <- ll - c(theta012%*%Ahat) # first stage estimates for residuals
            sigsq_1 <- c(resids_1%*%W%*%resids_1) / M # the first stage estimate of sigma^2
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
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
            teststats <- sapply(null.value,
                function(x) {
                    v <- c(R_1%*%c(x, -1/2))
                    (M-3)*(sum(((diag(rep(1,M-1))-outer(v,v)/sum(v*v))%*%sqrtQ_1%*%C%*%ll)^2)/(M-1)/sigsqhat - 1)
                })
            pval <- pf(teststats, 1, M-3, lower.tail=FALSE)
            dfout <- data.frame(
                parameter_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_parameter=c(parameter=thetastarhat, error_variance=sigsqhat),
                Hypothesis_Tests=dfout,
                pval_cubic=pval_cubic
            )
            return(out)
        }
    }
}
# TODO: debugging required.
