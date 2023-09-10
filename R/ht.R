#' @export
ht <- function(x, ...) {
    UseMethod("ht")
}

#' Hypothesis tests based on simulation log likelihoodsbased log likelihood estimates
#'
#' `ht` outputs results of hypothesis tests carried out using simulation log likelihoods. See Park (2023) for more information.
#'
#' @name ht
#' @param simll A class `simll` object, containing simulation log likelihoods, the parameter values at which simulations are made (optional), and the weights for those simulations for regression (optional). See help(simll).
#' @param null.value The null value(s) for the hypothesis test. Either a numeric vector (for running a test for a single null value, which can have one or more components) or a list of numeric vectors (for running tests for multiple null values).
#' @param test A character string indicating which is to be tested about. One of "moments", "MESLE", or "parameter". See Details.
#' @param case When `test` is "parameter", `case` needs to be either "iid" or "stationary". `case` = "iid" means that the observations are iid, and `case` = "stationary" means that the observations form a stationary sequence. The `case` argument affects how the variance of the slope of the mean function (=K_1 in Park (2023)) is estimated. The default value is "stationary".
#' @param type When `test` is "moments", the `type` argument needs to be specified. `type` = "point" means that the test about the mean and the variance of simulation log likelihoods at a given parameter point is considered. `type` = "regression" means that the test about the mean function and the variance of simulation log likelihoods at various parameter values is considered. See Details.
#' @param weights An optional argument. The un-normalized weights of the simulation log likelihoods for regression. A numeric vector of length equal to the `params` attribute of the `simll` object. See Details below.
#' @param ncores An optional argument indicating the number of CPU cores to use for computation. Used only when `test`="parameter". The `mclapply` function in the `parallel` package is used. The `parallel` package needs to be installed unless `ncores`=1. In Windows, `ncores` greater than 1 is not supported (see ?mclapply for more information.)
#'
#' @details
#' This is a generic function, taking a class `simll` object as the first argument.
#' Hypothesis tests are carried out under a normal meta model--that is, the simulation log likelihoods (whose values are given in the `simll` object) are normally distributed.
#'
#' When `null.value` is a list, a hypothesis test is carried out for each null value specified in the list. For example, in order to run tests for more than one null values for a single-component parameter, you can let `null.value=as.list(vector_of_null_values)`.
#'
#' If `test` = "moments", the `type` argument needs to be either "point" or "regression".
#' If `type` = "point", a test about the mean and the variance of the simulation log likelihood at a single parameter value is conducted.
#' The `null.value` should be a numeric vector of length two (the first component being the mean and the second being the variance), or a list of numeric vectors of length two.
#' If `type` = "regression", the `simll` object should contain simulation log likelihoods obtained at more than one parameter values, specified by the `params` attribute of the `simll` object. A (weighted) quadratic regression for the simulation log likelihoods will be used for hypothesis tests, where the x-axis values are given by the `params` values of the `simll` object and the y-axis values are the corresponding simulation log likelihoods.
#' The test is about the quadruple \eqn{a, b, c, sigma^2} where \eqn{a, b, c} are coefficients of the polynomial describing the mean of the simulation log likelihood (i.e., \eqn{l(\theta) = a + b \theta + c \theta^2}) and \eqn{\sigma^2} is the variance of the simulation log likelihood.
#' If `test` = "moments" and `type` is not specified, `type` defaults to "point" if the `params` attribute of the `simll` object is not supplied or has length one, and defaults to "regression" otherwise.
#'
#' When `test` = "MESLE" or "parameter", the `simll` object should have the `params` attribute.
#'
#' If `test` = "MESLE", the test is about the location of the maximum expected simulation log likelihood estimate.
#'
#' If `test` = "parameter", inference on the simulation based surrogate will be carried out under the local asymptotic normality for simulation log likelihood (see Park (2023) for more information.)
#'
#' The default value for `test` is "parameter".
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified. The length of `weights` should be equal to that of the `params` attribute of the `simll`, which is equal to the number of rows in the simulation log likelihood matrix in the `simll` object. It is important to note that the weights are not supposed to be normalized (i.e., sum to one). Multiplying all weights by the same constant changes the estimation outputs. If not supplied, the `weights` attribute of the `simll` object is used. If neither is supplied, `weights` defaults to the vector of all ones.
#'
#' @return A list consisting of the following components are returned.
#' \itemize{
#' \item{meta_model_MLE_for_*: point estimate for the tested quantity under a normal meta model,}
#' \item{Hypothesis_Tests: a data frame of the null values and the corresponding p-values,}
#' \item{pvalue_numerical_error_size: When `test`="moments", approximate size of error in numerical evaluation of p-values (automatically set to approximately 0.01 or 0.001). For these case, p-values are found using the SCL distributions, whose cumulative distribution functions are numerically evaluated using random number generations. Thus p-values have some stochastic error. The size of the numerical error is automatically set to approximately 0.01, but if p-value found is less than 0.01 for any of the provided null values, more computations are carried out to reduce the numerical error size to approximately 0.001. Note that when `test`="MESLE", "information", or "parameter", the (standard) F distribution is used, so this list component is omitted.}
#' \item{pval_cubic: The p-value of the test about whether the cubic term in the cubic polynomial regression is significant. If so, the result of the ht function may be biased.}
#' }
#'
#' @references Park, J. (2023). On simulation based inference for implicitly defined models
#' @export
ht.simll <- function(simll, null.value, test=NULL, case=NULL, type=NULL, weights=NULL, ncores=1) {
    validate_simll(simll)
    if (is.null(test)) {
        test <- "parameter"
        message("The `test` argument is not supplied. Defaults to `parameter`.")
    }
    if (!is.null(test)) {
        match.arg(test, c("moments", "MESLE", "parameter"))
    }
    if (test=="parameter") {
        if (is.null(case)) {
            case <- "stationary"
        match.arg(case, c("iid", "stationary"))
        }
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
                "`null.value` should be numeric or a list of numeric values",
                call. = FALSE
            )
        }
        if (is.list(null.value) && !all(sapply(null.value, is.numeric))) {
            stop(
                "If `null.value` is a list, all of its components should be numeric.",
                call. = FALSE
            )
        }
    }

    if (is.numeric(null.value)) {
        if (test=="moments") {
            if (type=="point" && length(null.value)!=2) {
                stop(
                    "If `test` is 'moments', `type` is 'point', and `null.value` is not a list, the length of `null.value` should be 2 (mean and variance of the simulation log likelihoods).",
                    call. = FALSE
                )
            }
            if (type=="regression" && length(null.value)!=4) {
                ## TODO: change to multiparameter case
                stop(
                    "If `test` is 'moments', `type` is 'regression', and `null.value` is not a list, the length of `null.value` should be 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the simulation log likelihood).",
                    call. = FALSE
                )
            }
        }
        if (test=="MESLE" && length(null.value)!=1) {
         stop(
                "If `test` is 'MESLE' and `null.value` is not a list, the length of `null.value` should be 1.",
                call. = FALSE
            )
        }
        if (test=="parameter" && length(null.value)!=1) {
            stop(
                "If `test` is 'parameter' and `null.value` is not a list, the length of `null.value` should be 1.",
                call. = FALSE
            )
        }
    }
    if (is.list(null.value)) {
        if (test=="moments") {
            if (type=="point" && !all(sapply(null.value, length)==2)) {
                stop(
                    "If `null.value` is a list, `test` is 'moments', and `type` is 'point', all components of `null.value` should be a numeric vector of length 2 (mean and variance of simulation log likelihoods).",
                    call. = FALSE
                )
            }
            if (type=="regression" && !all(sapply(null.value, length)==4)) {
                stop(
                    "If `null.value` is a list, `test` is 'moments', and `type` is 'regression', all components of `null.value` should be a numeric vector of length 4 (the three coefficients of a quadratic polynomial for the mean function, and the variance of the simulation log likelihoods).",
                    call. = FALSE
                )
            }
        }
        if (test=="MESLE" && !all(sapply(null.value, length)==1)) {
            stop(
                "If `null.value` is a list and `test` is 'MESLE', all list components of `null.value` should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
        if (test=="parameter" && !all(sapply(null.value, length)==1)) {
            stop(
                "If `null.value` is a list and `test` is 'parameter', all list components of `null.value` should be a numeric vector of length 1.",
                call. = FALSE
            )
        }
    }

    if (test=="moments" && type=="point") {
        llmat <- unclass(simll)
        ll <- apply(llmat, 2, sum) # simulation log likelihood for y_{1:n}
        muhat <- mean(ll)
        Ssq <- var(ll)
        M <- length(ll)
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
                    "When `type` = 'regression' and the `weights` argument is given, `weights` have to be a numeric vector.",
                    call. = FALSE
                )
            }
            if (length(weights) != dim(simll)[2]) {
                stop(
                    "When `type` = 'regression' and the `weights` argument is given, the length of `weights` should be equal to the number of rows in the simulation log likelihood matrix in `simll`.",
                    call. = FALSE
                )
            }
            w <- weights
        }
        if (is.null(weights)) {
            if (!is.null(attr(simll, "weights"))) {
                if (!is.numeric(attr(simll, "weights"))) {
                    stop(
                        "When the `simll` object has `weights` attribute, it has to be a numeric vector.",
                        call. = FALSE
                    )
                }
                if (dim(simll)[2] != length(attr(simll, "weights"))) {
                    stop(
                        "When the `simll` object has `weights` attribute, the length of `weights` should be the same as the number of rows in the simulation log likelihood matrix in `simll`.",
                        call. = FALSE
                    )
                }
                w <- attr(simll, "weights")
            } else {
                w <- rep(1, dim(simll)[2])
            }
        }
        ## weighted quadratic regression
        vech <- function(mat) { # half-vectorization
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
        W  <- diag(w)
        theta <- cbind(attr(simll, "params")) # coerce into a matrix
        llmat <- unclass(simll)
        ll <- apply(llmat, 2, sum)
        M <- length(ll)
        Theta012 <- t(apply(theta, 1, vec012))
        Ahat <- c(solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%ll))
        ahat <- Ahat[1]
        d <- dim(theta)[2]
        bindex <- 2:(d+1) # the positions in A that correspond to b
        bhat <- Ahat[bindex]
        cindex <- (d+2):((d^2+3*d+2)/2) # the positions in A that correspond to vech(c)
        vech_chat <- Ahat[cindex]
        chat <- unvech(vech_chat)
        resids <- ll - c(Theta012%*%Ahat)
        sigsqhat <- c(resids%*%W%*%resids) / M
        MESLEhat <- unname(-solve(chat,bhat)/2)
        ## cubic test: not included because it can be cumbersome for d>1.
        if (d==1) {
            Theta0123 <- cbind(1, theta, theta^2, theta^3) # cubic regression to test whether the cubic coefficient = 0
            Ahat_cubic <- c(solve(t(Theta0123)%*%W%*%Theta0123, t(Theta0123)%*%W%*%ll))
            resids_cubic <- ll - c(Theta0123%*%Ahat_cubic)
            sigsqhat_cubic <- c(resids_cubic%*%W%*%resids_cubic) / M
            pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-4), 1, sum(w>0)-4, lower.tail=FALSE)
        }
        ## test about moments
        if (test=="moments") {
            ## TODO: check if the moments test are correct
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            teststats <- sapply(null.value,
                function(x) {
                    err <- ll - c(Theta012%*%x[1:((d^2+3*d+2)/2)])
                    .5*M*log(sigsqhat/x[(d^2+3*d+2)/2+1]) - .5*c(err%*%W%*%err)/x[(d^2+3*d+2)/2+1] + M/2
                })
            num.error.size <- 0.01
            pvalout <- pscl(teststats, M, (d^2+3*d+2)/2, num_error_size=num.error.size)
            if (length(pvalout)==0) { # execution of pscl stopped by user input
                stop("Hypothesis tests stopped by user input", call. = FALSE)
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$numerical_error_size
            if (any(pval < .01)) {
                prec <- 0.001
                pvalout <- pscl(teststats, M, (d^2+3*d+2)/2, num_error_size=num.error.size)
                if (length(pvalout)==0) { # execution of pscl stopped by user input
                    stop("Hypothesis tests stopped by user input", call. = FALSE)
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$numerical_error_size
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- data.frame(
                a_null=sapply(null.value, function(x) x[1]),
                b_null=sapply(null.value, function(x) x[2:(d+1)]),
                c_null=sapply(null.value, function(x) x[(d+2):((d^2+3*d+2)/2)]),
                sigma_sq_null=sapply(null.value, function(x) x[(d^2+3*d+2)/2+1]),
                pvalue=round(pval, digits=precdigits)
            )
            out <- list(meta_model_MLE_for_moments=c(a=Ahat[1], b=Ahat[2:(d+1)], c=Ahat[(d+2):((d^2+3*d+2)/2)], sigma_sq=sigsqhat),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
            )
            return(out)
        }
        ## test about MESLE
        if (test=="MESLE") {
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            U <- t(Theta012)%*%W%*%Theta012
            V <- U[-1,-1] - U[-1,1,drop=FALSE]%*%U[1,-1,drop=FALSE]/U[1,1]
            #mtheta1 <- sum(w*theta)/sum(w)
            #mtheta2 <- sum(w*theta*theta)/sum(w)
            #mtheta3 <- sum(w*theta*theta*theta)/sum(w)
            #mtheta4 <- sum(w*theta*theta*theta*theta)/sum(w)
            #v11 <- sum(w)*(mtheta2 - mtheta1^2)
            #v12 <- sum(w)*(mtheta3 - mtheta1*mtheta2)
            #v22 <- sum(w)*(mtheta4 - mtheta2^2)
            teststats <- sapply(null.value,
                function(x) {
                    theta0mat <- matricize(null.value)
                    Vq <- solve(cbind(diag(d), 2*theta0mat) %*% solve(V) %*% rbind(diag(d), 2*t(theta0mat)))
                    xi <- t(bhat + 2*theta0mat %*% vech_chat) %*% Vq %*% (bhat + 2*theta0mat %*% vech_chat)
                    (M-(d^2+3*d+2)/2)*xi/(M*d*sigsqhat)
                })
            pval <- pf(teststats, d, M-(d^2+3*d+2)/2, lower.tail=FALSE)
            dfout <- data.frame(
                MESLE_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(meta_model_MLE_for_MESLE=c(MESLE=MESLEhat),
                Hypothesis_Tests=dfout
            )
            return(out)
        }
        ## test on the simulation based surrogate under LAN
        if (test=="parameter") {
            Winv <- diag(1/w)
            nobs <- dim(llmat)[1] # number of observations
            if (all(sapply(1:d, function(k) { min(theta[,k]) <= MESLEhat[k] && MESLEhat[k] <= max(theta[,k]) }))) {
                ## if MESLEhat is within the range of the given theta matrix componentwise, find the slope at MESLEhat
                slope_at <- MESLEhat
            } else {
                ## otherwise, find the slope at the componentwise mean of the theta matrix
                slope_at <- apply(MESLEhat, 2, mean)
            }
            if (ncores>1) {
                require(parallel)
                slope <- simplify2array(mclapply(1:nobs, function(i) {
                    ## quadratic regression for each observation piece (each row of simll)
                    ll_i <- llmat[i,]
                    Ahat_i <- c(solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%ll_i))
                    return(Ahat_i[bindex]+2*unvech(Ahat_i[cindex])%*%slope_at) # estimated slope
                }, mc.cores=ncores)
                ) # matrix of dimension d times nobs
            } else {
                slope <- sapply(1:nobs, function(i) {
                    ## quadratic regression for each observation piece (each row of simll)
                    ll_i <- llmat[i,]
                    Ahat_i <- c(solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%ll_i))
                    return(Ahat_i[bindex]+2*unvech(Ahat_i[cindex])%*%slope_at) # estimated slope
                } ) # matrix of dimension d times nobs
            }
            slope <- rbind(slope) # if d = 1, make `slope` into a 1 X nobs matrix
            var_slope_vec <- var(t(slope)) # estimate of the variance of the slope of the fitted quadratic
            E_condVar_slope <- c(cbind(0,diag(d),2*matricize(slope_at))%*%solve(t(Theta012)%*%W%*%Theta012, t(cbind(0,diag(d),2*matricize(slope_at))))) * sigsqhat / nobs # an estimate of the expected value of the conditional variance of the estimated slope given Y (expectation taken with respect to Y)
            K1hat <- var_slope_vec - E_condVar_slope # estimate of K1
            if (any(svd(K1hat)$d <= 0)) {
                stop("The estimate of K1 is not positive definite. Hypothesis test stopped.")
            }
            resids_1 <- ll - c(Theta012%*%Ahat) # first stage estimates for residuals
            sigsq_1 <- c(resids_1%*%W%*%resids_1) / M # the first stage estimate of sigma^2
            if (!is.list(null.value)) {
                null.value <- list(null.value)
            }
            C <- cbind(-1, diag(rep(1,M-1)))
            Theta12 <- Theta012[,-1]
            Ctheta <- C%*%theta
            Q_1 <- solve(C%*%Winv%*%t(C) + nobs/sigsq_1*Ctheta%*%K1hat%*%t(Ctheta)) 
            svdQ_1 <- svd(Q_1)
            sqrtQ_1 <- svdQ_1$u %*% diag(sqrt(svdQ_1$d)) %*% t(svdQ_1$v)
            R_1 <- sqrtQ_1%*%C%*%Theta12
            estEq_2 <- c(solve(t(R_1)%*%R_1, t(R_1)%*%(sqrtQ_1%*%C%*%ll))) # estimating equation for K2 and theta_star. (thetastarhat // -I/2) * n * vech(K2hat) = (R_1^T R_1)^{-1} R_1^T (Q_1^{1/2} C lS)
            vech_K2hat <- -2*estEq_2[(d+1):((d^2+3*d)/2)]/nobs # second stage estimate of vech(K2)
            K2hat <- unvech(vech_K2hat)
            thetastarhat <- solve(K2hat,estEq_2[1:d])/nobs # maximum meta model likelihood estimate for theta_star (simulation based surrogate)
            sigsqhat_lan <- 1/(M-1)*sum((sqrtQ_1%*%C%*%ll - R_1%*%estEq_2)^2)
            ##varGamma <- 1; varLogGamma <- 1.645; covGammaLogGamma <- 1; # variance and covariance of Gamma(1,1) and log(Gamma(1,1)) ## TODO: remove this line
            ##theta_true <- 1 ## TODO: remove this line
            ##sigsq_true <- theta_true^(-2)*varGamma + sum(y^2)*varLogGamma - 2/theta_true*sum(y)*covGammaLogGamma ## true sigma^2, for this gamma-Poisson model. TODO: remove this line
            ##K2_true <- 1 # true K2 for this model.  TODO: remove this line
            teststats <- sapply(null.value,
                function(x) {
                    desgmat <- rbind(-2*matricize(x), diag((d^2+d)/2))
                    G <- R_1%*%desgmat%*%solve(t(desgmat)%*%t(R_1)%*%R_1%*%desgmat, t(desgmat)%*%t(R_1))
                    (M-(d^2+3*d+2)/2)/d*(sum(((diag(M-1)-G)%*%sqrtQ_1%*%C%*%ll)^2)/(M-1)/sigsqhat_lan - 1)
                })
            pval <- pf(teststats, d, M-(d^2+3*d+2)/2, lower.tail=FALSE)
            dfout <- data.frame(
                parameter_null=unlist(null.value),
                pvalue=round(pval, digits=3)
            )
            out <- list(
                meta_model_MLE_for_parameter=c(parameter=thetastarhat, K1=K1hat, K2=K2hat, error_variance=sigsqhat_lan, error_variance_non_lan=sigsqhat),
                teststats=teststats,
                Hypothesis_Tests=dfout
            )
            return(out)
        }
    }
}
    
### TODO: for test="parameter", `type` should be `iid` or `stationary`.
### The stationary case should be implemented.  A test of stationarity can be carried out, and a warning should be given when the test is positive.
### TODO: use mclapply for regression for each observation piece.
