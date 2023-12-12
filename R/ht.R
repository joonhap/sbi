#' @export
ht <- function(x, ...) {
    UseMethod("ht")
}

#' Hypothesis tests based on simulation based log likelihood estimates
#'
#' `ht` carries out hypothesis tests for models defined implicitly by a random simulator. It takes as input estimates of the log likelihood obtained via simulations of the model. Tests are carried out using a simulation meta model. See Park (2023) for more details on the method.
#'
#' @name ht
#' @param simll A class `simll` object, containing simulation log likelihoods, the parameter values at which simulations are made (may be omitted if all simulations are made at the same parameter value), and the weights for those simulations for regression (optional). See help(simll).
#' @param null.value The null value(s) for the hypothesis test. The expected format depends on which teset will be carried out. See the Details section for more information.
#' @param test A character string indicating which is to be tested about. One of "moments", "MESLE", or "parameter". See Details.
#' @param case When `test` is "parameter", `case` needs to be either "iid" or "stationary". `case` = "iid" means that the observations are iid, and `case` = "stationary" means that the observations form a stationary sequence. The `case` argument affects how the variance of the slope of the mean function (=K_1 in Park (2023)) is estimated. The default value is "stationary".
#' @param type When `test` is "moments", the `type` argument needs to be specified. `type` = "point" means that the test about the mean and the variance of simulation log likelihoods at a given parameter point is considered. `type` = "regression" means that the test about the mean function and the variance of simulation log likelihoods at various parameter values is considered. See Details.
#' @param weights An optional argument. The un-normalized weights of the simulation log likelihoods for regression. A numeric vector of length equal to the `params` attribute of the `simll` object. See Details below.
#' @param K1_est_method Either "batch" or "autocov". Used when `test` is "parameter" and `case` is "stationary". The default is "batch". See Details for more information.
#' @param batch_size Numeric. The size of the batch when `K1_est_method` is "batch". If not supplied, the default value is `round(n^0.4)` where `n` is the number of observations in the data.
#' @param max_lag When `test` is "parameter" and `case` is "stationary", the value of `max_lag` gives the truncation point for lagged autocovariance when estimating K1 as a sum of lagged autocovariances of estimates slopes. If not supplied, default is the maximum lag for which at least one of the entries of the matrix of lagged autocorrelation has absolute value greater than 4/sqrt(nobs), where the lagged autocorrelation is found up to lag `10*log10(nobs/d)`. Here `nobs` is the number of observations and `d` is the dimension of the parameter space.
#' @param plot_acf Logical. Should the autocorrelation plot be generated when estimating K1 for the case where `test` is "parameter" and `case` is "stationary"?
#' @param MCcorrection For tests on the simulation based parameter surrogate (`test`="parameter"), `MCcorrection` determines if and how the sampling distribution of the test statistic will be corrected by a Monte Carlo method to account for the variability in the estimate of K1. Possible values are "none" (default) and "Wishart". See the Details section and Park (2023) for more details.
#'
#' @details
#' This is a generic function, taking a class `simll` object as the first argument.
#' Hypothesis tests are carried out under a normal meta model--that is, the simulation log likelihoods (whose values are given in the `simll` object) are normally distributed.
#'
#' If `test` = "moments", the `type` argument needs to be either "point" or "regression".
#' If `type` = "point", a test about the mean and the variance of the simulation log likelihood at a single parameter value is conducted.
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
#' When `test` is "moments" and `type` is "point", `null.value` is either a vector of length two (one entry for the mean and the other for the variance of the simulation log likelihoods), a matrix of two columns (one for the mean and the other for the variance), or a list of vectors of length two (each entry of the list gives a null value consisting of the mean and the variance.)
#' When `test` is "moments" and `type` is "regression", `null.value` can be a list of length four, or a list of lists of length four. The first case corresponds to when a single null hypothesis is tested. The four components are a) the constant term in the quadratic mean function (scalar), b) the linear coefficient term in the mean function (vector of length \eqn{d} where \eqn{d} is the dimension of the parameter vector), c) the quadratic coefficient term in the mean function (symmetric matrix of dimension \eqn{d \times d}), and d) the variance of the simulation log likelihood (scalar). The second case is when more than one null values are tested. In this case each component of the list is a list having four entries as described for the case of a single null value.
#' When `test` is "MESLE" or "parameter", `null.value` is a vector of length \eqn{d} (a single null value), a matrix having \eqn{d} columns (each row giving a vector for a null value), or a list of vectors of length \eqn{d} (more than one null values).
#' 
#' @return A list consisting of the following components are returned.
#' \itemize{
#' \item{regression_estimates: point estimates for the meta model parameters, a, b, c, and sigma^2. Given only when test="MESLE" or "parameter".}
#' \item{meta_model_MLE_for_*: point estimate for the tested quantity under a normal meta model}
#' \item{Hypothesis_Tests: a data frame of the null values and the corresponding p-values. When `test`="moments" and `type`="regression", each null value is given in the form of c(a,b,c,sigma^2) where a, b, c, sigma^2 are first, second, third, and fourth entries of the given null value.}
#' \item{pvalue_numerical_error_size: When `test`="moments", approximate size of error in numerical evaluation of p-values (automatically set to approximately 0.01 or 0.001). For these case, p-values are found using the SCL distributions, whose cumulative distribution functions are numerically evaluated using random number generations. Thus p-values have some stochastic error. The size of the numerical error is automatically set to approximately 0.01, but if any of the p-values found is less than 0.01, more computations are carried out to reduce the numerical error size to approximately 0.001. Note that when `test`="MESLE" or "parameter", the (standard) F distribution is used, so this list component is omitted.}
#' \item{max_lag: if `test`="parameter" and `case`="stationary", the maximum lag for computing the autocovariance in estimating K1 is shown.}
#' \item{pval_cubic: The p-value of the test about whether the cubic term in the cubic polynomial regression is significant. If so, the result of the ht function may be biased. The test on the cubic term is carried out only when the number of simulation log likelihoods is greater than \eqn{(d+1)*(d+2)*(d+3)/6} where \eqn{d} is the dimension of the parameter vector.}
#' }
#'
#' @references Park, J. (2023). On simulation-based inference for implicitly defined models (https://arxiv.org/abs/2311.09446)
#' @export
ht.simll <- function(simll, null.value, test=c("parameter","MESLE","moments"), case=NULL, type=NULL, weights=NULL, K1_est_method="batch", batch_size=NULL, max_lag=NULL, plot_acf=FALSE, MCcorrection="none", ncores=1) {
    validate_simll(simll)
    match.arg(test, c("moments", "MESLE", "parameter"))
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

    if (test=="moments") {
        if (type=="point") {
            if (is.numeric(null.value)) {
                if (is.matrix(null.value)) {
                    if (ncol(null.value)!=2) {
                        stop("If `test` is 'moments', `type` is 'point', and `null.value` should be a vector of length 2, a matrix having two columns, or a list whose entries are all vector of lengths 2.")
                    }
                    null.value <- apply(null.value, 1, identity, simplify=FALSE) # make null.value into a list
                } else {
                    if (!is.null(dim(null.value))) { # if null.value is a 3 or higer dimensional array
                        stop("If `test` is 'moments', `type` is 'point', and `null.value` should be a vector of length 2, a matrix having two columns, or a list whose entries are all vector of lengths 2.")
                    }
                    if (length(null.value)!=2) {
                        stop("If `test` is 'moments', `type` is 'point', and `null.value` should be a vector of length 2, a matrix having two columns, or a list whose entries are all vector of lengths 2.")
                    }
                    null.value <- list(null.value)
                }
            } else if (is.list(null.value)) { 
                if (!all(sapply(null.value, function(n) { is.numeric(n) && length(n)==2 }))) {
                    stop("If `test` is 'moments', `type` is 'point', and `null.value` should be a vector of length 2, a matrix having two columns, or a list whose entries are all vector of lengths 2.")
                }
            } else { # if null.value is not numeric or a list
                stop("`null.value` should be a numeric vector, a matrix, or a list.")
            }
        }
        if (type=="regression") {
            if (is.null(attr(attr(simll, "params"), "dim"))) {
                d <- 1
            } else {
                d <- dim(attr(simll, "params"))[2]
            }
            if (!is.list(null.value)) {
                stop("If `test` is 'moments' and `type` is 'regression', `null.value` should be a list.")
            }
            if (length(null.value)!=4 || !all(c(is.numeric(null.value[[1]]), is.numeric(null.value[[2]]), is.numeric(null.value[[3]]), is.numeric(null.value[[4]]), length(null.value[[1]])==1, length(null.value[[2]])^2==length(null.value[[3]]), length(null.value[[4]])==1))) {
                if (!all(sapply(null.value, is.list))) {
                    stop("If `test` is 'moments' and `type` is 'regression', `null.value` should be either a list of length four or a list of lists of length four. In the first case, the first entry should be a numeric scalar, the second a numeric vector of length d where d is the dimension of the parameter space, the third a symmetric matrix of dimension d X d, and the fourth a numeric scalar. In the second case, each entry of the list should be of the same form as described for the first case.")
                } 
                if (!all(sapply(null.value, function(n) { all(c(is.numeric(n[[1]]), is.numeric(n[[2]]), is.numeric(n[[3]]), is.numeric(n[[4]]), length(n[[1]])==1, length(n[[2]])==d, length(n[[3]])==d^2, length(n[[4]])==1)) }))) {
                stop("If `test` is 'moments' and `type` is 'regression', `null.value` should be either a list of length four or a list of lists of length four. In the first case, the first entry should be a numeric scalar, the second a numeric vector of length d where d is the dimension of the parameter space, the third a symmetric matrix of dimension d X d, and the fourth a numeric scalar. In the second case, each entry of the list should be of the same form as described for the first case.")
                } 
            }
        }
    }
    if (test=="MESLE" || test=="parameter") {
        if (is.null(attr(attr(simll, "params"), "dim"))) {
            d <- 1
        } else {
            d <- dim(attr(simll, "params"))[2]
        }
        if (is.numeric(null.value)) {
            if (is.matrix(null.value)) {
                if (ncol(null.value)!=d) {
                    stop("If `test` is 'MESLE' or 'parameter', `null.value` should be a numeric vector of length d, a matrix having d columns, or a list whose entries are vector of length d, where d is the dimension of the parameter space.")
                }
                null.value <- apply(null.value, 1, identity, simplify=FALSE) # make null.value into a list
            } else {
                if (!is.null(dim(null.value))) { # if null.value is a 3 or higher dimensional array
                    stop("If `test` is 'MESLE' or 'parameter', `null.value` should be a numeric vector of length d, a matrix having d columns, or a list whose entries are vector of length d, where d is the dimension of the parameter space.")
                }
                if (length(null.value)!=d) {
                    stop("If `test` is 'MESLE' or 'parameter', `null.value` should be a numeric vector of length d, a matrix having d columns, or a list whose entries are vector of length d, where d is the dimension of the parameter space.")
                }
                null.value <- list(null.value)
            }
        } else if (is.list(null.value)) {
            if (!all(sapply(null.value, function(n) { is.numeric(n) && length(n)==d }))) {
                stop("If `test` is 'MESLE' or 'parameter', `null.value` should be a numeric vector of length d, a matrix having d columns, or a list whose entries are vector of length d, where d is the dimension of the parameter space.")
            }
        } else { # if null.value is not numeric or a list
            stop("If `test` is 'MESLE' or 'parameter', `null.value` should be a numeric vector, a matrix, or a list.")
        }
    }
    if (test=="moments" && type=="point") {
        ll <- unclass(simll)
        muhat <- mean(ll)
        Ssq <- var(ll)
        M <- length(ll)
        if (any(sapply(null.value, function(x) x[2]<=0))) {
            stop("The second component of null.value (the variance of simulation based log likelihood estimator) should be positive.")
        }
        teststats <- sapply(null.value, function(x) -.5*M*(muhat - x[1])^2/x[2] - (M-1)/2*Ssq/x[2] + M/2*log((M-1)*Ssq/(M*x[2])) + M/2)
        num.error.size <- 0.01
        pvalout <- pscl(teststats, M, 1, num_error_size=num.error.size)
        if (length(pvalout)==0) { # execution of pscl stopped by user input
            stop("Hypothesis tests stopped by user input")
        }
        pval <- pvalout$probs
        num.error.size <- pvalout$numerical_error_size
        if (any(pval < .01)) {
            num.error.size <- 0.001
            pvalout <- pscl(teststats, M, 1, num_error_size=num.error.size)
            if (length(pvalout)==0) { # execution of pscl stopped by user input
                stop("Hypothesis tests stopped by user input")
            }
            pval <- pvalout$probs
        }
        precdigits <- max(-floor(log10(num.error.size)), 1) + 1
        dfout <- data.frame(
            mu_null=sapply(null.value, function(x) x[1]),
            sigma_sq_null=sapply(null.value, function(x) x[2]),
            pvalue=pval
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
                stop("When `type` = 'regression' and the `weights` argument is given, `weights` have to be a numeric vector.")
            }
            if (length(weights) != dim(simll)[2]) {
                stop("When `type` = 'regression' and the `weights` argument is given, the length of `weights` should be equal to the number of rows in the simulation log likelihood matrix in `simll`.")
            }
            w <- weights
        }
        if (is.null(weights)) {
            if (!is.null(attr(simll, "weights"))) {
                if (!is.numeric(attr(simll, "weights"))) {
                    stop("When the `simll` object has `weights` attribute, it has to be a numeric vector.")
                }
                if (dim(simll)[2] != length(attr(simll, "weights"))) {
                    stop("When the `simll` object has `weights` attribute, the length of `weights` should be the same as the number of rows in the simulation log likelihood matrix in `simll`.")
                }
                w <- attr(simll, "weights")
            } else {
                w <- rep(1, dim(simll)[2])
            }
        }
        ## weighted quadratic regression
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
        W  <- diag(w, nrow=length(w))
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
        Ahat <- c(solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%(W%*%ll)))
        ahat <- Ahat[1]
        d <- dim(theta)[2]
        bindex <- 2:(d+1) # the positions in A that correspond to b
        bhat <- Ahat[bindex]
        cindex <- (d+2):((d^2+3*d+2)/2) # the positions in A that correspond to vech(c)
        vech_chat <- Ahat[cindex]
        chat <- unvech(vech_chat)
        ahat_b <- c(Ahat[1] - bhat%*%diag(1/theta_sd, nrow=length(theta_sd))%*%theta_mean + theta_mean%*%diag(1/theta_sd, nrow=length(theta_sd))%*%chat%*%diag(1/theta_sd, nrow=length(theta_sd))%*%theta_mean) # the constant term ahat on the original scale (transformed back)
        bhat_b <- diag(1/theta_sd, nrow=length(theta_sd))%*%bhat - 2*diag(1/theta_sd, nrow=length(theta_sd))%*%chat%*%diag(1/theta_sd, nrow=length(theta_sd))%*%theta_mean # bhat on the original scale
        chat_b <- diag(1/theta_sd, nrow=length(theta_sd))%*%chat%*%diag(1/theta_sd, nrow=length(theta_sd)) # chat on the original scale
        resids <- ll - c(Theta012%*%Ahat)
        sigsqhat <- c(resids%*%W%*%resids) / M
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
            Ahat_cubic <- c(solve(t(Theta0123)%*%W%*%Theta0123, t(Theta0123)%*%(W%*%ll)))
            resids_cubic <- ll - c(Theta0123%*%Ahat_cubic)
            sigsqhat_cubic <- c(resids_cubic%*%W%*%resids_cubic) / M
            pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-(d+1)*(d+2)*(d+3)/6)/(d*(d+1)*(d+2)/6), d*(d+1)*(d+2)/6, sum(w>0)-(d+1)*(d+2)*(d+3)/6, lower.tail=FALSE)
        } else {
            cubic_test <- FALSE
        }
        ## test about moments
        if (test=="moments") {
            if (!is.list(null.value[[1]])) { # if the first element of null.value is not a list, it should be a list of length four (the null values for a, b, c, and sigma^2). This corresponds to the case where only a single quadruple is tested. If this is the case, coerce `null.value` into a list of a list of length four to be consistent with the other case where multiple quadruples are tested.
                null.value <- list(null.value)
            }
            teststats <- sapply(null.value,
                function(x) {
                    a_null <- c(x[[1]] + x[[2]]%*%theta_mean + theta_mean %*% x[[3]] %*% theta_mean) # a in the transformed scale
                    b_null <- diag(theta_sd, nrow=length(theta_sd))%*%(x[[2]] + 2*x[[3]]%*%theta_mean) # b in the transformed scale
                    c_null <- diag(theta_sd, nrow=length(theta_sd))%*%x[[3]]%*%diag(theta_sd, nrow=length(theta_sd))
                    err <- ll - c(Theta012%*%c(a_null, b_null, vech(c_null)))
                    .5*M*log(sigsqhat/x[[4]]) - .5*c(err%*%W%*%err)/x[[4]] + M/2
                })
            num.error.size <- 0.01
            pvalout <- pscl(teststats, M, (d^2+3*d+2)/2, num_error_size=num.error.size)
            if (length(pvalout)==0) { # execution of pscl stopped by user input
                stop("Hypothesis tests stopped by user input")
            }
            pval <- pvalout$probs
            num.error.size <- pvalout$numerical_error_size
            if (any(pval < .01)) {
                num.error.size <- 0.001
                pvalout <- pscl(teststats, M, (d^2+3*d+2)/2, num_error_size=num.error.size)
                if (length(pvalout)==0) { # execution of pscl stopped by user input
                    stop("Hypothesis tests stopped by user input")
                }
                pval <- pvalout$probs
                num.error.size <- pvalout$numerical_error_size
            }
            precdigits <- max(-floor(log10(num.error.size)), 1)
            dfout <- lapply(1:length(null.value), function(ii) {
                list(
                    a_null=null.value[[ii]][[1]],
                    b_null=null.value[[ii]][[2]],
                    c_null=null.value[[ii]][[3]],
                    sigma_sq_null=null.value[[ii]][[4]],
                    pvalue=pval[ii]
                )}
            )
            out <- list(meta_model_MLE_for_moments=list(a=ahat_b, b=bhat_b, c=chat_b, sigma_sq=sigsqhat),
                Hypothesis_Tests=dfout,
                pvalue_numerical_error_size=num.error.size
            )
            if (cubic_test) {
                out <- c(out, pval_cubic=pval_cubic)
            }
            return(out)
        }
        ## test about MESLE
        if (test=="MESLE") {
            null.value_n <- lapply(null.value, trans_n) # transform the null values
            U <- t(Theta012)%*%W%*%Theta012
            V <- U[-1,-1] - U[-1,1,drop=FALSE]%*%U[1,-1,drop=FALSE]/U[1,1]
            teststats <- sapply(null.value_n,
                function(x) {
                    theta0mat <- matricize(x)
                    Vq <- solve(cbind(diag(d), 2*theta0mat) %*% solve(V) %*% rbind(diag(d), 2*t(theta0mat)))
                    xi <- t(bhat + 2*theta0mat %*% vech_chat) %*% Vq %*% (bhat + 2*theta0mat %*% vech_chat)
                    (M-(d^2+3*d+2)/2)*xi/(M*d*sigsqhat)
                })
            pval <- pf(teststats, d, M-(d^2+3*d+2)/2, lower.tail=FALSE)
            MESLE_null <- sapply(null.value, identity)
            if (!is.null(attr(MESLE_null,"dim"))) {
                MESLE_null <- t(MESLE_null) # if MESLE_null is a matrix, transpose it
            }
            dfout <- data.frame(
                MESLE_null,
                pvalue=pval
            )
            out <- list(regression_estimates=list(a=ahat_b, b=bhat_b, c=chat_b, sigma_sq=sigsqhat),
                meta_model_MLE_for_MESLE=c(trans_b(MESLEhat)),
                Hypothesis_Tests=dfout
            )
            if (cubic_test) {
                out <- c(out, pval_cubic=pval_cubic)
            }
            return(out)
        }
        ## test on the simulation based surrogate under LAN
        if (test=="parameter") {
            Winv <- diag(1/w, nrow=length(w))
            nobs <- dim(llmat)[1] # number of observations
            slope_at <- apply(theta_n, 2, mean)
            ## estimation of slopes
            if (case=="iid" || K1_est_method=="autocov") {
                ## quadratic regression for each observation piece (each row of simll)
                Ahat_i <- solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%t(llmat)) # each column of Ahat_i is the regression estimate for a single i
                slope <- Ahat_i[bindex,,drop=FALSE]+2*matricize(slope_at)%*%Ahat_i[cindex,,drop=FALSE] # estimated slope, a d X nobs matrix
            } else if (K1_est_method=="batch") {   
                if (is.null(batch_size)) {
                    batch_size <- round(nobs^0.4)
                }
                bsizes <- rep(batch_size, floor(nobs/batch_size))
                if (nobs != batch_size*floor(nobs/batch_size)) {
                    bsizes <- c(bsizes, nobs-batch_size*floor(nobs/batch_size)) # actual batch sizes
                }
                tllmat_b <- sapply(seq(1,nobs,by=batch_size), function(bst) {
                    apply(llmat[bst:min(bst+batch_size-1,nobs),,drop=FALSE], 2, sum)
                })
                Ahat_b <- solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%tllmat_b) # each column of Ahat_b is the regression estimate for a single batch
                slope_b <- Ahat_b[bindex,,drop=FALSE]+2*matricize(slope_at)%*%Ahat_b[cindex,,drop=FALSE] # estimated slope, a d X (num.of.batches) matrix
                if (any(abs(acf(t(slope_b), plot=plot_acf)$acf[2,,,drop=FALSE])>2*sqrt(d*batch_size/nobs))) {
                    warning("The slope estimates at consecutive batches had significant correlation. Consider manualy increasing the batch size for estimation of K1.")
                }
            } else {
                stop("`K1_est_method` should be 'autocov' or 'batch' when `case` is 'stationary'.")
            }
            if (case=="iid") {
                var_slope_vec <- var(t(slope)) # estimate of the variance of the slope of the fitted quadratic
            } else if (case=="stationary" && K1_est_method=="batch") {
                mean_slope_vec <- apply(slope_b, 1, sum)/nobs
                var_slope_vec <- 1/(length(bsizes)-1) * (slope_b - outer(mean_slope_vec, bsizes)) %*% diag(1/bsizes) %*% t(slope_b - outer(mean_slope_vec, bsizes))
            } else if (case=="stationary" && K1_est_method=="autocov") {
                if (is.null(max_lag)) {
                    acfmat <- acf(t(slope),plot=plot_acf)$acf
                    insig_corr <- which(apply(acfmat[-1,,,drop=FALSE], 1, function(x) { all(abs(x) < 2*sqrt(d/nobs)) })) # lags for which the autocorrelation is not significant
                    if (length(insig_corr)==0) {
                        max_lag <- dim(acfmat)[1] - 1
                    } else {
                        max_lag <- min(insig_corr) - 1
                    }
                }
                var_slope_vec <- matrix(0, d, d) # estimate of the variance of the slope of the fitted quadratic
                for (i1 in 1:d) {
                    for (i2 in 1:d) {
                        for (lag in max_lag:-max_lag) {
                            var_slope_vec[i1, i2] <- var_slope_vec[i1, i2] + cov(slope[i1,max(1,1+lag):min(nobs,nobs+lag)], slope[i2,max(1,1-lag):min(nobs,nobs-lag)])
                        }
                    }
                }
            }
            E_condVar_slope <- cbind(0,diag(d),2*matricize(slope_at))%*%solve(t(Theta012)%*%W%*%Theta012, t(cbind(0,diag(d),2*matricize(slope_at)))) * sigsqhat / nobs # an estimate of the expected value of the conditional variance of the estimated slope given Y (expectation taken with respect to Y)
            K1hat <- var_slope_vec - E_condVar_slope
            if (any(eigen(K1hat)$values <= 0)) {
                warning("The estimate of K1 is not positive definite. The result of the hypothesis test may not be reliable.")
            }
            null.value_n <- lapply(null.value, trans_n) # transform the null values
            C <- cbind(-1, diag(rep(1,M-1)))
            Theta12 <- Theta012[,-1]
            Ctheta <- C%*%theta_n
            Q_1 <- solve(C%*%Winv%*%t(C) + nobs/sigsqhat*Ctheta%*%K1hat%*%t(Ctheta)) 
            svdQ_1 <- svd(Q_1)
            sqrtQ_1 <- svdQ_1$u %*% diag(sqrt(svdQ_1$d), nrow=M-1) %*% t(svdQ_1$v)
            invsqrtQ_1 <- svdQ_1$v %*% diag(sqrt(1/svdQ_1$d), nrow=M-1) %*% t(svdQ_1$u)
            R_1 <- sqrtQ_1%*%C%*%Theta12
            estEq_2 <- c(solve(t(R_1)%*%R_1, t(R_1)%*%(sqrtQ_1%*%(C%*%ll)))) # estimating equation for K2 and theta_star. (thetastarhat // -I/2) * n * vech(K2hat) = (R_1^T R_1)^{-1} R_1^T (Q_1^{1/2} C lS)
            vech_K2hat <- -2*estEq_2[(d+1):((d^2+3*d)/2)]/nobs # second stage estimate of vech(K2)
            K2hat <- unvech(vech_K2hat)
            thetastarhat <- solve(K2hat,estEq_2[1:d])/nobs # maximum meta model likelihood estimate for theta_star (simulation based surrogate)
            sigsqhat_lan <- 1/(M-1)*sum((sqrtQ_1%*%(C%*%ll) - R_1%*%estEq_2)^2)
            teststats <- sapply(null.value_n,
                function(x) {
                    desgmat <- rbind(-2*matricize(x), diag((d^2+d)/2))
                    G <- R_1%*%desgmat%*%solve(t(desgmat)%*%t(R_1)%*%R_1%*%desgmat, t(desgmat)%*%t(R_1))
                    (M-(d^2+3*d+2)/2)/d*(sum(((diag(M-1)-G)%*%(sqrtQ_1%*%(C%*%ll)))^2)/(M-1)/sigsqhat_lan - 1)
                })
            if (MCcorrection=="none") {
                pval <- pf(teststats, d, M-(d^2+3*d+2)/2, lower.tail=FALSE)
            } else if (MCcorrection=="Wishart") {
                MCsize <- 500 # Monte Carlo sample size
                QhCllMC <- outer(c(R_1%*%estEq_2), rep(1,MCsize)) + sqrt(sigsqhat_lan)*matrix(rnorm((M-1)*MCsize), M-1, MCsize) # Q^{1/2}C l^S
                llMC <- rbind(0, invsqrtQ_1%*%(QhCllMC)) # Monte Carlo draws for the simulation log likelihood vector
                resMC <- llMC - Theta012%*%solve(t(Theta012)%*%W%*%Theta012, t(Theta012)%*%W%*%llMC)
                sigsqhatMC <- 1/M*apply(resMC, 2, function(v) sum(v*v))
                if (case=="iid") {
                    df <- nobs-1
                } else if (K1_est_method=="autocov") {
                    df <- nobs/(2*max_lag+1)-1
                } else if (K1_est_method=="batch") {
                    df <- nobs/batch_size-1
                }
                tau1MC <- rWishart(MCsize, df=df, Sigma=var_slope_vec/df)
                tau2pre <- cbind(0,diag(d),2*matricize(slope_at))%*%solve(t(Theta012)%*%W%*%Theta012, t(cbind(0,diag(d),2*matricize(slope_at)))) / nobs
                #sqrtQ_1MC <- array(NA, dim=c(M-1, M-1, MCsize))
                #R_1MC <- array(NA, dim=c(M-1, d*(d+3)/2, MCsize))
                #tRRMC <- array(NA, dim=c(d*(d+3)/2, d*(d+3)/2, MCsize))
                #sigsqhat_lanMC <- rep(NA, MCsize)
                desgmat <- rbind(-2*matricize(thetastarhat), diag((d^2+d)/2))
                teststatsMC <- sapply(1:MCsize, function(i_mc) {
                    tau2MC <- tau2pre * sigsqhatMC[i_mc]
                    K1hatMC <- tau1MC[,,i_mc] - tau2MC
                    Q_1MC <- solve(C%*%Winv%*%t(C) + nobs/sigsqhatMC[i_mc]*Ctheta%*%K1hatMC%*%t(Ctheta))
                    svdQ_1MC <- svd(Q_1MC)
                    sqrtQ_1MC <- svdQ_1MC$u %*% diag(sqrt(svdQ_1MC$d), nrow=M-1) %*% t(svdQ_1MC$v)
                    R_1MC <- sqrtQ_1MC%*%(C%*%Theta12)
                    sigsqhat_lanMC <- 1/(M-1)* sum(((diag(M-1) - R_1MC%*%solve(t(R_1MC)%*%R_1MC, t(R_1MC)))%*%QhCllMC[,i_mc])^2)
                    DMC <- R_1MC%*%desgmat
                    GMC <- DMC%*%solve(t(DMC)%*%DMC,t(DMC))
                    (M-(d^2+3*d+2)/2)/d*(sum(((diag(M-1)-GMC)%*%QhCllMC[,i_mc])^2)/(M-1)/sigsqhat_lanMC - 1)
                })
                pval <- sapply(1:length(null.value_n), function(i_nv) {
                    nv <- null.value_n[[i_nv]]
                    mean(teststats[i_nv]<=teststatsMC)
                })
            }
            parameter_null <- sapply(null.value, identity)
            if (!is.null(attr(parameter_null,"dim"))) {
                parameter_null <- t(parameter_null) # if surrogate_null is a matrix, transpose it
            }
            dfout <- data.frame(
                parameter_null,
                pvalue=pval
            )
            out <- list(regression_estimates=list(a=ahat_b, b=bhat_b, c=chat_b, sigma_sq=sigsqhat),
                meta_model_MLE_for_parameter=c(trans_b(thetastarhat)),
                K1=diag(1/theta_sd, nrow=length(theta_sd))%*%K1hat%*%diag(1/theta_sd, nrow=length(theta_sd)),
                K2=diag(1/theta_sd, nrow=length(theta_sd))%*%K2hat%*%diag(1/theta_sd, nrow=length(theta_sd)),
                error_variance=sigsqhat_lan,
                Hypothesis_Tests=dfout
            )
            if (case=="stationary") {
                out <- c(out, max_lag=max_lag)
            }
            if (cubic_test) {
                out <- c(out, pval_cubic=pval_cubic)
            }
            return(out)
        }
    }
}


