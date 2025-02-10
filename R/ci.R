#' @export
ci <- function(simll, ...) {
    UseMethod("ci")
}

#' Confidence interval for scalar parameter constructed using simulated log likelihoods
#'
#' `ci` constructs confidence intervals for a scalar (one-dimensional) parameter using simulated log likelihoods. See Park (2025) for more information.
#'
#' @name ci
#' @param simll A class 'simll' object, containing simulated log likelihoods, the parameter values at which simulations are made, and the weights for those simulations for regression (optional). See help(simll).
#' @param level A numeric vector of confidence levels.
#' @param ci A character string indicating the quantity for which a confidence interval is to be constructed. Either "MESLE" or "parameter". See Details.
#' @param case When `ci` is "parameter", `case` is either "iid" or "stationary" (default). `case` = "iid" means that the observations are iid, and `case` = "stationary" means that the observations form a stationary sequence. The `case` argument affects how the variance of the slope of the mean function (=K_1 in Park (2025)) is estimated.
#' @param weights An optional argument. The un-normalized weights of the simulated log likelihoods for regression. A numeric vector of length equal to the 'params' attribute of the 'simll' object. See Details below.
#' @param autoScaling logical. If TRUE, simulation points at which the third order term is statistically significant in the cubic approximation to the simulated log-likelihooods have discounted weights for metamodel fitting. The weights of the points relatively far from the estimated MESLE are more heavily discounted. These weight discount factors are multiplied to the originally given weights for parameter estimation. See Park (2025) for more details. If `autoScaling` is FALSE, the weight discount step is skipped. Defaults to FALSE.
#' @param K1_est_method Either "autocov" or "batch"
#' @param batch_size Numeric
#' @param max_lag When `test` is "parameter" and `case` is "stationary", the value of `max_lag` gives the truncation point for lagged autocovariance when estimating K1 as a sum of lagged autocovariances of estimates slopes. If not supplied, default is the maximum lag for which the lagged autocorrelation has absolute value greater than 4/sqrt(nobs), where the lagged autocorrelation is found up to lag `10*log10(nobs)`. Here `nobs` is the number of observations.
#' @param plot_acf Logical.  When `test` is "parameter" and `case` is "stationary", If `plot_acf` is TRUE, the autocorrelation plot of the estimated slopes of the quadratic fit to the simulated log likelihoods is shown.
#' @param ... Other optional arguments, not currently used.
#'
#' @details
#' This is a generic function, taking a class 'simll' object as the first argument.
#' Confidence intervals are constructed under a normal, locally quadratic meta model where the simulated log likelihoods given in the 'simll' object are normally distributed.
#'
#' When 'level' has length greater than one, a confidence interval is constructed for each value in the vector.
#'
#' Quadratic regression for the simulated log likelihoods is carried out to construct confidence intervals, where the x-axis values are the 'params' values of the 'simll' object and the y-axis values are the corresponding simulated log likelihoods.
#' In the case where 'ci' = "parameter", inference on the simulation based surrogate will be carried out under the local asymptotic normality for simulated log likelihoods (see Park (2025) for more information.)
#' The default value of 'ci' is "parameter".
#'
#' If 'ci' = "MESLE", confidence intervals are constructed for the maximum expected simulation likelihood estimate given the observed data.
#'
#' When quadratic regression is carried out, the weights for the simulation based likelihood estimates can be specified. The length of 'weights' should be equal to that of the 'params' attribute of the 'simll', which is equal to the number of rows in the simulated log likelihood matrix in the 'simll' object. It is important to note that the weights are not normalized (i.e., not sum to one). Multiplying all weights by the same constant changes the estimation outputs. If not supplied, the 'weights' attribute of the 'simll' object is used. If neither is supplied, 'weights' defaults to the vector of all ones.
#'
#' @return A list consisting of the followings are returned.
#' \itemize{
#' \item{regression_estimates: point estimates for the meta model parameters, a, b, c, and sigma^2.}
#' \item{meta_model_MLE_for_*: point estimate for the quantity for which confidence intervals are constructed under a normal meta model}
#' \item{confidence_interval: a data frame of the lower and upper bounds of the confidence intervals and the corresponding confidence levels. Note that in some unfortunate cases (especially if the quadratic coefficient of the estimated quadratic fit of the log likelihood estimates is close to zero or nonnegative), the confidence interval may be inverted, meaning that it is of the form (-infty, bound1) U (bound2, infty). This case can happen if the signal-to-noise ratio in simulated log likelihoods is too small. The inverted confidence interval will be indicated by the additional column "inverted" in the data frame taking values of 0 or 1.}
#' \item{max_lag: if `test`="parameter" and `case`="stationary", the maximum lag for computing the autocovariance in estimating K1 is shown.}
#' \item{pval_cubic: The p-value of the test about whether the cubic term in the cubic polynomial regression is significant. If so, the constructed confidence interval may be biased.}
#' }
#'
#' @references Park, J. (2025). Scalable simulation based inference for implicitly defined models using a metamodel for log-likelihood estimator <https://doi.org/10.48550/arxiv.2311.09446>
#' @export
ci.simll <- function(simll, level, ci=NULL, case=NULL, weights=NULL, autoScaling=FALSE, K1_est_method="batch", batch_size=NULL, max_lag=NULL, plot_acf=FALSE, ...) {
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
    if (!is.null(attr(attr(simll, "params"), "dim")) && dim(attr(simll, "params"))[2] != 1) {
        stop(
            "Simulation based confidence interval is only constructed for scalar parameters.",
            call. = FALSE
        )
    }
    d <- 1
    if (ci=="parameter" && is.null(case)) {
        case <- "stationary"
    }
    if (ci=="parameter") {
        match.arg(case, c("iid", "stationary"))
    }
    if (!is.numeric(level)) {
        stop(
            "'level' should be a numeric vector (length >= 1).",
            call. = FALSE
        )
    }
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
                    "When the 'simll' object has 'weights' attribute, the length of 'weights' should be the same as the number of columns in the simulated log likelihood matrix in 'simll'.",
                    call. = FALSE
                )
            }
            w <- attr(simll, "weights")
        } else {
            w <- rep(1, dim(simll)[2])
        }
    }
    ## weighted quadratic regression
    theta <- cbind(attr(simll, "params"))
    llmat <- unclass(simll)
    ll <- apply(llmat, 2, sum)
    M <- length(ll)
    theta012 <- cbind(1, theta, theta^2)
    Ahat <- c(solve(t(theta012)%*%(outer(w,rep(1,3))*theta012), t(theta012)%*%(w*ll))) # Ahat=(ahat,bhat,chat)
    resids <- ll - c(theta012%*%Ahat)
    sigsqhat <- c(resids%*%(w*resids)) / M
    MESLEhat <- unname(-Ahat[2]/(2*Ahat[3])) # estimate of the MESLE
    ## cubic polynomial test below
    theta0123 <- cbind(1, theta, theta^2, theta^3) # cubic regression to test whether the cubic coefficient = 0
    Ahat_cubic <- c(solve(t(theta0123)%*%(outer(w,rep(1,4))*theta0123), t(theta0123)%*%(w*ll)))
    resids_cubic <- ll - c(theta0123%*%Ahat_cubic)
    sigsqhat_cubic <- c(resids_cubic%*%(w*resids_cubic)) / M
    pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-4), 1, sum(w>0)-4, lower.tail=FALSE)
    if (autoScaling) {
        if (M <= (d+1)*(d+2)*(d+3)/6) { # carry out cubic test if this condition is met
            stop("The number of simulations is not large enough to carry out cubic polynomial fitting (should be greater than (d+1)*(d+2)*(d+3)/6)")
        }
        qa <- function(x) { sum(c(1, x, x^2)*Ahat) }
        logwpen <- function(point) { -(qa(MESLEhat)-qa(point))/refgap } # penalizaing weight (weight discount factor)
        vec3 <- function(val) {
            val^3
        }
        theta0123 <- cbind(theta012, theta^3) # design matrix for cubic regression to test whether the cubic coefficient = 0
        refgap <- Inf # reference value for the gap qa(MESLEhat)-qa(theta) where qa is the quadratic approximation
        exit_upon_sufficient_ESS <- FALSE
        repeat{
            ## Weight points appropriately to make the third order term insignificant
            wadj <- w * exp(apply(theta, 1, logwpen)) # adjusted weights
            WadjTheta012 <- outer(wadj,rep(1,3))*theta012
            Ahat <- c(solve(t(theta012)%*%WadjTheta012, t(theta012)%*%(wadj*ll)))
            resids <- ll - c(theta012%*%Ahat)
            sigsqhat <- c(resids%*%(wadj*resids)) / M
            MESLEhat <- unname(-Ahat[2]/(2*Ahat[3])) 
            Ahat_cubic <- c(solve(t(theta0123)%*%(outer(wadj,rep(1,4))*theta0123), t(theta0123)%*%(wadj*ll)))
            resids_cubic <- ll - c(theta0123%*%Ahat_cubic)
            sigsqhat_cubic <- c(resids_cubic%*%(wadj*resids_cubic)) / M
            pval_cubic <- pf((sigsqhat-sigsqhat_cubic)/sigsqhat_cubic*(sum(w>0)-4), 1, sum(w>0)-4, lower.tail=FALSE)
            ESS <- sum(wadj)^2/sum(wadj^2) # effective sample size (ESS)
            if (ESS <= (d+1)*(d+2)*(d+3)/6) { # if the ESS is too small, increase refgap
                exit_upon_sufficient_ESS <- TRUE # break from loop as soon as the ESS is large enough
                refgap <- refgap * 1.5
                next
            }
            if (exit_upon_sufficient_ESS) {
                break
            }
            if (pval_cubic < .01) {
                if (refgap==Inf) {
                    refgap <- qa(MESLEhat) - min(apply(theta, 1, qa))
                } else {
                    refgap <- refgap / 1.8
                }
            } else {
                break
            }
        }
        w <- wadj
    }
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
        out <- list(regression_estimates=list(a=Ahat[1], b=Ahat[2], c=Ahat[3], sigma_sq=sigsqhat),
            meta_model_MLE_for_MESLE=c(MESLE=MESLEhat),
            confidence_interval=t(lub)
        )
        if (autoScaling) {
            out[["updated_weights"]] <- w
        }
        return(out)
    }
    ## ci for the simulation surrogate under LAN
    if (ci=="parameter") {
        nobs <- dim(llmat)[1] # number of observations
        slope_at <- mean(theta)
        ## estimation of slopes
        if (case=="iid" || K1_est_method=="autocov") {
            ## quadratic regression for each observation piece (each row of simll)
            Ahat_i <- solve(t(theta012)%*%(outer(w,rep(1,3))*theta012), t(outer(w,rep(1,3))*theta012)%*%t(llmat)) # each column of Ahat_i is the regression estimate for a single i
            slope <- Ahat_i[2,]+2*slope_at*Ahat_i[3,] # estimated slope, a d X nobs matrix
        } else if (K1_est_method=="batch") {
            if (is.null(batch_size)) {
                batch_size <- round(nobs^0.4)
            }
            tllmat_b <- sapply(seq(1,nobs,by=batch_size), function(bst) {
                apply(llmat[bst:min(bst+batch_size-1,nobs),,drop=FALSE], 2, sum)
            })
            Ahat_b <- solve(t(theta012)%*%(outer(w,rep(1,3))*theta012), t(outer(w,rep(1,3))*theta012)%*%tllmat_b) # each column of Ahat_b is the regression estimate for a single batch
            slope_b <- Ahat_b[2,]+2*slope_at*Ahat_b[3,] # estimated slope, a d X (num.of.batches) matrix
            if (any(abs(acf(slope_b, plot=plot_acf)$acf[2,,,drop=FALSE])>2*sqrt(batch_size/nobs))) {
                warning("The slope estimates at consecutive batches had significant correlation. Consider manualy increasing the batch size for estimation of K1.")
            }
        } else {
            stop("`K1_est_method` should be 'autocov' or 'batch' when `case` is 'stationary'.")
        }
        if (case=="iid") {
            var_slope <- var(slope) # estimate of the variance of the slope of the fitted quadratic
        } else if (case=="stationary" && K1_est_method=="batch") {
            var_slope <- var(slope_b) / batch_size
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
            var_slope <- 0 # estimate of the variance of the slope of the fitted quadratic
            for (lag in max_lag:-max_lag) {
                var_slope <- var_slope + cov(slope[max(1,1+lag):min(nobs,nobs+lag)], slope[max(1,1-lag):min(nobs,nobs-lag)])
            }
        }
        E_condVar_slope <- c(cbind(0,1,2*slope_at)%*%solve(t(theta012)%*%(outer(w,rep(1,3))*theta012), rbind(0,1,2*slope_at))) * sigsqhat / nobs # an estimate of the expected value of the conditional variance of the estimated slope given Y (expectation taken with respect to Y)
        K1hat <- var_slope - E_condVar_slope
        if (K1hat <= 0) {
            warning("The estimate of K1 is not positive definite. The constructed confidence interval will not be reliable.")
        }
        theta12 <- theta012[,-1]
        Wbartheta12 <- outer(w,rep(1,2))*theta12 - 1/sum(w)*cbind(w)%*%(w%*%theta12)
        Wbartheta <- Wbartheta12[,1]
        Wbarll <- w*ll - 1/sum(w)*w*sum(w*ll)
        ## P_1 = Wbar - Wbar %*% theta %*% solve(sigsqhat/K1hat/nobs + t(theta) %*% Wbar %*% theta, t(theta) %*% Wbar)
        Ptheta12 <- Wbartheta12 - Wbartheta %*% solve(sigsqhat/nobs*solve(K1hat) + t(theta) %*% Wbartheta, t(theta)%*%Wbartheta12)
        Pll <- Wbarll - Wbartheta %*% solve(sigsqhat/nobs*solve(K1hat) + t(theta) %*% Wbartheta, t(theta)%*%Wbarll)
        estEq_2 <- solve(t(theta12)%*%Ptheta12, t(Ptheta12)%*%ll) # estimating equation for K2 and theta_star. (thetastarhat // -I/2) * n * vech(K2hat) = (R_1^T R_1)^{-1} R_1^T (Q_1^{1/2} C lS)
        K2hat <- -2*estEq_2[2]/nobs # second stage estimate of K2
        thetastarhat <- estEq_2[1]/(nobs*K2hat) # maximum meta model likelihood estimate for theta_star (simulation based surrogate)
        sigsqhat_lan <- 1/(M-1)*c(t(ll-theta12%*%estEq_2)%*%(Pll-Ptheta12%*%estEq_2))
        theta12Ptheta12 <- t(theta12)%*%Ptheta12
        rho11 <- theta12Ptheta12[1,1]
        rho12 <- theta12Ptheta12[1,2]
        rho22 <- theta12Ptheta12[2,2]
        lub <- sapply(1:length(level), function(i) {
            lvl <- level[i]
            q <- qf(lvl, 1, M-3)
            zeta0 <- c(t(ll)%*%Pll) - (M-1)*sigsqhat_lan*(q/(M-3)+1)
            zeta12 <- t(theta12)%*%Pll
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
        out <- list(regression_estimates=list(a=Ahat[1], b=Ahat[2], c=Ahat[3], sigma_sq=sigsqhat),
            meta_model_MLE_for_parameter=c(parameter=thetastarhat, K1=K1hat, K2=K2hat, error_variance=sigsqhat_lan),
            confidence_interval=t(lub)
        )
        if (case=="stationary") {
            out <- c(out, max_lag=max_lag)
        }
        if (autoScaling) {
            out[["updated_weights"]] <- w
        }
        out <- c(out, pval_cubic=pval_cubic)
        return(out)
    }
}
# When `case`="stationary", a test of stationarity can be carried out, and a warning should be given when the test is positive.

