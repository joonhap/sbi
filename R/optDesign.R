#' @export
optDesign <- function(simll, ...) {
    UseMethod("optDesign")
}

#' Find the next optimal design point for simulation-based inference
#'
#' `optDesign` finds the next design point at which simulation should be carried out for approximately best efficiency in a metamodel-based inference. See Park (2025) for more details on this method. It takes a class `simll` object.
#'
#' @name optDesign
#' @param s A class `simll` object, containing simulation log likelihoods, the parameter values at which simulations are made, and the weights for those simulations for regression (optional). See help(simll).
#' @param init (optional) An initial parameter vector at which a search for optimal point starts.
#' @param weight (optional) A positive real number indicating the user-assigned weight for the new design point. The default value is 1. This value should be chosen relative to the weights in the provided simll object.
#' @param refgap (optional) A positive real number that determines the weight discount factor for the significance of the third order term in Taylor approximation. The weight of a point `theta` is discounted by a factor of exp(-(qa(theta)-qa(MESLEhat))/refgap), where MESLEhat is the estimated MESLE and qa is the quadratic approximation to the simulated log-likelihoods (before weight discount). Defaults to `refgap=Inf`.
#' @param ... Other optional arguments, not currently used.
#'
#' @details
#' This is a generic function, taking a class `simll` object as the first argument.
#' Parameter inference for implicitly defined simulation models can be carried out under a metamodel for the distribution of the log-likelihood estimator.
#' See function `ht` for hypothesis testing and `ci` for confidence interval construction for a one-dimensional parameter.
#' This function `optDesign` finds the next point at which a simulation is to be carried out such that the variance of the parameter estimate is reduced approximately the most.
#' Points far from the estimated MESLE have discounted weights such that the third order term in the Taylor approximation is not statistically significant.
#' The weight discount factor the a point at theta is given by exp(-(qa(theta)-qa(MESLEhat))^2/g^2), where qa is the quadratic approximation, MESLEhat is the estimated MLE, and g is a scaling parameter such that the resulting test for the third order approximation term is not statistically significant.
#' These discount factors are multipled to the original `weights` given to the simulation points specified in the `s` object.
#'
#' @return A matrix of parameter values at which next simulations are to be carried out for approximately best efficiency.
#' Each row gives a parameter vector.
#'
#' @references Park, J. (2025). Scalable simulation-based inference for implicitly defined models using a metamodel for log-likelihood estimator <https://doi.org/10.48550/arxiv.2311.09446>
#' @export
optDesign.simll <- function(s, init=NULL, weight=1, refgap=Inf, ...) {
    #validate_simll(s)
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
    ## weighted quadratic regression
    if (is.null(attr(attr(s, "params"), "dim"))) {
        d <- 1
    } else {
        d <- dim(attr(s, "params"))[2]
    }
    if (!is.null(attr(s, "weights"))) {
        if (!is.numeric(attr(s, "weights"))) {
            stop("When the `simll` object `s` has `weights` attribute, it has to be a numeric vector.")
        }
        if (dim(s)[2] != length(attr(s, "weights"))) {
            stop("When the `simll` object `s` has `weights` attribute, the length of `weights` should be the same as the number of rows in the simulated log likelihood matrix in `s`.")
        }
        w <- attr(s, "weights")
    } else {
        w <- rep(1, dim(s)[2])
    }
    theta <- cbind(attr(s, "params")) # coerce into a matrix
    theta_mean <- apply(theta, 2, mean)
    theta_sd <- apply(theta, 2, sd)
    trans_n <- function(vec) { (vec-theta_mean)/theta_sd } # normalize by centering and scaling
    trans_b <- function(vec) { vec*theta_sd + theta_mean } # transform back to the original scale
    theta_n <- apply(theta, 1, trans_n) |> rbind() |> t() # apply trans_n rowwise (result:M-by-d matrix)
    llmat <- unclass(s)
    ll <- apply(llmat, 2, sum)
    M <- length(ll)
    Theta012 <- t(apply(theta_n, 1, vec012))
    dim012 <- 1 + d + (d^2+d)/2
    ## first stage approximation of MESLEhat
    WTheta012 <- outer(w,rep(1,dim012))*Theta012
    Ahat <- c(solve(t(Theta012)%*%WTheta012, t(Theta012)%*%(w*ll)))
    bhat <- Ahat[2:(d+1)]
    vech_chat <- Ahat[(d+2):((d^2+3*d+2)/2)]
    chat <- unvech(vech_chat)
    resids <- ll - c(Theta012%*%Ahat)
    sigsqhat <- c(resids%*%(w*resids)) / M
    MESLEhat <- unname(-solve(chat,bhat)/2)
    qa <- function(x) { sum(vec012(x)*Ahat) }
    logwpen <- function(point) { -(qa(MESLEhat)-qa(point))/refgap } # penalizaing weight (weight discount factor)
    ## cubic test
    if (M <= (d+1)*(d+2)*(d+3)/6) { # carry out cubic test if this condition is met
        stop("The number of simulations is not large enough to carry out cubic polynomial fitting (should be greater than (d+1)*(d+2)*(d+3)/6)")
    }
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
    dim0123 <- dim(Theta0123)[2]
    refgap_init <- refgap # initial value of refgap
    exit_upon_sufficient_ESS <- FALSE
    repeat{
        ## Weight points appropriately to make the third order term insignificant
        wadj <- w * exp(apply(theta_n, 1, logwpen)) # adjusted weights
        WadjTheta012 <- outer(wadj,rep(1,dim012))*Theta012
        Ahat <- c(solve(t(Theta012)%*%WadjTheta012, t(Theta012)%*%(wadj*ll)))
        bhat <- Ahat[2:(d+1)]
        vech_chat <- Ahat[(d+2):((d^2+3*d+2)/2)]
        chat <- unvech(vech_chat)
        resids <- ll - c(Theta012%*%Ahat)
        sigsqhat <- c(resids%*%(wadj*resids)) / M
        MESLEhat <- unname(-solve(chat,bhat)/2)
        Ahat_cubic <- c(solve(t(Theta0123)%*%(outer(wadj,rep(1,dim0123))*Theta0123), t(Theta0123)%*%(wadj*ll)))
        resids_cubic <- ll - c(Theta0123%*%Ahat_cubic)
        sigsqhat_cubic <- c(resids_cubic%*%(wadj*resids_cubic)) / M
        sigratio <- (sigsqhat-sigsqhat_cubic)/sigsqhat_cubic
        multiplier <- (sum(w>0)-(d+1)*(d+2)*(d+3)/6)/(d*(d+1)*(d+2)/6)
        fstat <- sigratio * multiplier
        pval_cubic <- pf(fstat, d*(d+1)*(d+2)/6, sum(w>0)-(d+1)*(d+2)*(d+3)/6, lower.tail=FALSE)
        ESS <- sum(wadj)^2/sum(wadj^2) # effective sample size (ESS)
        if (ESS <= (d+1)*(d+2)*(d+3)/6) { # if the ESS is too small, increase refgap
            ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
            exit_upon_sufficient_ESS <- TRUE # break from loop as soon as the ESS is large enough
            refgap <- refgap * 1.5
            next
        }
        if (exit_upon_sufficient_ESS) {
            ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
            ##cat("break, sufficient ESS reached.\n")
            break
        }
        if (pval_cubic < .01) {
            ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
            if (refgap==Inf) {
                refgap <- qa(MESLEhat) - min(apply(theta_n, 1, qa))
            } else {
                refgap <- refgap / 1.8
            }
        } else if (pval_cubic > .3) {
            if (refgap >= 10*refgap_init) { # if refgap has been increased a lot already, stop. Note: in order to account for the case where pval_cubic > .3 even with refgap = Inf, the comparison should be ">=" rather than ">".
                ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
                ##cat("break, maximally increased refgap.\n")
                break
            }
            ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
            refgap <- refgap * 1.3
        } else {
            ##cat("refgap:",refgap,"pval_cubic:",pval_cubic,"sigsqhat:",sigsqhat,"sigsqhat_cubic:",sigsqhat_cubic,"sigratio:",sigratio,"multiplier:",multiplier,"fstat:",fstat,"ESS:",ESS,"\n")
            ##cat("break, p-val suitable.\n")
            break
        }
    }
    ## Compute the gradient of the variance of MESLEhat with respect to new design point
    chatinv <- solve(chat)
    pMpchat <- matrix(0, d, (d^2+d)/2) # partial MESLEhat / partial vech(chat)
    n <- 0
    for (i in 1:d) {
        l <- d+1-i
        pMpchat[,(n+1):(n+l)] <- -chatinv[,i:d]*MESLEhat[i]
        n <- n+l
    }
    pMpAhat <- cbind(0, -1/2*chatinv, pMpchat) # partial MESLEhat / partial Ahat
    ## partial MESLEhat / partial ahat = 0
    ## partial MESLEhat / partial bhat = -1/2*solve(chat)
    ## partial MESLEhat / partial vech(chat) = (-solve(chat)*MESLEhat[1], -solve(chat)[,2:d]*MESLEhat[2], ..., -solve(chat)[,d]*MESLEhat[d])
    CubicCoef <- function(i,j,k) { # Supposing that the third order term in the cubic approximation is given by Sum_{j,k,l} e_{j,k,l}*(theta_j-thetahat_j)*(theta_k-thetahat_k)*(theta_l-thetahat_l), this function gives the value of e_{i,j,k}
        sm <- min(i,j,k) # smallest
        lg <- max(i,j,k) # largest
        mi <- i+j+k-sm-lg # middle number
        entry <- Ahat_cubic[(d+1)*(d+2)/2 + lg*(lg-1)*(lg+1)/6+mi*(mi-1)/2+sm]
        nuniq <- length(unique(c(i,j,k))) # number of distinct (unique) values among i,j,k
        if (nuniq==1) { return(entry) }
        if (nuniq==2) { return(entry/3) } # Ahat_cubic entry is (e.g.) e_{112}+e_{121}+e_{211}
        if (nuniq==3) { return(entry/6) } # Ahat_cubic entry is the sum of six permutations of indices
    } ## TODO: remove this CubicCoef function, no longer necessary
    Dlogwpen <- function(point) { # derivative of logwpen
        (bhat + 2*c(chat%*%point))/refgap
    }
    pVinvpPti <- function(point, index) { # partial Var(Ahat)^{-1} / partial point[index], multiplied by sigma^2
        ei <- rep(0, d); ei[index] <- 1
        pPtsqpPti <- matrix(0, d, d) # partial point^2 / partial point[index], where point^2 is a d-by-d matrix (see Park(2025) for definition)
        pPtsqpPti[,index] <- 2*point; pPtsqpPti[index,] <- 2*point
        pPt012pPti <- c(0,ei,vech(pPtsqpPti)) # partial point^{0:2} / partial point[index]
        pt012 <- vec012(point)
        logwpenpt <- logwpen(point)
        out <- exp(logwpenpt) * (Dlogwpen(point)*outer(pt012,pt012) + outer(pPt012pPti,pt012)+outer(pt012,pPt012pPti))
        out * weight
    }
    VinvAhat <- function(point) {
        t(Theta012)%*%WadjTheta012 + exp(logwpen(point))*outer(vec012(point),vec012(point)) # inverse of (variance of Ahat / sigma^2)
    }
    pSTVpPti <- function(point, index) { # partial STV(MESLEhat) / partial point[index], where STV = trace(-chat^{-1}%*%Var(MESLEhat)) is the scaled total variation of MESLEhat
        V_pMpAhatT <- solve(VinvAhat(point), t(pMpAhat)) # Var(Ahat) %*% (partial MESLEhat / partial Ahat)^T (divided by sigma^2)
        sum(diag(solve(chat, t(V_pMpAhatT) %*% pVinvpPti(point, index) %*% V_pMpAhatT)))
    }
    pSTVpPt <- function(point) { # partial log(STV(MESLEhat)) / partial point
        sapply(1:d, function(i) { pSTVpPti(point, i) })
    }
    STV <- function(point) { # STV(MESLEhat) when a new simulation is conducted at `point`
        -sum(diag(solve(chat, pMpAhat %*% solve(VinvAhat(point), t(pMpAhat)))))
    }
    logSTV <- function(point) { log(STV(point)) }
    plogSTVpPt <- function(point) { # partial log(STV(MESLEhat)) / partial point
        sapply(1:d, function(i) { pSTVpPti(point, i) }) / STV(point)
    }
    ## original
    pTVpPti <- function(point, index) { # partial TV(MESLEhat) / partial point[index], where TV = trace(Var(MESLEhat)) is the total variation of MESLEhat
        Vinv_pMpAhatT <- solve(VarAhat(point), t(pMpAhat)) # Var(Ahat) %*% (partial MESLEhat / partial Ahat)^T
        pVpPti_eval <- pVpPti(point, index)
        out <- 0
        for (i in 1:d) {
            out <- out + c(Vinv_pMpAhatT[,i] %*% pVpPti_eval %*% Vinv_pMpAhatT[,i])
        }
        -out
    }
    TV <- function(point) { # TV(MESLEhat) when a new simulation is conducted at `point`
        VarAhatpt <- VarAhat(point)
        out <- 0
        for (i in 1:d) {
            out <- out + c(pMpAhat[i,] %*% solve(VarAhatpt, pMpAhat[i,]))
        }
        out
    }
    logTV <- function(point) { log(TV(point)) }
    plogTVpPt <- function(point) { # partial log(TV(MESLEhat)) / partial point
        sapply(1:d, function(i) { pTVpPti(point, i) }) / TV(point)
    }
    ## gradient descent search
    if (is.null(init)) {
        init_n <- MESLEhat
    } else if (!is.numeric(init) || length(init) != d) {
        init_n <- MESLEhat
        message("The `init` should be a numeric vector of the same length as the parameter vector. Defaults to the estimated MESLE.")
    } else {
        init_n <- trans_n(init)
    }
    if (logwpen(init_n) < log(1e-4)) { # if wpen is too small, move the initial point closer to the MESLEhat such that the issue of too small a gradient can be avoided.
        init_n <- MESLEhat + (init_n - MESLEhat) * (-logwpen(init_n))^(-1/2)
    }

    opt <- optim(trans_n(init_n), fn=logSTV, gr=plogSTVpPt, method="BFGS")

    if (!opt$convergence%in%c(0,1)) {
        print(opt)
        stop("The optimization procedure did not end properly.")
    }
    if (opt$convergence==1) {
        message("The optimization procedure did not end within the max number of iterations.")
    }

    list(par=trans_b(opt$par), logSTV=opt$value, w=exp(logwpen(opt$par)), Wadj=wadj, refgap=refgap)
}

## TODO: add explanations for the return values
## TODO: uncomment validate_simll
