#' Stochastic volatility model
#'
#' `SV_pomp` defines a stochastic volatility model used as an example in the package vignette. The model is defined using the `R` package `pomp`. 
#' @name SV_pomp
#' @param t0: time zero where the rinit is called
#' @param times: observation times
#' @param kappa: a parameter in the model (between 0 and 1, close to 1 means the volatility has high autocorrelation)
#' @param tau: a parameter in the model (positive, the standard deviation of volatility at a given time point)
#' @param sim_seed: random seed for data generation
#' 
#' @details
#' The model is described as follows.
#' 
#' Observations: r_t = exp(s_t)*W_t where W_t~t_5 (t-distribution with 5 degrees of freedom)
#' 
#' Latent stochastic volatility: s_t = kappa*s_t + tau*sqrt(1-kappa^2)*V_t,  V_t~N(0,1),  for t>1
#'        s_1 = tau*V_1,  V_1~N(0,1) for t=1.
#' @export
SV_pomp <- function(t0, times, kappa, gamma, tau, sim_seed=NULL) {
    require(pomp)

    rinit <- Csnippet("
        s = gamma + tau*rnorm(0,1);
    ")

    rproc <- Csnippet("
        s = gamma + kappa*(s-gamma) + tau*sqrt(1-kappa*kappa)*rnorm(0,1);
    ")

    rmeas <- Csnippet("
        r = exp(s)*rt(5);
    ")
    
    dmeas <- Csnippet("
        double logdmeas = dt(r*exp(-s), 5, 1) - s;
        lik = (give_log) ? logdmeas : exp(logdmeas);
    ")

    if (!is.null(sim_seed)) {
        set.seed(sim_seed)
    }
    
    simulate(t0=t0, times=times,
        rinit=rinit,
        rprocess=discrete_time(rproc,delta.t=1),
        rmeasure=rmeas,
        dmeasure=dmeas,
        params=c(kappa=kappa, gamma=gamma, tau=tau),
        statenames = "s",
        obsnames = "r",
        paramnames = c("kappa", "gamma", "tau")
    )
}

