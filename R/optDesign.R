#' @export
optDesign <- function(simll, ...) {
    UseMethod("optDesign")
}

#' Find the next optimal design point for simulation-based inference
#'
#' `optDesign` finds the next design point at which simulation should be carried out for approximately best efficiency in a metamodel-based inference. See Park (2025) for more details on this method. It takes a class `simll` object.
#'
#' @name optDesign
#' @param simll A class `simll` object, containing simulation log likelihoods, the parameter values at which simulations are made, and the weights for those simulations for regression (optional). See help(simll).
#' @param n An integer indicating the number of next Design points to be found
#' @param penalty A positive real that determines how heavily deviation from local quadratic approximation should be penalized. Specifically, a design point is assigned a penalization weight of `exp(-penalty*abs(cubic_term)/quadratic_term)`, where the `quadratic_term` and `cubic_term` respectively indicate the estimated second-order and third-order Taylor expansion terms with respect to the estimated maximum of the expected simulated log-likelihood function. The larger this `penalty` coefficient, the closer the next optimal design points are found near the currently estimated maximizer of the expected simulated log-likelihood. The default value is 2.
#' @param ... Other optional arguments, not currently used.
#'
#' @details
#' This is a generic function, taking a class `simll` object as the first argument.
#' Parameter inference for implicitly defined simulation models can be carried out under a metamodel for the distribution of the log-likelihood estimator.
#' See function `ht` for hypothesis testing and `ci` for confidence interval construction for a one-dimensional parameter.
#' This function `optDesign` founds the next points at which simulations are to be carried out such that the variance of the parameter estimate is reduced approximately the most.
#' This function computes penalization weights given by `exp(-penalty*abs(cubic_term)/quadratic_term)` as described in the explanation for the `penalty` parameter to find design points that are within the scope where local quadratic approximation has relative low error.
#' These penalization weights are multiplied to the original `weights` given to the simulation points specified in the `simll` object.
#' Quadratic regression through the simulated log-likelihoods with these product weights is considered for the selection of next design points.
#'
#' @return A matrix of parameter values at which next simulations are to be carried out for approximately best efficiency.
#' Each row gives a parameter vector.
#'
#' @references Park, J. (2025). Scalable simulation-based inference for implicitly defined models using a metamodel for log-likelihood estimator <https://doi.org/10.48550/arxiv.2311.09446>
#' @export
optDesign.simll <- function(simll, n, penalty=2, ...) {
    validate_simll(simll)

}
