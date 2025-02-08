grid <- seq(-5, 5, by=.5)
set.seed(909852)
llfun <- function(pt) {
    -(pt-.2)^2 + min(4, .2*pt^3) - .1*(pt+1)^4
}
ll <- llfun(grid) + 5*rnorm(length(grid))
#w <- rep(1, length(grid))
s <- simll(rbind(ll), params=grid) #, weights=w)

newpt <- 0
optDesign.simll(s, init=newpt, penalty=2)

for (m in 1:100) {
    opt <- optDesign.simll(s, init=runif(1,-5,5))
    newpt <- unname(opt[["par"]])
    grid <- c(grid, newpt)
    #w <- c(w, unname(opt["w"]))
    ll <- c(ll, llfun(newpt) + 5*rnorm(1))
    s <- simll(rbind(ll), params=grid) #, weights=w)
}

htout <- ht(s, null.value=0, test="MESLE", weight_penalization=TRUE)
MESLEhat <- htout$meta_model_MLE_for_MESLE

estCoef <- c(htout$regression_estimates$a, htout$regression_estimates$b, htout$regression_estimates$c)

Ahat_cubic <- htout$first_stage_cubic_approx

plot(grid,ll)

plotgrid <- seq(range(grid)[1], range(grid)[2], length.out=151)

lines(plotgrid, sapply(plotgrid, function(x) { sum(estCoef*c(1,x,x^2)) }))

theta <- cbind(attr(s, "params")) # coerce into a matrix
theta_mean <- apply(theta, 2, mean)
theta_sd <- apply(theta, 2, sd)
trans_n <- function(vec) { (vec-theta_mean)/theta_sd } # normalize by centering and scaling
trans_b <- function(vec) { vec*theta_sd + theta_mean } # transform back to the original scale

ca <- function(x) { sum(c(vec012(trans_n(x)), vec3(trans_n(x))) * Ahat_cubic) } # cubic approx
Third <- function(x) { sum(vec3(trans_n(x)-trans_n(MESLEhat)) * Ahat_cubic[(dim012+1):dim0123]) } # third order term of the cubic approximation in the form of (multidimensional equivalent of) cubic_coeff*(theta-MESLEhat)^3
tca <- function(x) { ca(x) - Third(x) } # truncated cubic approximation where the third order term is dropped
wpen <- function(point) { # penalty weight due to the discrepancy between quadratic and cubic approx
    exp(-penalty*(ca(point) - tca(point))^2/(tca(point) - tca(MESLEhat))^2)
}
upd_w <- sapply(grid, wpen) # update w by multiplying penalizing weights
## The second stage approximation, accounting for the penalizing weights, follows next.

lines(plotgrid, sapply(plotgrid, ca), lty=2)

lines(plotgrid, sapply(plotgrid, tca), col='red')

plot(grid, sapply(grid, function(point) { (ca(point) - tca(point))^2/(tca(point) - tca(MESLEhat))^2 }))

plot(grid, sapply(grid, function(point) { (ca(point) - tca(point))^2 }))

plot(grid, sapply(grid, function(point) { (tca(point) - tca(MESLEhat))^2 }))

plot(grid, log(htout$updated_weights))

plot(grid)

plot(grid[1:(length(grid)-1)], log(opt[["allweights"]]))
