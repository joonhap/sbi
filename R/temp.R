grid <- seq(-5, 5, by=.5)
set.seed(909852)
llfun <- function(pt) {
    -(pt-.2)^2 + min(4, .2*pt^3) - .1*(pt+1)^4
}
ll <- llfun(grid) + 5*rnorm(length(grid))
s <- simll(rbind(ll), params=grid)

newpt <- 0
optDesign.simll(s, init=newpt, penalty=2)

for (m in 1:50) {
    newpt <- optDesign.simll(s, init=newpt)
    grid <- c(grid, newpt)
    ll <- c(ll, llfun(newpt) + 5*rnorm(1))
    s <- simll(rbind(ll), params=grid)
}


plot(grid,ll)

plot(grid)
