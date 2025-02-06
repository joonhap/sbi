grid <- seq(-5, 5, by=1)
set.seed(909852)
llfun <- function(pt) {
    -(pt-.2)^2 + .1*pt^3 + 3*rnorm(length(pt))
}
ll <- llfun(grid)
s <- simll(rbind(ll), params=grid)

optDesign.simll(s, init=0, penalty=1)

for (m in 1:100) {
    newpt <- optDesign.simll(s, init=0)
    grid <- c(grid, newpt)
    ll <- c(ll, llfun(newpt))
    s <- simll(rbind(ll), params=grid)
}
    

plot(grid,ll)

plot(grid)
