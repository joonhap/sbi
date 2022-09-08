rsclR <- function(n, M) {
    z <- rnorm(n)
    v <- rchisq(n, M-1)
    .5*(-z*z - v + M*log(v))
}

set.seed(123); rsclR(10,10)
set.seed(123); rscl(10,10)

microbenchmark(
    rsclR(100000,10),
    rscl(100000,10),
    rscl2(100000,10)
)

