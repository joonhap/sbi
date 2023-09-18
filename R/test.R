vec3 <- function(vec) {
    d <- length(vec)
    l <- 0
    out <- numeric((d^3+3*d^2+2*d)/6)
    for (k1 in 1:d) {
        for (k2 in 1:k1) {
            out[(l+1):(l+k2)] <- vec[k1]*vec[k2]*vec[1:k2]
            l <- l+k2
        }
    }
    out
}
