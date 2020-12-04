zsep.binom <- function(n, pstar, se = 1, sp = 1) {
    sep <- 1 - (1 - (se * pstar + (1 - sp) * (1 - pstar)))^n
    return(sep)
}
