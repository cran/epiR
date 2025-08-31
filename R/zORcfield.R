zORcfield <- function (dat, conf.level, interval = c(1e-08, 1e+08)){
  
#   dFNCHypergeo <- function(x, m1, m2, n, odds, precision = 1e-07){
#    stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2),
#    is.numeric(n), is.numeric(odds), is.numeric(precision))
#    .Call("dFNCHypergeo", as.integer(x), as.integer(m1), as.integer(m2),
#    as.integer(n), as.double(odds), as.double(precision),
#   PACKAGE = "BiasedUrn")
# }
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d
  
  cfOR.p <- (a / b) / (c / d)
  
  if (((a == 0) && (c == 0)) || ((a == N1) && (c == N0))) {
    ll <- 0
    ul <- Inf
  }
  
  else if (c == N0 || a == 0) {
    ll <- 0
    ul <- uniroot(function(or) {
      sum(sapply(max(0, a + c - N0):a, BiasedUrn::dFNCHypergeo, N1, N0, a + c, or)) - BiasedUrn::dFNCHypergeo(a, N1, N0, a + c, or)/2 - (1- conf.level)/2
    }, interval = interval)$root
  }
  else if (a == N1 || c == 0) {
    ll <- uniroot(function(or) {
      sum(sapply(a:min(N1, a + c), BiasedUrn::dFNCHypergeo, N1, N0, a + c, or)) - BiasedUrn::dFNCHypergeo(a, N1, N0, a + c, or)/2 - (1 - conf.level)/2
    }, interval = interval)$root
    ul <- Inf
  }
  else {
    ll <- uniroot(function(or) {
      sum(sapply(a:min(N1, a + c), BiasedUrn::dFNCHypergeo, N1, N0, a + c, or)) - BiasedUrn::dFNCHypergeo(a, N1, N0, a + c, or)/2 - (1 - conf.level)/2
    }, interval = interval)$root
    ul <- uniroot(function(or) {
      sum(sapply(max(0, a + c - N0):a, BiasedUrn::dFNCHypergeo, N1, N0, a + c, or)) - BiasedUrn::dFNCHypergeo(a, N1, N0, a + c, or)/2 - (1 - conf.level)/2
    }, interval = interval)$root
  }
  rval <- c(cfOR.p, ll, ul)
  return(rval)
}