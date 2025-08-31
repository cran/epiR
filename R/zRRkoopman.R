zRRkoopman <- function(dat, conf.level){
  
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  x2 <- qchisq(p = conf.level, df = 1)
  tol <- .Machine$double.eps^0.25
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  
  # a <- as.numeric(dat[1])
  # b <- as.numeric(dat[2])
  # c <- as.numeric(dat[3])
  # d <- as.numeric(dat[4])
  
  x <- a
  m <- a + b
  y <- c
  n <- c + d
  
  if(x == 0 & y == 0) {lower <- 0;  upper <- Inf; est = 0; varhat = NA} 
  else {
    
    a1 <- n * (n * (n + m) * x + m * (n + x) * (z^2))
    a2 <- -n * (n * m * (y + x) + 2 * (n + m) * y * x + m * (n + y + 2 * x) * (z^2))
    a3 <- 2 * n * m * y * (y + x) + (n + m) * (y^2) * x + n * m * (y + x) * (z^2)
    a4 <- -m * (y^2) * (y + x)
    
    b1 <- a2 / a1
    b2 <- a3 / a1
    b3 <- a4 / a1
    
    c1 <- b2 - (b1^2) / 3
    c2 <- b3 - b1 * b2 / 3 + 2 * (b1^3) / 27
    
    ceta <- suppressWarnings(acos(sqrt(27) * c2 / (2 * c1 * sqrt(-c1))))
    t1 <- suppressWarnings(-2 * sqrt(-c1 / 3) * cos(pi / 3 - ceta / 3))
    t2 <- suppressWarnings(-2 * sqrt(-c1 / 3) * cos(pi / 3 + ceta / 3))
    t3 <- suppressWarnings(2 * sqrt(-c1 / 3) * cos(ceta / 3))
    
    p01 <- t1 - b1 / 3
    p02 <- t2 - b1 / 3
    p03 <- t3 - b1 / 3
    
    p0sum <- p01 + p02 + p03
    p0up <- min(p01, p02, p03)
    p0low <- p0sum - p0up - max(p01, p02, p03)
    
    U <- function(a){
      p.hatf <- function(a){
        (a * (m + y) + x + n - ((a * (m + y) + x + n)^2 - 4 * a * (m + n) * (x + y))^0.5) / (2 * (m + n))
      }
      p.hat <- p.hatf(a)
      (((x - m * p.hat)^2) / (m * p.hat * (1 - p.hat))) * (1 + (m * (a - p.hat)) / (n * (1 - p.hat))) - x2
    }
    
    est <- (x / m) / (y / n)
    nrat <- (x / m) / (y / n)
    varhat <- (1 / x) - (1 / m) + (1 / y) - (1 / n) 
    
    if((x == 0) & (y != 0)) {
      nrat <- ((x + 0.5) / m) / (y / n)
      varhat <- (1 / (x + 0.5)) - (1/m) + (1 / y) - (1 / n)
      } 
    
    if((y == 0) & (x != 0)) {
      nrat <- (x / m) / ((y + 0.5) / n)
      varhat <- (1 / x) - (1 / m) + (1 / (y + 0.5)) - (1 / n)
      }
    
    if((y == n) & (x == m)) {
      nrat <- 1
      varhat <- (1 / (m - 0.5)) - (1 / m) + 1 / (n - 0.5) - (1 / n)
      }
    
    La <- nrat * exp(-1 * z * sqrt(varhat)) * 1 / 4
    Lb <- nrat
    Ha <- nrat
    Hb <- nrat * exp(z * sqrt(varhat)) * 4
    
    if((x != 0) & (y == 0)) {
      if(x == m){
        lower <- (1 - (m - x) * (1 - p0low) / (y + m - (n + m) * p0low)) / p0low
        upper <- Inf
      }
      else{         
        lower <- uniroot(U, c(La, Lb), tol=tol)$root
        upper <- Inf
      }
    }
    
    if((x == 0) & (y != n)) {
      upper <- uniroot(U, c(Ha, Hb), tol = tol)$root
      lower <- 0
    }
    
    if(((x == m)|(y == n)) & (y != 0)){
      
      
      if((x == m)&(y == n)){
        U.0 <- function(a){
          if(a <= 1) {m * (1 - a) / a - x2}
          else{(n * (a - 1)) - x2}
        }                                         
        lower <- uniroot(U.0, c(La, est), tol = tol)$root
        upper <- uniroot(U.0, c(est, Hb), tol = tol)$root
      }
      
      if((x == m) & (y != n)){ 
        
        phat1 <- x / m
        phat2 <- y / n
        phihat <- phat2 / phat1
        phiu <- 1.1 * phihat
        r <- 0
        while (r >= -z) {
          a <- (m + n) * phiu
          b <- -((x + n) * phiu + y + m)
          c <- x + y
          p1hat <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
          p2hat <- p1hat * phiu
          q2hat <- 1 - p2hat
          var <- (m * n * p2hat)/(n * (phiu - p2hat) + m * q2hat)
          r <- ((y - n * p2hat) / q2hat) / sqrt(var)
          phiu1 <- phiu
          phiu <- 1.0001 * phiu1
        }
        upper <- (1 - (m - x) * (1 - p0up) / (y + m - (n + m) * p0up)) / p0up
        lower <- 1 / phiu1
      }
      
      if((y == n) & (x != m)){
        p.hat2 <- y / n
        p.hat1 <- x / m
        phihat <- p.hat1 / p.hat2
        phil <- 0.95 * phihat
        r <- 0
        if(x != 0){  
          while(r <= z) {
            a <- (n + m) * phil
            b <- -((y + m) * phil + x + n)
            c <- y + x
            p1hat <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
            p2hat <- p1hat * phil
            q2hat <- 1 - p2hat
            var <- (n * m * p2hat) / (m * (phil - p2hat) + n * q2hat)
            r <- ((x - m * p2hat) / q2hat) / sqrt(var)
            lower <- phil
            phil <- lower / 1.0001
          }
        }
        
        phiu <- 1.1 * phihat
        
        if(x == 0){
          lower <- 0
          phiu <- ifelse(n < 100, 0.01, 0.001)}
        
        r = 0
        while(r >= -z) {
          a <- (n + m) * phiu
          b <- -((y + m) * phiu + x  + n)
          c <- y + x
          p1hat <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
          p2hat <- p1hat * phiu
          q2hat <- 1 - p2hat
          var <- (n * m * p2hat) / (m * (phiu - p2hat) + n * q2hat)
          r <- ((x  - m * p2hat) / q2hat) / sqrt(var)
          phiu1 <- phiu
          phiu <- 1.0001 * phiu1
        }
        upper <- phiu1
      }
    } 
    
    else if((y != n) & (x != m) & (x != 0) & (y != 0)){       
      lower <- uniroot(U, c(La, Lb), tol = tol)$root
      upper <- uniroot(U, c(Ha, Hb), tol = tol)$root
    }
  }
  rval <- c(est, lower, upper)
  return(rval)
}