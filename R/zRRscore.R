zRRscore <- function(dat, conf.level){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  # See https://users.stat.ufl.edu/~aa/cda/R/two-sample/R2/index.html and function PropCIs::riskscoreci
  a <- as.numeric(dat[1])
  b <- as.numeric(dat[3])
  c <- as.numeric(dat[2])
  d <- as.numeric(dat[4])
  
  # a <- as.numeric(dat[1]) 
  # b <- as.numeric(dat[2]) 
  # c <- as.numeric(dat[3]) 
  # d <- as.numeric(dat[4])
  
  N1 <- a + b
  N0 <- c + d
  
  scRR.p <- (a / N1) / (c / N0)
  
  if ((c == 0) && (a == 0)){
    scRR.ll = 0
    scRR.ul = Inf
  }
  
  else{  
    a1 <-  N0 * (N0 * (N0 + N1) * a + N1 * (N0 + a) * (z^2))
    
    a2 <- -N0 * (N0 * N1 * (c + a) + 2 * (N0 + N1) * c * a + N1 * (N0 + c + 2 * a) * (z^2))     
    a3 <- 2 * N0 * N1 * c * (c + a) + (N0 + N1) * (c^2) * a + N0 * N1 * (c + a) * (z^2)
    
    a4 <- -N1 * (c ^ 2) * (c + a)
    
    b1 <- a2 / a1
    b2 <- a3 / a1
    b3 <- a4 / a1
    c1 <- b2 - (b1^2) / 3
    c2 <- b3 - b1 * b2 / 3 + 2 * (b1^3) / 27
    
    # Argument x for acos must be between -1 and +1. Changed 280722.
    acos.x <- sqrt(27) * c2 / (2 * c1 * sqrt(-c1))
    acos.x <- ifelse(acos.x <= -1 | acos.x >= 1, 0, acos.x)
    ceta <- acos(acos.x)
    
    t1 <- -2 * sqrt(-c1 / 3) * cos(pi / 3 - ceta / 3)
    t2 <- -2 * sqrt(-c1 / 3) * cos(pi / 3 + ceta / 3)
    t3 <- 2 * sqrt(-c1 / 3) * cos(ceta / 3)
    p01 <- t1 - b1 / 3
    p02 <- t2 - b1 / 3
    p03 <- t3 - b1 / 3
    p0sum <- p01 + p02 + p03
    p0up <- min(p01, p02, p03)
    p0low <- p0sum - p0up - max(p01, p02, p03)
    
    if((c == 0) && (a != 0)){
      scRR.ll <- (1 - (N1 - a) * (1 - p0low) / (c + N1 - (N0 + N1) * p0low)) / p0low 
      scRR.ul <- Inf 
    }
    
    else if((c != N0) && (a == 0)){
      scRR.ll <- 0
      scRR.ul <- (1 - (N1 - a) * (1 - p0up) / (c + N1 - (N0 + N1) * p0up)) / p0up
    }
    
    else if((c == N0) && (a == N1)){
      scRR.ll <- N1 / (N1 + z^2)
      scRR.ul <- (N0 + z^2) / N0
    }
    
    else if((a == N1) || (c == N0)){
      if((c == N0) && (a == 0)) {
        scRR.ll <- 0
        }
      
      if((c == N0) && (a != 0)) {
        x1 <- a
        n1 <- N1
        x2 <- c
        n2 <- N0
        
        phat1  <- x2 / n2
        phat2  <-  x1 / n1
        phihat <- phat2 / phat1
        phil <- conf.level * phihat
        chi2 <- 0
        while (chi2 <= z){
          a. <- (n2 + n1) * phil
          b. <- -((x2 + n1) * phil + x1 + n2)
          c. <- x2 + x1
          p1hat <- (-b. - sqrt(b.^2 - 4 * a. * c.)) / (2 * a.)
          p2hat <- p1hat * phil
          q2hat <- 1 - p2hat
          var <- (n2 * n1 * p2hat) / (n1 * (phil - p2hat) + n2 * q2hat)
          chi2 <- ((x1 - n1 * p2hat) / q2hat) / sqrt(var)
          scRR.ll <- phil
          phil <- scRR.ll / 1.0001}
        } 
      
      i <- c
      j <- a
      ni <- N0 
      nj <- N1 
      
      if(a == N1){               
        i <- a
        j <- c
        ni <- N1 
        nj <- N0
      } 
      phat1  <- i / ni
      phat2  <-  j / nj
      phihat <- phat2 / phat1
      phiu <- 1.1 * phihat
      
      if((c == N0) && (a == 0)) { 
        if(N0 < 100) {
          phiu = 0.01
          }
        else {phiu = 0.001
        }
      } 
      
      chi1 <- 0
      while (chi1 >= -z){
        a. <- (ni + nj) * phiu
        b. <- -((i + nj) * phiu + j + ni)
        c. <- i + j
        p1hat <- (-b. - sqrt(b.^2 - 4 * a. * c.)) / (2 * a.)
        p2hat <- p1hat * phiu
        q2hat <- 1 - p2hat
        var <- (ni * nj * p2hat) / (nj * (phiu - p2hat) + ni * q2hat)
        chi1 <- ((j - nj * p2hat) / q2hat) / sqrt(var)
        phiu1 <- phiu
        phiu <- 1.0001 * phiu1
      }
      
      if(a == N1) {
        scRR.ll <- 1 / phiu1 
        scRR.ul <- (1 - (N1 - a) * (1 - p0up) / (c + N1 - (N0 + N1) * p0up)) / p0up  
      }
      
      else{scRR.ul <- phiu1}                        
    }   
    
    else{
      scRR.ll <- (1 - (N1 - a) * (1 - p0low) / (c + N1 - (N0 + N1) * p0low)) / p0low 
      scRR.ul <- (1 - (N1 - a) * (1 - p0up) / (c + N1 - (N0 + N1) * p0up)) /p0up
    }
  }  
  
  rval <- c(scRR.p, scRR.ll, scRR.ul)
  return(rval)
}