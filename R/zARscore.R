zARscore <- function(dat, conf.level, units){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d
  
  sARisk.p <- ((a / N1) - (c / N0))       
  px = a / N1
  py = c / N0
  z = qchisq(conf.level, 1)
  proot = px - py
  dp = 1 - proot
  niter = 1
  while(niter <= 50){
    dp = 0.5 * dp
    up2 = proot + dp
    score = zz2stat(px, N1, py, N0, up2)
    if(score < z){proot = up2}
    niter = niter + 1
    if((dp < 0.0000001) || (abs(z - score) < 0.000001)){
      niter = 51
      ul = up2
    }
  } 
  
  proot = px - py
  dp = 1 + proot
  niter = 1
  while(niter <= 50){
    dp = 0.5 * dp
    low2 = proot - dp
    score = zz2stat(px, N1, py, N0, low2)
    if(score < z){proot = low2}
    niter = niter + 1
    if((dp < 0.0000001) || (abs(z - score) < 0.000001)){
      ll = low2
      niter = 51
    }
  }
  c(sARisk.p * units, ll * units, ul * units)
}