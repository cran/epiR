zORscore <- function(dat, conf.level){
  x1 <- dat[1]; n1 <- dat[1] + dat[3]
  x2 <- dat[2]; n2 <- dat[2] + dat[4]
  
  px = x1 / n1
  py = x2 / n2
  
  scOR.p <- (dat[1] / dat[3]) / (dat[2] / dat[4])
  
  if(((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2 == n2))){
    ul = 1/0
    ll = 0   
  } 
  
  else if((x1 == 0) || (x2 == n2)){
    ll = 0
    theta = 0.01 / n2 
    ul = zlimit(x1, n1, x2, n2, conf.level, theta, 1)      
  }
  
  else if((x1 == n1) || (x2 == 0)){
    ul = 1 / 0
    theta = 100 * n1
    ll = zlimit(x1, n1, x2, n2, conf.level, theta, 0)       
  }
  
  else{
    theta = px / (1 - px) / (py / (1 - py)) / 1.1
    ll = zlimit(x1, n1, x2, n2, conf.level, theta, 0)       
    theta = px / (1 - px) / (py / (1 - py)) * 1.1
    ul = zlimit(x1, n1, x2, n2, conf.level, theta, 1)      
  }
  rval <- c(scOR.p, ll,ul)  
  return(rval)
}