epi.ssdxtest <- function(pi, se, sp, epsilon.api, epsilon.ase, epsilon.asp, epsilon.asesp, r = 1, nfractional = FALSE, verbose = FALSE, conf.level = 0.95){
  
  se1 <- se[1]; se2 <- se[2]
  sp1 <- sp[1]; sp2 <- sp[2]
  pi1 <- pi[1]; pi2 <- pi[2]

  alpha <- (1 - conf.level) / 2
  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
  
  # pbars:
  pbar111 <- ((pi1 * se1 * se2) + ((1 - pi1) * (1 - sp1) * (1 - sp2)))
  pbar211 <- ((pi2 * se1 * se2) + ((1 - pi2) * (1 - sp1) * (1 - sp2)))
  pbar112 <- ((pi1 * se1 * (1 - se2)) + ((1 - pi1) * (1 - sp1) * sp2))
  pbar121 <- ((pi1 * se2 * (1 - se1)) + ((1 - pi1) * (1 - sp2) * sp1))
  pbar212 <- ((pi2 * se1 * (1 - se2)) + ((1 - pi2) * (1 - sp1) * sp2))
  pbar221 <- ((pi2 * se2 * (1 - se1)) + ((1 - pi2) * (1 - sp2) * sp1))
  pbar122 <- ((pi1 * (1 - se1) * (1 - se2)) + ((1 - pi1) * sp1 * sp2))
  pbar222 <- ((pi2 * (1 - se1) * (1 - se2)) + ((1 - pi2) * sp1 * sp2))
  
  # Matrix A:
  matrixA <- matrix(rep(0, times = 36), nrow = 6)
  dimnames(matrixA) <- list(c("pi1","pi2","b1","b2","a1","a2"), c("pi1","pi2","b1","b2","a1","a2"))
  
  matrixA[1,1] <- ((((((1 - sp1) * (1 - sp2)) - (se1 * se2))^2) / pbar111) + (((((1 - sp1) * sp2) - (se1 * (1 - se2)))^2) / pbar112) + ((((sp1 * (1 - sp2)) - ((1 - se1) * se2))^2) / pbar121) + ((((sp1 * sp2) - ((1 - se1) * (1 - se2)))^2) / pbar122))
  
  matrixA[2,1] <- 0 
  
  matrixA[3,1] <- ((((1 - sp1) * (1 - sp2) * se2) / pbar111) + ((sp2 * (1 - sp1) * (1 - se2)) / pbar112) + (-((1 - sp2) * (sp1) * se2) / pbar121) + (-(sp1 * (sp2) * (1 - se2)) / pbar122))
  
  matrixA[4,1] <- ((((1 - sp1) * (1 - sp2) * se1) / pbar111) + (((1 - sp2) * (sp1) * (1 - se1)) / pbar121) + (-((sp2) * (1 - sp1) * se1) / pbar112) + (-(sp1 * (sp2) * (1 - se1)) / pbar122))
  
  matrixA[5,1] <- ((((1 - sp2) * se1 * se2) / pbar111) + ((sp2 * se1 * (1 - se2)) / pbar112) + (-( (1 - sp2) * (1 - se1) * se2) / pbar121) + (-(sp1 * (1 - se1) * (1 - se2)) / pbar122))
  
  matrixA[6,1] <- ((((1 - sp1) * se1 * se2) / pbar111) + ((sp1 * se2 * (1 - se1)) / pbar121) + (-( (1 - sp1) * (1 - se2) * se1) / pbar112) + (-(sp2 * (1 - se1) * (1 - se2)) / pbar122))
  
  matrixA[2,2] <- r * ((((((1 - sp1) * (1 - sp2)) - (se1 * se2))^2) / pbar211) + (((((1 - sp1) * sp2) - (se1 * (1 - se2)))^2) / pbar212) + ((((sp1 * (1 - sp2)) - ((1 - se1) * se2))^2) / pbar221) + ((((sp1 * sp2) - ((1 - se1) * (1 - se2)))^2) / pbar222))
  
  matrixA[3,2] <- r * ((((1 - sp1) * (1 - sp2) * se2) / pbar211) + ((sp2 * (1 - sp1) * (1 - se2)) / pbar212) + (-((1 - sp2) * (sp1) * se2) / pbar221) + (-(sp1 * (sp2) * (1 - se2)) / pbar222))
  
  matrixA[4,2] <- r * ((((1 - sp1) * (1 - sp2) * se1) / pbar211) + (((1 - sp2) * (sp1) * (1 - se1)) / pbar221) + (-((sp2) * (1 - sp1) * se1) / pbar212) + (-(sp1 * (sp2) * (1 - se1)) / pbar222))
  
  matrixA[5,2] <- r * ((((1 - sp2) * se1 * se2) / pbar211) + ((sp2 * se1 * (1 - se2)) / pbar212) + (-((1 - sp2) * (1 - se1) * se2) / pbar221) + (-(sp1 * (1 - se1) * (1 - se2)) / pbar222))
  
  matrixA[6,2] <- r * ( ( ( (1-sp1)*se1*se2)/pbar211) + ( (sp1*se2*(1-se1))/pbar221) + ( - ( (1-sp1)*(1-se2)*se1)/pbar212) + ( - ( sp2*(1-se1)*(1-se2))/pbar222))
  
  matrixA[3,3] <- ( ( (pi1^2)*( ((se2^2)/pbar111) + ( ( (1-se2)^2)/pbar112) + ( (se2^2)/pbar121) + (((1-se2)^2)/pbar122))) + ( r * ((pi2)^2) * ( ( (se2^2)/pbar211) + ( ( (1-se2)^2)/pbar212) + ( (se2^2)/pbar221) + ( ( (1-se2)^2)/pbar222))))
  
  matrixA[4,3] <- ( ( (pi1^2)* (((se1*se2)/pbar111) + ( - (se1*(1-se2))/pbar112) + (-((1-se1)*se2)/pbar121) + (((1-se1)*(1-se2))/pbar122))) + ( r * (pi2^2)* (((se1*se2)/pbar211) + ( - (se1*(1-se2))/pbar212) + (-((1-se1)*se2)/pbar221) + (((1-se1)*(1-se2))/pbar222))))
  
  matrixA[5,3] <- - ( ( (pi1*(1-pi1)) * ( (((1-sp2)*se2)/pbar111) + ((sp2*(1-se2))/pbar112) + (((1-sp2)*se2)/pbar121) + ((sp2*(1-se2))/pbar122))) + ( r * (pi2*(1-pi2)) * ( (((1-sp2)*se2)/pbar211) + ((sp2*(1-se2))/pbar212) + (((1-sp2)*se2)/pbar221) + ((sp2*(1-se2))/pbar222))) )
  
  matrixA[6,3] <- ( ( (pi1*(1-pi1)) * ( ( - ( (1-sp1)*se2)/pbar111) + (sp1*se2/pbar121) + ( (1-sp1)*(1-se2)/pbar112) + ( - (sp1*(1-se2))/pbar122))) + ( r * (pi2*(1-pi2)) * ( ( - ( (1-sp1)*se2)/pbar211) + (sp1*se2/pbar221) + ( (1-sp1)*(1-se2)/pbar212) + ( - (sp1*(1-se2))/pbar222))) )
  
  matrixA[4,4] <- ( ( (pi1^2)*( ((se1^2)/pbar111) + ( ( (1-se1)^2)/pbar121) + ( (se1^2)/pbar112) + (((1-se1)^2)/pbar122))) + ( r * ((pi2)^2) * ( ( (se1^2)/pbar211) + ( ( (1-se1)^2)/pbar221) + ( (se1^2)/pbar212) + ( ( (1-se1)^2)/pbar222))))
  
  matrixA[5,4] <- ( ( (pi1*(1-pi1)) * ( ( - ( (1-sp2)*se1)/pbar111) + (sp2*se1/pbar112) + ( (1-sp2)*(1-se1)/pbar121) + ( - (sp2*(1-se1))/pbar122))) + (r * (pi2*(1-pi2)) * ( ( - ( (1-sp2)*se1)/pbar211) + (sp2*se1/pbar212) + ( (1-sp2)*(1-se1)/pbar221) + ( - (sp2*(1-se1))/pbar222))) )
  
  matrixA[6,4] <- - ( ( (pi1*(1-pi1)) * ( (((1-sp1)*se1)/pbar111) + ((sp1*(1-se1))/pbar121) + (((1-sp1)*se1)/pbar112) + ((sp1*(1-se1))/pbar122))) + ( r * (pi2*(1-pi2)) * ( (((1-sp1)*se1)/pbar211) + ((sp1*(1-se1))/pbar221) + (((1-sp1)*se1)/pbar212) + ((sp1*(1-se1))/pbar222))) )
  
  matrixA[5,5] <- ((((1-pi1)^2) * ((((1-sp2)^2)/pbar111) + ((sp2^2)/pbar112) + (((1-sp2)^2)/pbar121) + ((sp2^2)/pbar122))) + (r* ((1-pi2)^2) * ((((1-sp2)^2)/pbar211) + ((sp2^2)/pbar212) + (((1-sp2)^2)/pbar221) + ((sp2^2)/pbar222))))
  
  matrixA[6,5] <- ( ( ((1-pi1)^2) * ( ( (1-sp1)*(1-sp2)/pbar111) + ( - ((1-sp1)*sp2/pbar112)) + ( - (sp1*(1-sp2)/pbar121)) + (sp1*sp2/pbar122))) + ( r * ((1-pi2)^2) * ( ( (1-sp1)*(1-sp2)/pbar211) + ( - ((1-sp1)*sp2/pbar212)) + ( - (sp1*(1-sp2)/pbar221)) + (sp1*sp2/pbar222))) ) 
  
  matrixA[6,6] <- ((((1-pi1)^2) * ((((1-sp1)^2)/pbar111) + ((sp1^2)/pbar121) + (((1-sp1)^2)/pbar112) + ((sp1^2)/pbar122))) + (r * ((1-pi2)^2) * ((((1-sp1)^2)/pbar211) + ((sp1^2)/pbar221) + (((1-sp1)^2)/pbar212) + ((sp1^2)/pbar222))))
  
  matrixA[1,2] <- matrixA[2,1]  
  matrixA[1,3] <- matrixA[3,1]  
  matrixA[1,4] <- matrixA[4,1]  
  matrixA[1,5] <- matrixA[5,1]  
  matrixA[1,6] <- matrixA[6,1]  
  
  matrixA[2,3] <- matrixA[3,2]  
  matrixA[2,4] <- matrixA[4,2]  
  matrixA[2,5] <- matrixA[5,2]  
  matrixA[2,6] <- matrixA[6,2]
  
  matrixA[3,4] <- matrixA[4,3]  
  matrixA[3,5] <- matrixA[5,3]  
  matrixA[3,6] <- matrixA[6,3]  
  
  matrixA[4,5] <- matrixA[5,4]  
  matrixA[4,6] <- matrixA[6,4]
  
  matrixA[5,6] <- matrixA[6,5]    
  
  matrixC <- solve(matrixA)
  
  var.se1se2 <- matrixC[3,3] + matrixC[4,4] - 2 * matrixC[4,3]
  var.sp1sp2 <- matrixC[5,5] + matrixC[6,6] - 2 * matrixC[6,5]
  
  # Confidence interval width equals twice the absolute error:
  wpi <- epsilon.api * 2 
  wse <- epsilon.ase * 2
  wsp <- epsilon.asp * 2
  wsesp <- epsilon.asesp * 2
  w <- c(wpi, wse, wsp, wsesp)
  
  matrixC.d <- c(as.numeric(diag(matrixC)), var.se1se2, var.sp1sp2)
  n1 <- (2 * z.alpha * sqrt(matrixC.d) / w)^2
  n2 <- n1 * r
  
  sample.size <- data.frame(p1 = c(n1[1], n2[1]),
     p2 = c(n1[2], n2[2]),
     se1 = c(n1[3], n2[3]),
     se2 = c(n1[4], n2[4]),
     sp1 = c(n1[5], n2[5]),
     sp2 = c(n1[6], n2[6]),
     se1.se2 = c(n1[7], n2[7]),
     sp1.sp2 = c(n1[8], n2[8]))
  row.names(sample.size) <- c("pop1","pop2")
  
  if(nfractional == TRUE & verbose == TRUE){
    rval <- sample.size
  }
  
  if(nfractional == FALSE & verbose == TRUE){
    rval <- ceiling(sample.size)
  }
  
  if(nfractional == TRUE & verbose == FALSE){
    rval <- data.frame(n = c(max(sample.size[1,]), max(sample.size[2,])))
    row.names(rval) <- c("pop1","pop2")
  }
  
  if(nfractional == FALSE & verbose == FALSE){
    rval <- data.frame(n = c(ceiling(max(sample.size[1,])), ceiling(max(sample.size[2,]))))
    row.names(rval) <- c("pop1","pop2")
  }
  
  return(rval)
}
