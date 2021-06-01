epi.blcm.paras <- function(ntest.dep = 2, ntest.indep = 1, npop = 2){
  
  ntest.total <- ntest.indep + ntest.dep
  df <- (2^ntest.total - 1) * npop
  npar <- (2 * ntest.total + npop) + (ntest.dep^2 - ntest.dep)
  
  if(npar < df){
    ninf.prior <- 0
    }
  
  if(npar >= df){
    ninf.prior <- npar - df
    }
  
  return(list(ntest.dep = ntest.dep, ntest.indep = ntest.indep, npop = npop,
     df = df, npar = npar, ninf.prior = ninf.prior))
}