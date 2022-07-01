epi.herdtest <- function(se, sp, P, N, n, k){
   
   # Probability of testing positive:
   APpos <- P * se + (1 - P) * (1 - sp)
   APneg <- (1 - sp)
 
   # Binomial distribution recommended if n / N < 0.2
   
   # Binomial distribution:
   bHSe <- 1 - pbinom(q = k - 1, size = n, prob = P)
   bHSp <- pbinom(q = k - 1, size = n, prob = 1 - sp)
   dbinom <- data.frame(APpos = APpos, APneg = APneg, HSe = bHSe, HSp = bHSp)
   
   # Hypergeometric distribution:
   hHSe <- 1 - phyper(k - 1, N * APpos, N - N * APpos, n)
   hHSp <- phyper(k - 1, N * APneg, N - N * APneg, n)
   dhyper <- data.frame(APpos = APpos, APneg = APneg, HSe = hHSe, HSp = hHSp)

   # Sampling fraction:
   sfraction <- n / N
   
   rval <- list(sfraction = sfraction, dbinom = dbinom, dhyper = dhyper)
   rval
}

