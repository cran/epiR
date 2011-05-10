"epi.conf" <- function(dat, ctype = "mean.single", method, N, design = 1, conf.level = 0.95){
   
    # Define function to calculate confidence interval for a single proportion.
    # This is used on several occasions in this function:
    
    .propsingle <- function(tmp.dat, conf.level = conf.level){
         if (is.matrix(tmp.dat) == FALSE) 
           stop("Error: dat must be a two-column matrix")
        # The method implemented here follows Altman et al (2000) p 46:
        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)
        
        r <- tmp.dat[,1]
        n <- tmp.dat[,1] + tmp.dat[,2] 
        
        p <- r/n
        q <- 1 - r/n
        se <- sqrt((p * q) / n)
        A <- (2 * r) + (z * z)
        B <- z * sqrt((z * z) + (4 * r * q))
        C <- 2 * (n + (z * z))
        low <- (A - B) / C
        up <- (A + B) / C
        tmp.rval <- as.data.frame(cbind(p, se, low, up))
        names(tmp.rval) <- c("est", "se", "lower", "upper")
        tmp.rval
        }
        
     if(ctype == "mean.single"){
        if (is.vector(dat) == FALSE) 
           stop("Error: dat must be a vector")
        mean <- mean(dat)
        n <- length(dat)
        var <- var(dat)
        sd <- sqrt(var)
        se <- sd/sqrt(n)
        
        P <- (1 - conf.level)/2
        t <- abs(qt(P, n - 1))
          
        low <- mean - (t * se)
        up <- mean + (t * se)
        rval <- as.data.frame(cbind(mean, se, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
     }
     
     if(ctype == "mean.unpaired"){
     if (is.data.frame(dat) == FALSE) 
           stop("Error: dat must be a two-column data frame")
        n <- as.vector(by(dat[,2], dat[,1], length))
     
     if (length(n) > 2) 
           stop("Error: there must be only two groups")
     if (is.factor(dat[,1] == FALSE)) 
           stop("Error: the first column of the data frame must be factor")
        sum <- as.vector(by(dat[,2], dat[,1], sum))
        mean <- as.vector(by(dat[,2], dat[,1], mean))
        mean.diff <- mean[1] - mean[2]   
        var <- as.vector(by(dat[,2], dat[,1], var))
        s <- sqrt((((n[1] - 1) * var[1]) + ((n[2] - 1) * var[2])) / (n[1] + n[2] - 2))
        se.diff <- s * sqrt(1/n[1] + 1/n[2])
        
        P <- (1 - conf.level)/2
        t <- abs(qt(P, (n[1] + n[2] - 2)))
        
        low <- mean[1] - mean[2] - (t * se.diff)
        up <- mean[1] - mean[2] + (t * se.diff)
        rval <- as.data.frame(cbind(mean[1] - mean[2], se.diff, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
     }
     
     if(ctype == "mean.paired"){
     if (is.data.frame(dat) == FALSE) 
           stop("Error: dat must be a two-column data frame")
        diff <- as.vector(dat[,2] - dat[,1])
        n <- length(dat[,1])
        mean.diff <- mean(diff)
        sd.diff <- sqrt(var(diff))
        se.diff <- mean.diff/sqrt(n)
        
        P <- (1 - conf.level)/2
        t <- abs(qt(P, (n - 1)))
        
        low <- mean.diff - (t * se.diff)
        up <- mean.diff + (t * se.diff)
        rval <- as.data.frame(cbind(mean.diff, se.diff, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
     }
          
     if(ctype == "prop.single"){
        rval <- .propsingle(tmp.dat = dat, conf.level = conf.level)
        rval
        }
     
     if(ctype == "prop.unpaired"){
     if (is.matrix(dat) == FALSE) 
           stop("Error: dat must be a four-column matrix")
        
        # Work out the confidence interval for each proportion:
        prop.1 <- .propsingle(tmp.dat = matrix(dat[,1:2], ncol = 2), conf.level = conf.level)
        n1 <- dat[,1] + dat[,2]
        p1 <- prop.1[,1]
        l1 <- prop.1[,3]
        u1 <- prop.1[,4]
        
        prop.2 <- .propsingle(tmp.dat = matrix(dat[,3:4], ncol = 2), conf.level = conf.level)
        n2 <- dat[,3] + dat[,4]
        p2 <- prop.2[,1]
        l2 <- prop.2[,3]
        u2 <- prop.2[,4]
        
        # Altman's recommended method (p 48 - 49):
        D <- p1 - p2
        se.D <- sqrt(((p1*(1 - p1))/n1) + ((p2*(1 - p2))/n2))
        low <- D - sqrt((p1 - l1)^2 + (u2 - p2)^2)
        up <- D + sqrt((p2 - l2)^2 + (u1 - p1)^2)
        
        rval <- as.data.frame(cbind(D, se.D, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
     }
     
     if(ctype == "prop.paired"){
        if (is.matrix(dat) == FALSE) 
           stop("Error: dat must be a four-column matrix")
       
        n <- dat[,1] + dat[,2] + dat[,3] + dat[,4]
        r <- dat[,1]
        s <- dat[,2]
        t <- dat[,3]
        u <- dat[,4]
        
        p1 <- (r + s) / n
        p2 <- (r + t) / n
        D <- (s - t) / n
        A <- (r + s) * (t + u) * (r + t) * (s + u)
        B <- (r * u) - (s * t)
        
        se.D <- 1/n * sqrt(s + t - ((s - t)^2 / n))
 
        # Select an appropriate value for C:
        if(B > n/2) C <- B - n/2
        if(B >= 0 & B <= n/2) C <- 0
        if(B < 0) C <- B

        # Calculate phi:
        phi <- C / sqrt(A)
        
        # Set phi to zero if one of the following conditions are true:
        if(r + s == 0) phi <- 0
        if(t + u == 0) phi <- 0
        if(r + t == 0) phi <- 0
        if(s + u == 0) phi <- 0        

        # Calculate confidence intervals for the raw proportions:
        tmp.dat <- matrix(c((r + s), (n - (r + s))), ncol = 2)
        prop.1 <- .propsingle(tmp.dat, conf.level = conf.level)
        l1 <- prop.1[,3]
        u1 <- prop.1[,4]
        
        tmp.dat <- matrix(c((r + t), (n - (r + t))), ncol = 2)
        prop.2 <- .propsingle(tmp.dat, conf.level = conf.level)
        l2 <- prop.2[,3]
        u2 <- prop.2[,4]
        
        # Altman's recommended method (p 52):
        low <- D - sqrt((p1 - l1)^2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2)^2)
        up <- D + sqrt((p2 - l2)^2 - 2 * phi * (p2 - l2) * (u1 - p1) + (u1 - p1)^2)  
                
        rval <- as.data.frame(cbind(D, se.D, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
        }
     
     if(ctype == "inc.risk" | ctype == "prevalence"){
        if (is.matrix(dat) == FALSE) 
           stop("Error: dat must be a two-column matrix")
          
          if(method == "exact"){
          # Exact method (see http://www.folkesundhed.au.dk/uddannelse/software):
          N. <- 1 - ((1 - conf.level) / 2)
          a <- dat[,1]
          n <- dat[,2]
          b <- n - a
          p <- a / n
          
          # Exact binomial confidence limits (D. Collett (1999): Modelling binary data. Chapman & Hall/CRC, Boca Raton Florida, p. 24).
          a. <- ifelse(a == 0, a + 1, a); b. <- ifelse(b == 0, b + 1, b) 
          low <- a. /(a. + (b. + 1) * (1 / qf(1 - N., 2 * a., 2 * b. + 2)))
          up <- (a. + 1) / (a. + 1 + b. / (1 / qf(1 - N., 2 * b., 2 * a. + 2)))
          low <- ifelse(a == 0, 0, low)
          up <- ifelse(a == n, 1, up)

          rval <- as.data.frame(cbind(p, low, up))
          names(rval) <- c("est", "lower", "upper")
          rval
          }
        
          if(method == "wilson"){
          # Wilson's method (see Rothman, Epidemiology An Introduction, page 132): 
          N. <- 1 - ((1 - conf.level) / 2)
          z <- qnorm(N., mean = 0, sd = 1)
          a <- dat[,1]
          n <- dat[,2]
          p <- dat[,1] / dat[,2]
        
          a. <- n/(n + z^2)
          b. <- a/n
          c. <- z^2/(2 * n)
          d. <- (a * (n - a)) / n^3
          e. <- z^2 / (4 * n^2)
          var.wil <- sqrt(d. + e.)
          
          # Design effect equals [var.obs] / [var.srs]. 
          # var.wil has been computed assuming simple random sampling so if an argument for design effect is provided we need to adjust se.wil accordingly:
          se.wil <- sqrt(design * var.wil)
          low <- a. * (b. + c. - (z * se.wil))
          up <- a. * (b. + c. + (z * se.wil))

          rval <- as.data.frame(cbind(p, se.wil, low, up))
          names(rval) <- c("est", "se", "lower", "upper")
          rval
          }
          
          if(method == "fleiss"){
          # Sampling for Epidemiologists, Kevin M Sullivan
          a <- dat[,1]
          n <- dat[,2]
          p <- dat[,1] / dat[,2]
          q <- (1 - p)
          var.fl <- ((p * q) / (n - 1)) * ((N - n) / N)
          # Design effect equals [var.obs] / [var.srs]. 
          # var.fl has been computed assuming simple random sampling so if an argument for design effect is provided we need to adjust se.wil accordingly:
          se.fl <- sqrt(design * var.fl)
          
          df <- n - 1
          N. <- 1 - ((1 - conf.level) / 2)
          t <- abs(qt(p = N., df = df))
          low <- p - (t * se.fl)
          up <- p + (t * se.fl)
          
          rval <- as.data.frame(cbind(p, se.fl, low, up))
          names(rval) <- c("est", "se", "lower", "upper")
          rval
          }
     }
     
     if(ctype == "inc.rate"){
     if (is.matrix(dat) == FALSE) 
           stop("Error: dat must be a two-column matrix")
          
          if(method == "exact"){
          # Exact method (see http://www.folkesundhed.au.dk/uddannelse/software):          
          N. <- 1 - ((1 - conf.level) / 2)
          a <- dat[,1]
          n <- dat[,2]
          p <- a / n
          # If numerator equals zero set lower bound of confidence limit to zero:
          low <- ifelse(a == 0, 0, (0.5 * qchisq(p = N., df = 2 * a + 2, lower.tail = FALSE) / n))
          up <- 0.5 * qchisq(p = 1 - N., df = 2 * a + 2, lower.tail = FALSE) / n
          
          rval <- as.data.frame(cbind(p, low, up))
          names(rval) <- c("est", "lower", "upper")
          rval
          }          
          
          if(method == "byar"){
          # Byar's method (see Rothman, Epidemiology An Introduction, page 134): 
          N. <- 1 - ((1 - conf.level) / 2)
          z <- qnorm(N., mean = 0, sd = 1)
          a.prime <- dat[,1] + 0.5
          p <- dat[,1]/dat[,2]
          PT <- dat[,2]
          low <- (a.prime * (1 - (1/(9 * a.prime)) - (z/3 * sqrt(1/a.prime)))^3)/PT
          up <- (a.prime * (1 - (1/(9 * a.prime)) + (z/3 * sqrt(1/a.prime)))^3)/PT
        
          rval <- as.data.frame(cbind(p, low, up))
          names(rval) <- c("est", "lower", "upper")
          rval
          }
     }
     
     else if(ctype == "smr"){
     if (is.matrix(dat) == FALSE) 
           stop("Error: dat must be a two-column matrix")
        
        # After Dobson et al. 1991. Adapted from Excel code written by Iain Buchan
        # Public Health Informatics at the University of Manchester (www.phi.man.ac.uk)
        # buchan@man.ac.uk
        # dat[,1] = obs; dat[,2] = pop
        
        N. <- 1 - ((1 - conf.level) / 2)
        obs <- dat[,1]
        exp <- (sum(dat[,1]) / sum(dat[,2])) * dat[,2]
        
        smr <- obs / exp
        se.smr <- sqrt(dat[,2]) / exp
        
        low <- ifelse(dat[,1] > 0, 
           ((qchisq(N., df = 2 * dat[,1], lower.tail = FALSE) / 2) / exp), 0)
        up <- ifelse(dat[,1] > 0, 
           ((qchisq(1 - N., df = 2 * (dat[,1] + 1), lower.tail = FALSE) / 2) / exp), 
           ((qchisq(1 - N., df = 2, lower.tail = FALSE) / 2) / exp))
        
        rval <- as.data.frame(cbind(smr, se.smr, low, up))
        names(rval) <- c("est", "se", "lower", "upper")
        rval
     }
     
   return(rval)
}
    
