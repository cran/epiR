## ---- echo = FALSE, message = FALSE-------------------------------------------

# If you want to create a PDF document paste the following after line 9 above:
#   pdf_document:
#     toc: true
#     highlight: tango
#     number_sections: no
#     latex_engine: xelatex    
# header-includes: 
#    - \usepackage{fontspec}

knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----ssrs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------
library(pander)
panderOptions('table.split.table', Inf)
# panderOptions('table.alignment.default', function(df) ifelse(sapply(df, is.numeric), 'right', 'left'))

set.caption("Functions to estimate sample size using representative population sampling data.")

ssrs.tab <- " 
Sampling       | Outcome               | Details                    | Function
Representative | Prob disease freedom  | Imperfect Se, perfect Sp   | `rsu.sspfree.rs`
Representative | SSe                   | Imperfect Se, perfect Sp   | `rsu.sssep.rs`
Two stage representative | SSe         | Imperfect Se, perfect Sp   | `rsu.sssep.rs2st`
Representative | SSe                   | Imperfect Se, imperfect Sp, known N | `rsu.sssep.rsfreecalc`
Pooled representative    | SSe                   | Imperfect Se, imperfect Sp | `rsu.sssep.rspool`"

ssrs.df <- read.delim(textConnection(ssrs.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(ssrs.df) <- unname(as.list(ssrs.df[1,])) # put headers on
ssrs.df <- ssrs.df[-1,] # remove first row
row.names(ssrs.df) <- NULL
pander(ssrs.df, style = 'rmarkdown')

## ----message = FALSE----------------------------------------------------------
library(epiR)
rsu.sssep.rs(N = NA, pstar = 0.05, se.p = 0.95, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sspfree.rs(N = NA, prior = 0.50, p.intro = 0.01, pstar = 0.05, pfree = 0.95, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 500, pstar = 0.05, se.p = 0.95, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rsfreecalc(N = 5000, pstar = 0.05, mse.p = 0.95, 
   msp.p = 0.95, se.u = 0.95, sp.u = 0.98, method = "hypergeometric", 
   max.ss = 32000)$summary

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rsfreecalc(N = 5000, pstar = 0.10, mse.p = 0.95, 
   msp.p = 0.95, se.u = 0.95, sp.u = 0.98, method = "hypergeometric", 
   max.ss = 32000)$summary

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 20000, pstar = 0.005, se.p = 0.95, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = NA, pstar = 0.05, se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs2st(H = 20000, N = NA, pstar.c = 0.005, pstar.u = 0.05, se.p = 0.95, se.c = 0.95, se.u = 0.90)

## ----seprs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----

set.caption("Functions to estimate surveillance system sensitivity (SSe) using representative population sampling data.")

seprs.tab <- " 
Sampling       | Outcome      | Details                    | Function
Representative | SSe          | Imperfect Se, perfect Sp   | `rsu.sep.rs`
Two stage representative      | SSe          | Imperfect Se, perfect Sp   | `rsu.sep.rs2st`
Representative | SSe          | Imperfect Se, perfect Sp, multiple components   | `rsu.sep.rsmult`
Representative | SSe          | Imperfect Se, imperfect Sp | `rsu.sep.rsfreecalc`
Pooled representative         | SSe   | Imperfect Se, perfect Sp   | `rsu.sep.rspool`
Representative | SSe          | Imperfect Se, perfect Sp   | `rsu.sep.rsvarse`
Representative | SSp          | Imperfect Sp               | `rsu.spp.rs`"

seprs.df <- read.delim(textConnection(seprs.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(seprs.df) <- unname(as.list(seprs.df[1,])) # put headers on
seprs.df <- seprs.df[-1,] # remove first row
row.names(seprs.df) <- NULL
pander(seprs.df, style = 'rmarkdown')

## ----message = FALSE----------------------------------------------------------
rsu.sep.rs(N = NA, n = 300, pstar = 0.01, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
N <- seq(from = 80, to = 100, by = 5)
n <- rep(30, times = length(N))

herd.sep <- rsu.sep.rs(N = N, n = n, pstar = 0.01, se.u = 0.95)
sort(round(herd.sep, digits = 2))

## ----message = FALSE----------------------------------------------------------
# Diagnostic test sensitivities and the number of samples tested at each laboratory:
se.t1 <- 0.80; se.t2 <- 0.70
n.lab1 <- 50; n.lab2 <- 23

# Create a vector of test sensitivities for each sample:
se.all <- c(rep(se.t1, times = n.lab1), rep(se.t2, times = n.lab2))
rsu.sep.rsvarse(N = n.lab1 + n.lab2, pstar = 0.05, se.u = se.all)

## ----pfreers.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----

set.caption("Functions to estimate the probability of disease freedom using representative population sampling data.")

pfreers.tab <- " 
Sampling       | Outcome                             | Details                    | Function
Representative | Prob disease of freedom | Imperfect Se, perfect Sp   | `rsu.pfree.rs`
Representative | Equilibrium prob of disease freedom | Imperfect Se, perfect Sp   | `rsu.pfree.equ`"

pfreers.df <- read.delim(textConnection(pfreers.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(pfreers.df) <- unname(as.list(pfreers.df[1,])) # put headers on
pfreers.df <- pfreers.df[-1,] # remove first row
row.names(pfreers.df) <- NULL
pander(pfreers.df, style = 'rmarkdown')

## ----message = FALSE----------------------------------------------------------
library(ggplot2); library(lubridate); library(scales)

# Define a vector disease incursion probabilities (January to December):
p.intro <- c(0.01,0.01,0.01,0.02,0.04,0.10,0.10,0.10,0.08,0.06,0.04,0.02)

rval.df <- rsu.pfree.rs(se.p = rep(0.65, times = 12), p.intro = p.intro, prior = 0.50, by.time = TRUE)

# Re-format rval.df ready for for ggplot2:
dat.df <- data.frame(mnum = rep(1:12, times = 2),
   mchar = rep(seq(as.Date("2020/1/1"), by = "month", length.out = 12), times = 2),                 
   class = c(rep("Disease introduction", times = length(p.intro)), 
             rep("Disease freedom", times = length(p.intro))),
   prob = c(rval.df$PIntro, rval.df$PFree))

# Plot the results:
ggplot(data = dat.df, aes(x = mchar, y = prob, group = class, col = class)) +
  theme_bw() +
  geom_point() + 
  geom_line() +
  scale_colour_manual(values = c("red", "dark blue")) + 
  scale_x_date(breaks = date_breaks("1 month"), labels = date_format("%b"),
     name = "Month") +
  scale_y_continuous(limits = c(0,1), name = "Probability") +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", col = "blue") +
  guides(col = guide_legend(title = "")) +
  theme(legend.position = c(0.8, 0.5))

## ----ssrb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------
set.caption("Functions to estimate sample size using risk based sampling data.")

ssrb.tab <- " 
Sampling       | Outcome               | Details                  | Function
Risk-based     | SSe                   | Single Se for risk groups, perfect Sp        | `rsu.sssep.rbsrg`
Risk-based     | SSe                   | Multiple Se within risk groups, perfect Sp   | `rsu.sssep.rbmrg`
Risk-based     | SSe                   | Two stage sampling, 1 risk factor  | `rsu.sssep.rb2st1rf`
Risk-based     | SSe                   | Two stage sampling, 2 risk factors | `rsu.sssep.rb2st2rf`"

ssrb.df <- read.delim(textConnection(ssrb.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(ssrb.df) <- unname(as.list(ssrb.df[1,])) # put headers on
ssrb.df <- ssrb.df[-1,] # remove first row
row.names(ssrb.df) <- NULL
pander(ssrb.df, style = 'rmarkdown')

## ----message = FALSE----------------------------------------------------------
# Matrix listing the proportions of samples for each test in each risk group (the number of rows equal the number of risk groups, the number of columns equal the number of tests):
m <- rbind(1,1,1)

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
# Matrix listing the proportions of samples for each test in each risk group (the number of rows equal the number of risk groups, the number of columns equal the number of tests):
m <- rbind(c(0.8,0.2), c(0.5,0.5), c(0.7,0.3))

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = c(0.92,0.80))

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rbsrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.10,0.10,0.80), 
   spr = c(0.50,0.30,0.20), se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rb2st2rf(
   rr.c = c(5,3,1), ppr.c = c(0.10,0.20,0.70), spr.c = c(0.40,0.40,0.20),
   pstar.c = 0.02,
   rr.u = c(4,1), ppr.u = c(0.1, 0.9), spr.u = c(1,0),
   pstar.u = 0.10, 
   se.p = 0.95, se.c = 0.95, se.u = 0.95)

## ----seprb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----
set.caption("Functions to estimate surveillance system sensitivity using risk based sampling data.")

ssrb.tab <- " 
Sampling       | Outcome        | Details                  | Function
Risk-based     | SSe            | Varying Se, perfect Sp   | `rsu.sep.rb`
Risk-based     | SSe            | Varying Se, perfect Sp, one risk factor   | `rsu.sep.rb1rf`
Risk-based     | SSe            | Varying Se, perfect Sp, two risk factors  | `rsu.sep.rb2rf`"

ssrb.df <- read.delim(textConnection(ssrb.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(ssrb.df) <- unname(as.list(ssrb.df[1,])) # put headers on
ssrb.df <- ssrb.df[-1,] # remove first row
row.names(ssrb.df) <- NULL
pander(ssrb.df, style = 'rmarkdown')

## ----message = FALSE----------------------------------------------------------
# There are 1784 herds in the study area:
H <- 1784

# Twenty of the 1784 herds are sampled. Generate 20 herds of varying size:
set.seed(1234)
hsize <- rlnorm(n = 20, meanlog = log(10), sdlog = log(8))
hsize <- round(hsize + 20, digits = 0)

# Generate a matrix listing the number of growers and finishers in each of the 20 sampled herds. 
# Assume that anywhere between 80% and 95% of the pigs in each herd are growers:
set.seed(1234)
pctg <- runif(n = 20, min = 0.80, max = 0.95)
ngrow <- round(pctg * hsize, digits = 0)
nfini <- hsize - ngrow
N <- cbind(ngrow, nfini)

# Generate a matrix listing the number of grower and finisher pigs sampled from each herd. Fifteen pigs from each herd are sampled. If there's less than 15 pigs we sample the entire herd:
nsgrow <- rep(0, times = 20)
nsfini <- ifelse(nfini <= 15, nfini, 15)
n <- cbind(nsgrow, nsfini)

# The herd-level design prevalence is 0.01 and the individual pig-level design prevalence is 0.05: 
pstar.c <- 0.01
pstar.u <- 0.05

# For herds in the high-risk area the probability being A. hyopneumoniae positive is 3.5 times that of herds in the low-risk area. Ninety percent of herds are in the low risk area and 10% are in the high risk area:
rr.c <- c(3.5,1)
ppr.c <- c(0.1,0.9) 

# We've sampled 15 herds from the high risk area and 5 herds from the low risk area. Above, for vector rr.c, the relative risk for the high risk group is listed first so the vector rg follows this order:
rg <- c(rep(1, times = 15), rep(2, times = 5))

# The probability being A. hyopneumoniae positive for finishers is 5 times that of growers. For the matrices N and n growers are listed first then finishers. Vector rr.u follows the same order:
rr.u <- c(1,5)

# The diagnostic sensitivity of the A. hyopneumoniae ELISA is 0.95:
se.u <- 0.95

rsu.sep.rb2st(H = H, N = N, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = NA,
   se.u = se.u)

## ----message = FALSE----------------------------------------------------------
# Generate a matrix listing the proportion of growers and finishers in each of the 20 sampled herds:

ppr.u <- cbind(rep(0.9, times = 20), rep(0.1, times = 20))

# Set H (the number of clusters) and N (the number of surveillance units within each cluster) to NA:
rsu.sep.rb2st(H = NA, N = NA, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = ppr.u,
   se.u = se.u)

## ----message = FALSE----------------------------------------------------------
rsu.sep.pass(step.p = c(0.10,0.20,0.90,0.99), pstar.c = 0.01,
   p.inf.u = 0.98, N = 1000, n = 5, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.adjrisk(rr = c(5,1), ppr = c(0.10,0.90))

## ----message = FALSE----------------------------------------------------------
rsu.pstar(N = 193000, n = 7764, se.p = 0.95, se.u = 0.85) * 100000

## ----message = FALSE----------------------------------------------------------
rsu.sep(N = 193000, n = 7764, pstar = 10 / 100000, se.u = 0.85)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 193000, pstar = 10 / 100000, se.p = 0.95, se.u = 0.85)

