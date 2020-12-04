\name{rsu.sep.rb2rf}

\alias{rsu.sep.rb2rf}

\title{
Surveillance system sensitivity assuming risk-based sampling on two risk factors
}

\description{
Calculates risk-based surveillance system (population-level) sensitivity with a two risk factors, assuming [one-stage] risk-based sampling and allowing unit sensitivity to vary among risk strata.
}

\usage{
rsu.sep.rb2rf(N, n, rr1, ppr1, rr2, ppr2, pstar, se.u, method = "binomial")
}

\arguments{
\item{N}{matrix of population sizes for each risk group. Rows = levels of \code{rr1}, columns = levels of \code{rr2}.}
\item{n}{matrix of the number of surveillance units tested in each risk group. Rows = levels of \code{rr1}, columns = levels of \code{rr2}.}
\item{rr1}{scalar or vector defining the first set of relative risk values.}
\item{ppr1}{scalar or vector of the same length as that vector of \code{rr1} defining the population proportions in each of the first risk strata. Proportions must sum to one. Ignored if \code{method = "hypergeometric"}.}
\item{rr2}{matrix of relative risks for the second risk factor. Rows = levels of \code{rr1}, columns = levels of \code{rr2}.}
\item{ppr2}{scalar or vector of the same length as that vector of \code{rr1} defining the population proportions in each of the second risk strata. Row proportions must sum to one. Ignored if \code{method = "hypergeometric"}.}
\item{pstar}{scalar, defining the design prevalence.}
\item{se.u}{scalar or vector of the same length as that vector of \code{rr1} defining the unit sensitivity (which can vary across strata).}
\item{method}{character string indicating the method to be used. Options are \code{binomial} or \code{hypergeometric}. See details, below.}
}

\details{
If \code{method = binomial} \code{N} is ignored and values for \code{ppr} need to be entered. Conversely, if \code{method = hypergeometric}, \code{ppr} is ignored and calculated from \code{N}.
}

\value{
A list comprised of two elements:

\item{se.p}{scalar, surveillance system (population-level) sensitivity estimates.}
\item{epi}{vector, effective probability of infection estimates.}
\item{adj.risk1}{vector, adjusted relative risk estimates for the first risk factor.}
\item{adj.risk2}{vector, adjusted relative risk estimates for the second risk factor.}
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## disease using risk based sampling. Assume a design prevalence of 0.01 
## at the surveillance unit level. Surveillance units are categorised as 
## being either high or low risk with the probability of disease for 
## high risk surveillance units 3 times the probability of disease for low 
## risk units. The proportion of units in each risk group is 0.20 and 0.80,
## respectively.

## Within each of the two risk categories the probability of disease varies 
## with age with younger age groups having four times the risk of disease 
## as older age groups. In the high risk area 10\% of the population are young 
## and 90\% are old. In the low risk area 30\% of the population are young and 
## 70\% are old. 

## The total number of surveillance units in the population is unknown. The 
## numbers of young and old surveillance units tested in the high and low risk 
## groups are 40, 20, 20 and 10, respectively. You intend to use a test with 
## diagnostic sensitivity of 0.80. What is the surveillance system sensitivity?

rsu.sep.rb2rf(N = NA, n = rbind(c(40,20), c(20,10)), 
   rr1 = c(3,1), 
   ppr1 = c(0.20,0.80), 
   rr2 = rbind(c(4,1), c(4,1)), 
   ppr2 = rbind(c(0.10,0.90), c(0.30,0.70)),
   pstar = 0.01, 
   se.u = 0.80, method = "binomial")$se.p

## The surveillance system sensitivity is 0.93.


## EXAMPLE 2:
## This example shows the importance of sampling high risk groups. Take the 
## same scenario as above but switch the relative proportions sampled by 
## risk group --- taking a greater number of samples from the low risk group
## compared with the high risk group:

rsu.sep.rb2rf(N = NA, n = rbind(c(10,20), c(20,40)), 
   rr1 = c(3,1), 
   ppr1 = c(0.20,0.80), 
   rr2 = rbind(c(4,1), c(4,1)), 
   ppr2 = rbind(c(0.10,0.90), c(0.30,0.70)),
   pstar = 0.01, 
   se.u = 0.80, method = "binomial")$se.p
 
## The surveillance system sensitivity is 0.69. Here we've taken exactly the 
## same number of samples as Example 1, but there's a substantial decrease
## in surveillance system sensitivity because we've concentrated sampling on 
## a low risk group (decreasing our ability to detect disease).

}

\keyword{methods}

