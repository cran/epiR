\name{rsu.pfree.equ}

\alias{rsu.pfree.equ}

\title{
Equilibrium probability of disease freedom assuming representative or risk based sampling
}

\description{
Calculates the long-term equilibrium probability of disease freedom and equilibrium prior probability of freedom, after discounting for the probability that disease has been introduced into the population and assuming population sensitivity and probability of introduction are constant over time. It does not specify how long it might take to reach equilibrium.
}

\usage{
rsu.pfree.equ(se.p, p.intro)
}

\arguments{
\item{se.p}{scalar or vector, the surveillance system (population-level) sensitivity for the given time period.}
\item{p.intro}{scalar or vector of the same length as \code{sep} representing the probability of disease introduction for time period.}
}

\value{
A list comprised of two elements: 

\item{epfree}{a vector listing the equilibrium probability of disease freedom.}
\item{depfree}{a vector listing the discounted equilibrium probability of disease freedom.}
}

\examples{
## EXAMPLE 1:
## The current (ongoing) surveillance system for a given disease in your 
## country has been estimated to have a population sensitivity of 0.60 per 
## time period (one year). Assuming the probability of disease introduction 
## per unit time is 0.02, what is the eventual plateau level for confidence 
## of freedom and how long will it take to reach this level, assuming a 
# prior (starting) confidence of freedom of 0.50?

## Firstly, estimate the equilibrium (plateau) confidence of freedom:

conf.eq <- rsu.pfree.equ(se.p = 0.60, p.intro = 0.02)
conf.eq

## The equilibrium discounted probability of disease freedom is 0.986.

## Next, calculate confidence of freedom over 20 time periods for se.p = 0.60
## and p.intro = 0.02:

rval.df <- rsu.pfree.rs (se.p = rep(0.6, times = 20), 
   p.intro = rep(0.02, times = 20), prior = 0.50)
head(rval.df)

## When does the confidence of freedom first reach the equilibrium value 
## (rounded to 3 digits)? 

rsep.p <- which(rval.df$pfree >= round(conf.eq$depfree, digits = 3))
rsep.p[1]

## It takes 9 time periods (years) to reach the equilibrium level of 0.986. 


## EXAMPLE 2:
## You have been asked to design a surveillance system to detect a given 
## disease in your country. If the probability of disease introduction per 
## unit time is 0.10, what surveillance system sensitivity do you need to 
## be 95\% certain that disease is absent based on the testing carried out as 
## part of your program?

## Generate a vector of candidate surveillance system sensitivity estimates 
## from 0.1 to 0.99:

se.p <- seq(from = 0.10, to = 0.99, by = 0.01)

## Calculate the probability of disease freedom for each of the candidate 
## surveillance system sensitivity estimates:

rval.df <- rsu.pfree.equ(se.p = se.p, p.intro = 0.10)
rval.df <- data.frame(se.p = se.p, depfree = rval.df$depfree)
head(rval.df)

## Which of the surveillance system sensitivity estimates returns a 
## probability of freedom greater than 0.95?

rsep.p <- rval.df$se.p[rval.df$depfree > 0.95]
rsep.p[1]

## The required surveillance system sensitivity for this program is 0.69. 
## Plot the results:

## Not run: 
library(ggplot2)

ggplot(data = rval.df, aes(x = se.p, y = depfree)) +
   geom_point() + 
   geom_line() + 
   scale_x_continuous(limits = c(0,1), 
      name = "Surveillance system sensitivity") + 
   scale_y_continuous(limits = c(0,1), 
      name = "Equilibrium discounted probability of disease freedom") +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed") +
  geom_vline(aes(xintercept = rsep.p[1]), linetype = "dashed") +
   theme_bw()
   
## End(Not run)
}
\keyword{methods}

