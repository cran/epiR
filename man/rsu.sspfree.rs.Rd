\name{rsu.sspfree.rs}

\alias{rsu.sspfree.rs}

\title{
Sample size to achieve a desired probability of disease freedom assuming representative sampling
}

\description{
Calculates the required sample size to achieve a desired (posterior) probability of disease freedom assuming representative sampling, imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.sspfree.rs(N = NA, prior, p.intro, pstar, pfree, se.u)
}

\arguments{
\item{N}{scalar integer or vector of integers the same length as \code{n}, representing the population size. Use \code{NA} if unknown.}
\item{prior}{scalar probability (0 to 1), representing the prior probability that the population is free of disease.}
\item{p.intro}{scalar or vector of the same length as \code{pfree}, representing the probability of disease introduction during the next time period.}
\item{pstar}{scalar numeric or vector of numbers the same length as \code{pfree} representing the design prevalence.}
\item{pfree}{scalar numeric or vector of numbers the same length as \code{pfree} representing the desired probability of disease freedom.}
\item{se.u}{scalar (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}
}

\value{
A list comprised of three elements: 

\item{n}{a vector listing the required sample sizes.}
\item{sep}{a vector listing the population sensitivity estimates.}
\item{adj.prior}{a vector listing the adjusted priors.}
}

\references{
Martin P, Cameron A, Greiner M (2007). Demonstrating freedom from disease using multiple complex data sources 1: A new methodology based on scenario trees. Preventive Veterinary Medicine 79: 71 - 97.

Martin P, Cameron A, Barfod K, Sergeant E, Greiner M (2007). Demonstrating freedom from disease using multiple complex data sources 2: Case study - Classical swine fever in Denmark. Preventive Veterinary Medicine 79: 98 - 115.
}

\note{

This function returns the sample size to achieve a desired (posterior) probability of disease freedom. Function \code{\link{rsu.sssep.rs}} returns the sample size to achieve a desired surveillance system sensitivity.
}

\examples{
## EXAMPLE 1:
## Prior surveillance activities and expert opinion lead you to believe that 
## there's a 75\% chance that your country is free of disease X. To confirm
## your country's disease freedom status you intend to use a test at the herd 
## level which has a diagnostic sensitivity of 0.95. The probability of 
## disease introduction during the time period of interest is relatively
## low, say 0.01. How many herds need to be sampled to be 95\% confident
## that the country is free of disease X assuming a design prevalence of 
## 0.01?

rsu.sspfree.rs(N = NA, prior = 0.75, p.intro = 0.01, pstar = 0.01, 
   pfree = 0.95, se.u = 0.95)

## A total of 198 herds need to be sampled to meet the requirements of the
## study.

}
\keyword{methods}

