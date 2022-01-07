\name{rsu.spp.rs}

\alias{rsu.spp.rs}

\title{
Surveillance system specificity assuming representative sampling
}

\description{
Calculates surveillance system (population level) specificity assuming representative sampling and imperfect test specificity.
}

\usage{
rsu.spp.rs(N, n, c = 1, sp.u)
}

\arguments{
\item{N}{scalar or vector of the same length as that vector of \code{n} defining the [cluster] population size. Use \code{NA} if the size of the population not known, or for a more general application see details, below.}
\item{n}{scalar or vector defining the sample size.}
\item{c}{scalar or vector of the same length as that vector of \code{n} defining the cut-point number of positives to classify a cluster as positive, if the number of positive samples is less than \code{c} the cluster is declared is negative, if the number of positive samples is greater than \code{c} the cluster is declared positive.}
\item{sp.u}{scalar (0 to 1) or vector of same length as \code{n}, the specificity of the diagnostic test at the surveillance unit level.}
}

\details{
This function calculates population specificity using the hypergeometric distribution if \code{N} and \code{c} are provided and the binomial distribution otherwise.

If \code{N} is provided the number of false positives is fixed, based on \code{N} and test specificity \code{sp.u}. This implies that test specificity is a fixed individual-level characteristic (e.g., due to specific cross-reacting infection). If \code{N} is not supplied, cluster (e.g., herd) specificity is a random binomial function based only on the number of samples and test specificity (i.e., specificity is a function of the test and independent of individual characteristics).  
}

\value{
A vector of population specificity estimates.
}

\references{
Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\examples{
## EXAMPLE 1:
## Calculate the surveillance system specificity (i.e., the probability that 
## an uninfected population will be correctly identified as negative) if 30 
## surveillance units have been tested from a population of 150 using a 
## diagnostic test with surveillance unit specificity of 0.90, using a 
## cut-point of one or more positives to consider the population positive.

## A specificity of 0.90 means that 9 out of 10 samples from disease-negative
## surveillance units will return a negative result (i.e., one of them will be
## a false positive).
 
rsu.spp.rs(N = 150, n = 30, c = 1, sp.u = 0.90)

## The surveillance system specificity is 0.03. There is a probability of 
## 0.03 that all 30 samples will be negative.   


## EXAMPLE 2:
## Now assume we set a cut-point of 6. That is, 6 or more samples have to 
## return a positive result for us to declare the population positive:

rsu.spp.rs(N = 150, n = 30, c = 6, sp.u = 0.90)

## The surveillance system specificity is 0.95.

}
\keyword{methods}

