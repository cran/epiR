\name{rsu.sep.pass}

\alias{rsu.sep.pass}

\title{
Surveillance system sensitivity assuming passive surveillance and representative sampling within clusters
}

\description{
Calculates the surveillance system (population-level) sensitivity for detection of disease for a passive surveillance system assuming comprehensive population coverage and sampling of clinical cases within diseased clusters.}

\usage{
rsu.sep.pass(N, n, step.p, pstar.c, p.inf.u, se.u)
}

\arguments{
\item{N}{scalar or vector of length equal to the number of rows in \code{step.p} representing the population size.}
\item{n}{scalar or vector of length equal to the number of rows in \code{step.p} representing the number of units tested per cluster.}
\item{step.p}{vector or matrix of detection probabilities (0 to 1) for each step in the detection process. If a vector each value represents a step probability for a single calculation. If a matrix, columns
are step probabilities and rows are simulation iterations.}
\item{pstar.c}{scalar (0 to 1) or vector of length equal to the number of rows in \code{step.p} representing the cluster-level design prevalence.}
\item{p.inf.u}{scalar (0 to 1) or vector of length equal to the number of rows in \code{step.p} representing the probability of disease in sampled and tested units. This is equivalent to the positive predictive value for a given prior probability of infection.}
\item{se.u}{scalar (0 to 1) or vector of length equal to the number of rows in \code{step.p}, representing the unit sensitivity.}
}

\value{
A list comprised of two elements: 

\item{se.p}{scalar or vector, the estimated surveillance system (population-level) sensitivity of detection.}
\item{se.c}{scalar or vector, the estimated cluster-level sensitivity of detection.}

If \code{step.p} is a vector, scalars are returned. If \code{step.p} is a matrix, values are vectors of length equal to the number of rows in \code{step.p}.
}

\references{
Lyngstad T, Hellberg H, Viljugrein H, Bang Jensen B, Brun E, Sergeant E, Tavornpanich S (2016). Routine clinical inspections in Norwegian marine salmonid sites: A key role in surveillance for freedom from pathogenic viral haemorrhagic septicaemia (VHS). Preventive Veterinary Medicine 124: 85 - 95. DOI:  10.1016/j.prevetmed.2015.12.008.
}


\examples{
## EXAMPLE 1:
## A passive surveillance system for disease X operates in your country.
## There are four steps to the diagnostic cascade with detection probabilities 
## for each process of 0.10, 0.20, 0.90 and 0.99, respectively. Assuming the 
## probability that a unit actually has disease if it is submitted for 
## testing is 0.98, the sensitivity of the diagnostic test used at the unit 
## level is 0.90, the population is comprised of 1000 clusters (herds), 
## five animals from each cluster (herd) are tested and the cluster-level 
## design prevalence is 0.01, what is the sensitivity of disease detection 
## at the cluster (herd) and population level?

rsu.sep.pass(N = 1000, n = 5, step.p = c(0.10,0.20,0.90,0.99), 
   pstar.c = 0.01, p.inf.u = 0.98, se.u = 0.90)

## The sensitivity of disease detection at the cluster (herd) level is 0.018.
## The sensitivity of disease detection at the population level is 0.16.

}
\keyword{methods}

