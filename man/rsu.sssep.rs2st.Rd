\name{rsu.sssep.rs2st}

\alias{rsu.sssep.rs2st}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming two-stage sampling
}

\description{
Calculates the required sample size to achieve a desired surveillance system sensitivity assuming two-stage sampling (sampling of clusters and sampling of units within clusters), imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.sssep.rs2st(H = NA, N = NA, pstar.c, se.c, pstar.u, se.u, se.p)
}

\arguments{
  \item{H}{scalar, integer representing the total number of clusters in the population. Use \code{NA} if unknown.}
  \item{N}{vector, integer representing the number of units within each cluster. Use \code{NA} if unknown.}
  \item{pstar.c}{scalar, numeric (0 to 1) representing the cluster level design prevalence.}
  \item{se.c}{scalar, numeric (0 to 1) representing the required cluster level sensitivity.}
  \item{pstar.u}{scalar, numeric (0 to 1) representing the surveillance unit level design prevalence.}
  \item{se.u}{scalar (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}
  \item{se.p}{scalar (0 to 1) representing the desired surveillance system (population-level) sensitivity.}
}

\value{
A list comprised of two data frames: \code{clusters} and \code{units}. Data frame \code{clusters} lists:
  
\item{H}{the total number of clusters in the population, as entered by the user.}
\item{nsample}{the number of clusters to be sampled.}

Data frame \code{units} lists:

\item{N}{the number of units within each cluster, as entered by the user.}
\item{nsample}{the number of units to be sampled.}
}

\references{
Cameron A, Baldock C (1998). A new probability formula for surveys to substantiate freedom from disease. Preventive Veterinary Medicine 34: 1 - 17.

Cameron A (1999). Survey Toolbox for Livestock Diseases --- A practical manual and software package for active surveillance of livestock diseases in developing countries. Australian Centre for International Agricultural Research, Canberra, Australia.

MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}


\examples{
## EXAMPLE 1:
## Sampling is to be carried out to support a claim that a country is free
## of bovine brucellosis. We are not certain of the total number of herds
## in the country and we are not certain of the number of cows within each
## herd. 

## The design prevalence for this study is set to 0.01 at the herd level and
## if a herd is positive for brucellosis the individual animal level 
## design prevalence is set to 0.10. The sensitivity of the diagnostic 
## test to be used is 0.95. 

## How many herds and how many animals from within each herd 
## need to be sampled to be 95\% confident of detecting disease at the 
## herd and individual animal level? 

rsu.sssep.rs2st(H = NA, N = NA, pstar.c = 0.01, se.c = 0.95, 
   pstar.u = 0.10, se.u = 0.95, se.p = 0.95)
   
## A total of 314 herds need to be sampled, 31 cows from each herd.


## EXAMPLE 2:
## Now lets say we know that there are 500 cattle herds in the country and
## we have the results of a recent livestock census providing counts of the 
## number of cattle in each herd. How many herds and how many animals from
## within each herd need to be sampled to be 95\% confident of detecting
## disease at the herd and individual animal level?

# Generate a vector of herd sizes. The minimum herd size is 25.

set.seed(1234)
hsize <- ceiling(rlnorm(n = 500, meanlog = 1.5, sdlog = 2)) + 25 

nsample <- rsu.sssep.rs2st(H = 500, N = hsize, pstar.c = 0.01, se.c = 0.95,  
   pstar.u = 0.10, se.u = 0.95, se.p = 0.95)
  
nsample$clusters
head(nsample$units)
 
## A total of 238 of the 500 herds need to be tested. The number of animals 
## to sample from the first herd (comprised of 26 animals) is 18.
}

\keyword{univar}
