\name{epi.popsize}

\alias{epi.popsize}

\title{
Estimate population size
}

\description{
Estimates population size on the basis of capture-recapture sampling. 
}

\usage{
epi.popsize(T1, T2, T12, conf.level = 0.95, verbose = FALSE)
}

\arguments{
  \item{T1}{an integer representing the number of individuals tested in the first round.}
  \item{T2}{an integer representing the number of individuals tested in the second round.}
  \item{T12}{an integer representing the number of individuals tested in both the first and second round.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
 \item{verbose}{logical indicating whether detailed or summary results are to be returned.}
}

\details{
}

\value{
Returns the estimated population size and an estimate of the numbers of individuals that remain untested.
}

\references{
Cannon RM, Roe RT (1982). Livestock Disease Surveys A Field Manual for Veterinarians. Australian Government Publishing Service, Canberra, pp. 34.  
}

\author{
}

\note{
}

\seealso{
}
\examples{
## EXAMPLE 1
## In a field survey 400 feral pigs are captured, marked and then released.
## On a second occassion 40 of the orignal capture are found when another 400    
## pigs are captured. Estimate the size of this feral pig population. Estimate     
## the number of feral pigs that have not been tested.

epi.popsize(T1 = 400, T2 = 400, T12 = 40, conf.level = 0.95, verbose = FALSE)

## Estimated population size: 4000 (95\% CI 3125 - 5557)
## Estimated number of untested pigs: 3240 (95\% CI 2365 - 4797) 

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
