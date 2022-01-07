\name{epi.ssclus1estc}

\alias{epi.ssclus1estc}

\title{
Sample size to estimate a continuous outcome using one-stage cluster sampling
}

\description{
Sample size to estimate a continuous outcome using one-stage cluster sampling.
}

\usage{
epi.ssclus1estc(b, N, xbar, xsigma, epsilon, error = "relative", 
   rho, nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{b}{scalar integer or vector of length two, the number of individual listing units in each cluster to be sampled. See details, below.}
  \item{N}{scalar integer, representing the total number of individual listing units in the population.}
  \item{xbar}{scalar number, the expected mean of the continuous variable to be estimated.}
  \item{xsigma}{scalar number, the expected standard deviation of the continuous variable to be estimated.}
  \item{epsilon}{scalar number, the maximum difference between the estimate and the unknown population value expressed in absolute or relative terms.}
  \item{error}{character string. Options are \code{absolute} for absolute error and \code{relative} for relative error.}
  \item{rho}{scalar number, the intracluster correlation.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar number, the level of confidence in the computed result.}
}

\value{
A list containing the following: 
  \item{n.psu}{the total number of primary sampling units (clusters) to be sampled for the specified level of confidence and relative error.}
  \item{n.ssu}{the total number of secondary sampling units to be sampled for the specified level of confidence and relative error.}
  \item{DEF}{the design effect.}
  \item{rho}{the intracluster correlation, as entered by the user.}  
}

\details{
In many situations it is common for sampling units to be aggregated into clusters. Typical examples include individuals within households, children within classes (within schools) and cows within herds. We use the term primary sampling unit (PSU) to refer to what gets sampled first (clusters) and secondary sampling unit (SSU) to refer to what gets sampled second (individual listing units within each cluster). In this documentation the terms primary sampling unit and cluster are used interchangeably. Similarly, the terms secondary sampling unit and individual listing units are used interchangeably.

\code{b} as a scalar integer represents the total number of individual listing units from each cluster to be sampled. If \code{b} is a vector of length two the first element represents the mean number of individual listing units to be sampled from each cluster and the second element represents the standard deviation of the number of individual listing units to be sampled from each cluster.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 258.

Machin D, Campbell MJ, Tan SB, Tan SH (2018). Sample Sizes for Clinical, Laboratory ad Epidemiological Studies, Fourth Edition. Wiley Blackwell, London, pp. 195 - 214.
}

\examples{
## EXAMPLE 1:
## A survey to estimate the average number of residents over 75 years of 
## age that require the services of a nurse in a given retirement village is 
## to be carried out using a one-stage cluster sampling strategy. 
## There are five housing complexes in the village with 25 residents in each.
## We expect that there might be an average of 34 residents meeting this 
## criteria (SD 5.5). We would like the estimated sample size to provide us 
## with an estimate that is within 10\% of the true value. Previous studies 
## report an intracluster correlation for the number of residents requiring the 
## services of a nurse in this retirement village housing complexes to 
## be 0.10. How many housing complexes (clusters) should be sampled?

epi.ssclus1estc(b = 25, N = 5 * 25, xbar = 34, xsigma = 5.5, 
   epsilon = 0.10, error = "relative", rho = 0.10, nfractional = FALSE, 
   conf.level = 0.95)

## A total of 2 housing complexes need to be sampled to meet the specifications 
## of this study.

}

\keyword{univar}

