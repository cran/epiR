\name{epi.ssclus2estb}

\alias{epi.ssclus2estb}

\title{Number of clusters to be sampled to estimate a binary outcome using two-stage cluster sampling
}

\description{
Number of clusters to be sampled to estimate a binary outcome using two-stage cluster sampling.
}

\usage{
epi.ssclus2estb(b, Py, epsilon, error = "relative", 
   rho, nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{b}{scalar integer or vector of length two, the number of individual listing units in each cluster to be sampled. See details, below.}
  \item{Py}{scalar number, an estimate of the unknown population proportion.}
  \item{epsilon}{scalar number, the maximum difference between the estimate and the unknown population value expressed in absolute or relative terms.}
  \item{error}{character string. Options are \code{absolute} for absolute error and \code{relative} for relative error.}
  \item{rho}{scalar number, the intracluster correlation.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
\code{b} as a scalar integer represents the total number of individual listing units from each cluster to be sampled. If \code{b} is a vector of length two the first element represents the mean number of individual listing units to be sampled from each cluster and the second element represents the standard deviation of the number of individual listing units to be sampled from each cluster.

The methodology used in this function follows closely the approach described by Bennett et al. (1991). At least 25 primary sampling units are recommended for two-stage cluster sampling designs. If less than 25 clusters are returned by the function a warning is issued.

As a rule of thumb, around 30 clusters will provide good estimates of the true population value with an acceptable level of precision (Binkin et al. 1992) when: (1) the true population value is between 10\% and 90\%; and (2) the desired absolute error is around 5\%. For a fixed number of individuals selected per cluster (e.g. 10 individuals per cluster or 30 individuals per cluster), collecting information on more than 30 clusters can improve the precision of the final population estimate, however, beyond around 60 clusters the improvement in precision is minimal.
}

\value{
A list containing the following: 
  \item{n.psu}{the total number of primary sampling units (clusters) to be sampled for the specified level of confidence and relative error.}
  \item{n.ssu}{the total number of secondary sampling units to be sampled for the specified level of confidence and relative error.}
  \item{DEF}{the design effect.}
  \item{rho}{the intracluster correlation, as entered by the user.}  
}

\references{
Bennett S, Woods T, Liyanage W, Smith D (1991). A simplified general method for cluster-sample surveys of health in developing countries. World Health Statistics Quarterly 44: 98 - 106.

Binkin N, Sullivan K, Staehling N, Nieburg P (1992). Rapid nutrition surveys: How many clusters are enough? Disasters 16: 97 - 103.

Machin D, Campbell MJ, Tan SB, Tan SH (2018). Sample Sizes for Clinical, Laboratory ad Epidemiological Studies, Fourth Edition. Wiley Blackwell, London, pp. 195 - 214.
}

\examples{
## EXAMPLE 1 (from Bennett et al. 1991 p 102):
## We intend to conduct a cross-sectional study to determine the prevalence 
## of disease X in a given country. The expected prevalence of disease is 
## thought to be around 20\%. Previous studies report an intracluster 
## correlation for this disease to be 0.02. Suppose that we want to be 95\%
## certain that our estimate of the prevalence of disease is within 5\% of 
## the true population value and that we intend to sample 20 individuals per
## cluster. How many clusters should be sampled to meet the requirements of the
## study?

epi.ssclus2estb(b = 20, Py = 0.20, epsilon = 0.05, error = "absolute", 
   rho = 0.02, nfractional = FALSE, conf.level = 0.95)

## A total of 17 clusters need to be sampled to meet the specifications 
## of this study. epi.ssclus2estb returns a warning message that the number of 
## clusters is less than 25.
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
