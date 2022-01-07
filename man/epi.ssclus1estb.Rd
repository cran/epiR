\name{epi.ssclus1estb}

\alias{epi.ssclus1estb}

\title{
Sample size to estimate a binary outcome using one-stage cluster sampling
}

\description{
Sample size to estimate a binary outcome using one-stage cluster sampling.
}

\usage{
epi.ssclus1estb(b, Py, epsilon, error = "relative", 
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
In many situations it is common for sampling units to be aggregated into clusters. Typical examples include individuals within households, children within classes (within schools) and cows within herds. We use the term primary sampling unit (PSU) to refer to what gets sampled first (clusters) and secondary sampling unit (SSU) to refer to what gets sampled second (individual listing units within each cluster). In this documentation the terms primary sampling unit and cluster are used interchangeably. Similarly, the terms secondary sampling unit and individual listing units are used interchangeably. 

\code{b} as a scalar integer represents the total number of individual listing units from each cluster to be sampled. If \code{b} is a vector of length two the first element represents the mean number of individual listing units to be sampled from each cluster and the second element represents the standard deviation of the number of individual listing units to be sampled from each cluster.

At least 25 PSUs (clusters) are recommended for one-stage cluster sampling designs. If less than 25 PSUs are returned by the function a warning is issued.
}


\value{
A list containing the following: 
  \item{n.psu}{the total number of primary sampling units (clusters) to be sampled for the specified level of confidence and relative error.}
  \item{n.ssu}{the total number of secondary sampling units to be sampled for the specified level of confidence and relative error.}
  \item{DEF}{the design effect.}
  \item{rho}{the intracluster correlation, as entered by the user.}  
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 258.

Machin D, Campbell MJ, Tan SB, Tan SH (2018). Sample Sizes for Clinical, Laboratory ad Epidemiological Studies, Fourth Edition. Wiley Blackwell, London, pp. 195 - 214.
}

\examples{
## EXAMPLE 1:
## An aid project has distributed cook stoves in a single province in a 
## resource-poor country. At the end of three years, the donors would like 
## to know what proportion of households are still using their donated  
## stove. A cross-sectional study is planned where villages in the province 
## will be sampled and all households (approximately 75 per village) will be 
## visited to determine whether or not the donated stove is still in use.
## A pilot study of the prevalence of stove usage in five villages 
## showed that 0.46 of householders were still using their stove. The 
## intracluster correlation for a study of this type is unknown, but thought
## to be relatively high (rho = 0.20).

# If the donor wanted to be 90\% confident that the survey estimate of stove
## usage was within 10\% of the true population value, how many villages 
## (i.e., clusters) would need to be sampled?

epi.ssclus1estb(b = 75, Py = 0.46, epsilon = 0.10, error = "relative", 
   rho = 0.20, nfractional = FALSE, conf.level = 0.90)

## A total of 67 villages need to be sampled to meet the specifications 
## of this study.

## Now imagine the situation where the number of households per village 
## varies. We are told that the average number of households per village is
## 75 with the 0.025 quartile 40 households and the 0.975 quartile 180 
## households. The expected standard deviation of the number of households
## per village is (180 - 40) / 4 = 35. How many villages need to be sampled?

epi.ssclus1estb(b = c(75,35), Py = 0.46, epsilon = 0.10, error = "relative", 
   rho = 0.20, nfractional = FALSE, conf.level = 0.90)

## A total of 81 villages need to be sampled to meet the specifications 
## of this study.

}

\keyword{univar}

