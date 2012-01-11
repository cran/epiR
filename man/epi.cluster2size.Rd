\name{epi.cluster2size}

\alias{epi.cluster2size}

\title{Sample size under under two-stage cluster sampling
}

\description{
Returns the required number of clusters to be sampled using a two-stage cluster sampling strategy.
}

\usage{
epi.cluster2size(nbar, n, mean, var, epsilon.r, method = "mean", 
   conf.level = 0.95)
}


\arguments{
  \item{nbar}{integer, representing the total number of primary sampling units to be selected from each cluster.}
  \item{n}{vector of length two, specifying the total number of clusters in the population and the total number of primary sampling units in each cluster, respectively.}
  \item{mean}{vector of length two, specifying the mean of the variable of interest at the cluster level and the mean of the variable of interest at the primary sampling unit level, respectively.}
  \item{var}{vector of length two, specifying the variance of the variable of interest at the cluster level and the variance of the variable of interest at the primary sampling unit level, respectively.}
  \item{epsilon.r}{the maximum relative difference between our estimate and the unknown population value.}
  \item{method}{a character string indicating the method to be used. Options are \code{total}, or \code{mean}.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
In simple two-stage cluster sampling the desired number of primary sampling units to be selected from each cluster is determined on the basis of cost and on the basis of the relative sizes of the first- and second-stage variance components. Once the number of primary sampling units is fixed we might then wish to determine the total number of clusters to be selected at the first stage of sampling in order to be confident of obtaining estimates that reflect the true population value.
}

\value{
Returns an integer defining the required number of clusters to be sampled.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 292.
}

\examples{
## We intend to conduct a survey of nurse practitioners to estimate the 
## average number of patients seen by each nurse. There are five health
## centres in the study area, each with three nurses. We intend to sample
## two nurses from each health centre. We would like to be 95\% confident
## of obtaining an estimate that is within 30\% of the true population value. 
## We expect that the mean number of patients seen at the health centre level 
## is 84 (var 567) and the mean number of patients seen at the nurse 
## level is 28 (var 160). How many health centres should we sample?

nbar <- 2; n <- c(5, 3); mean <- c(84, 28); var <- c(567, 160)
epi.cluster2size(nbar, n, mean, var, epsilon.r = 0.3, method = "mean", 
   conf.level = 0.95)

## We would need a simple two-stage cluster sample of 3 health centres to 
## meet these specifications.
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
