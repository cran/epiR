\name{epi.simplesize}

\alias{epi.simplesize}

\title{
Sample size under under simple random sampling
}

\description{
Estimates the required sample size under under simple random sampling. 
}

\usage{
epi.simplesize(N = 1E+06, Vsq, Py, epsilon.r, method = "mean", 
   conf.level = 0.95)
}

\arguments{
  \item{N}{scalar, representing the population size.}
  \item{Vsq}{scalar, if method is \code{total} or \code{mean} this is the relative variance of the variable to be estimated (i.e. \code{var/mean^2}).}
  \item{Py}{scalar, if method is \code{proportion} this is an estimate of the unknown population proportion.}
  \item{epsilon.r}{the maximum relative difference between our estimate and the unknown population value.}
  \item{method}{a character string indicating the method to be used. Options are \code{total}, \code{mean}, or \code{proportion}.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
}

\value{
Returns an integer defining the size of the sample is required.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 70 - 75.

Scheaffer RL, Mendenhall W, Lyman Ott R (1996). Elementary Survey Sampling. Duxbury Press, New York, pp. 95.

Otte J, Gumm I (1997). Intra-cluster correlation coefficients of 20 infections calculated from the results of cluster-sample surveys. Preventive Veterinary Medicine 31: 147 - 150.

}

\author{
}

\note{
If the calculated sample size is greater than 10\% of the population, an adjusted sample size is returned.

\code{epsilon.r} defines the maximum relative difference between our estimate and the unknown population value. The sample estimate should not differ in absolute value from the true unknown population parameter \code{d} by more than \code{epsilon.r * d}.
}

\seealso{
}
\examples{
## EXAMPLE 1
## A city contains 20 neighbourhood health clinics and it is desired to take a 
## sample of clinics to estimate the total number of persons from all these 
## clinics who have been given, during the past 12 month period, prescriptions 
## for a recently approved antidepressant. If we assume that the average number 
## of people seen at these clinics is 1500 per year with the standard deviation 
## equal to 300, and that approximately 5\% of patients (regardless of clinic) 
## are given this drug, how many clinics need to be sampled to yield an estimate 
## that is within 20\% of the true population value?

pmean <- 1500 * 0.05; pvar <- (300 * 0.05)^2
epi.simplesize(N = 20, Vsq = (pvar / pmean^2), Py = NA, epsilon.r = 0.20, 
   method = "total", conf.level = 0.95)

## Three clinics need to be sampled to meet the survey requirements. 

## EXAMPLE 2
## We want to estimate the mean bodyweight of deer on a farm. There are 278
## animals present. We anticipate the mean body weight to be around 200 kg
## and the standard deviation of body weight to be 30 kg. We would like to
## be 95\% certain that our estimate is within 10 kg of the true mean. How
## many deer should be sampled?

epi.simplesize(N = 278, Vsq = 30^2 / 200^2, Py = NA, epsilon.r = 10/200, 
   method = "mean", conf.level = 0.95)

## A total of 28 deer need to be sampled to meet the survey requirements.

## EXAMPLE 3
## We want to estimate the seroprevalence of Brucella abortus in a population 
## of cattle. An estimate of the unknown prevalence of B. abortus in this 
## population is 0.15. We would like to be 95\% certain that our estimate is 
## within 20\% of the true proportion of the population that is seropositive 
## to B. abortus. Calculate the required sample size.

n.crude <- epi.simplesize(N = 1E+06, Vsq = NA, Py = 0.15, epsilon.r = 0.20,
   method = "proportion", conf.level = 0.95)
n.crude

## A total of 544 cattle need to be sampled to meet the survey requirements.

## EXAMPLE 3 (continued)
## Being seropositive to brucellosis is likely to cluster within herds.
## Otte and Gumm (1997) cite the intraclass correlation coefficient of
## Brucella abortus to be in the order of 0.09. Adjust the sample size
## estimate to account for clustering at the herd level. Assume that, on
## average, herds in your area of interest are comprised of 100 animals.

## rho = (design - 1) / (nbar - 1)
## D <- rho * (nbar - 1) + 1

## Above, rho equals the intracless correlation coefficient and nbar equals
## the average number of individuals per cluster.

rho <- 0.09; nbar <- 100
D <- rho * (nbar - 1) + 1

n.adj <- ceiling(n.crude * D)
n.adj

## After accounting for the presence of clustering at the herd level we
## estimate that a total of 5392 cattle need to be sampled to meet
## the survey requirements.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
