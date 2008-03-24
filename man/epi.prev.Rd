\name{epi.prev}

\alias{epi.prev}

\title{Estimate true prevalence
}
\description{
Computes the true prevalence of a disease in a population on the basis of an imperfect test.
}

\usage{
epi.prev(pos, tested, N, se, sp, conf.level = 0.95)
}

\arguments{
  \item{pos}{the number of positives.}
  \item{tested}{the number tested.}
  \item{N}{the size of the population of interest.}
  \item{se}{test sensitivity (0 - 1).}
  \item{sp}{test specificity (0 - 1).}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
}

\value{
A data frame with four elements: \code{crude}: the apparent prevalence, \code{est}: the estimated true prevalence, \code{lower}: the lower bounds of the confidence interval of the estimated prevalence, and \code{upper}: the upper bounds of the confidence interval of the estimated prevalence.
}

\references{
Rogan W, Gladen B (1978). Estimating prevalence from results of a screening test. American Journal of Epidemiology 107: 71 - 76.

Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, pp. 130 - 143. 

}

\author{
}

\note{
This function uses the apparent prevalence, test sensitivity and test specificity to estimate true prevalence (after Rogan and Gladen, 1978). Confidence limits for the true prevalence estimate are calculated using Wilson's approximation (Rothman 2002, page 132).  

}

\seealso{
}

\examples{
## A simple random sample of 150 cows from a herd of 2560 is taken.
## Each cow is given a screening test for brucellosis which has a 
## sensitivity of 96\% and a specificity of 89\%. Of the 150 cows tested
## 23 were positive to the screening test. What is the estimated prevalence 
## of brucellosis in this herd (and its 95\% confidence interval)?

epi.prev(23, 150, 2560, se = 0.96, sp = 0.89, conf.level = 0.05)

## The estimated true prevalence of brucellosis in this herd is 5.1 cases per 
## 100 cows (95\% CI 0 -- 12 cases per 100 cows).
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
