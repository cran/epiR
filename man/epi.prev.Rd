\name{epi.prev}

\alias{epi.prev}

\title{Estimate true prevalence
}
\description{
Computes the true prevalence of a disease in a population on the basis of an imperfect test.
}

\usage{
epi.prev(pos, tested, se, sp, conf.level = 0.95)
}

\arguments{
  \item{pos}{the number of positives.}
  \item{tested}{the number tested.}
  \item{se}{test sensitivity (0 - 1).}
  \item{sp}{test specificity (0 - 1).}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Exact binomial confidence limits are calculated for apparent prevalence (see Collett 1999 for details). 
}

\value{
A list containing the following:
   \item{ap}{the point estimate of apparent prevalence, the standard error of the apparent prevalence, and the lower and upper bounds of the confidence interval around the apparent prevalence estimate.}
   \item{tp}{the point estimate of the true prevalence, the standard error of the true prevalence, and the lower and upper bounds of the confidence interval around the true prevalence estimate.}
}

\references{
Abel U (1993). DieBewertung Diagnostischer Tests. Hippokrates, Stuttgart.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, p. 24. 

Gardener IA, Greiner M (1999). Advanced Methods for Test Validation and Interpretation in Veterinary Medicince. Freie Universitat Berlin, ISBN 3-929619-22-9; 80 pp.

Messam L, Branscum A, Collins M, Gardner I (2008) Frequentist and Bayesian approaches to prevalence estimation using examples from Johne's disease. Animal Health Research Reviews 9: 1 - 23.

Rogan W, Gladen B (1978). Estimating prevalence from results of a screening test. American Journal of Epidemiology 107: 71 - 76.
}

\author{
}

\note{
This function uses apparent prevalence, test sensitivity and test specificity to estimate true prevalence (after Rogan and Gladen, 1978). The standard error of the Rogan Gladen true prevalence estimate is based on Abel (1993) and discussed in Messam et al. (2008). It is assumed that test sensitivity and specificity are known with certainty.

The Rogan Gladen true prevalence estimate is unreliable for small sample sizes and when true prevalence is believed to be close to zero. The algorithm implemented here makes no correction to the Rogan Gladen estimate of true prevalence if it is less than zero or greater than one (simply to remind the user that it provides unreliable estimates of true prevalence under these conditions). In this situation one is advised to adopt a Bayesian approach to true prevalence estimation. See Messam et al. (2008) for a very readable introduction.
}

\seealso{
}

\examples{
## A simple random sample of 150 cows from a herd of 2560 is taken.
## Each cow is given a screening test for brucellosis which has a 
## sensitivity of 96\% and a specificity of 89\%. Of the 150 cows tested
## 23 were positive to the screening test. What is the estimated prevalence 
## of brucellosis in this herd (and its 95\% confidence interval)?

epi.prev(pos = 23, tested = 150, se = 0.96, sp = 0.89, conf.level = 0.95)

## The estimated true prevalence of brucellosis in this herd is 5.1 cases per 
## 100 cows (95\% CI 0 -- 12 cases per 100 cows).
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
