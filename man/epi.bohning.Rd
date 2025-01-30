\name{epi.bohning}

\alias{epi.bohning}

\title{Bohning's test for overdispersion of Poisson data}

\description{
A test for overdispersion of Poisson data.
}

\usage{
epi.bohning(obs, exp, alpha = 0.05)
}

\arguments{
  \item{obs}{the observed number of cases in each area.}
  \item{exp}{the expected number of cases in each area.}
  \item{alpha}{alpha level to be used for the test of significance. Must be a single number between 0 and 1.}
}

\value{
A data frame with two elements: \code{test.statistic}, Bohning's test statistic and \code{p.value} the associated P-value.
}

\references{
Bohning D (2000). Computer-assisted Analysis of Mixtures and Applications. Chapman and Hall, Boca Raton.

Ugarte MD, Ibanez B, Militino AF (2006). Modelling risks in disease mapping. Statistical Methods in Medical Research 15: 21 - 35. 
}

\examples{
## EXAMPLE 1:
## Lip cancer in Scotland, 1975 to 1980:
data(epi.SClip)
obs <- epi.SClip$cases
pop <- epi.SClip$population

## Calculate the expected number of lip cancer cases per district:
exp <- (sum(obs) / sum(pop)) * pop

## Bohning's test:
epi.bohning(obs, exp, alpha = 0.05)

## The test statistic was 53.33. The associated P value was <0.01. We reject 
## the null hypothesis of no over dispersion and accept the null hypothesis
## concluding that the lip cancer data are over dispersed.

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
