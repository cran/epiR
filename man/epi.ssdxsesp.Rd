\name{epi.ssdxsesp}

\alias{epi.ssdxsesp}

\title{
Sample size to estimate the sensitivity or specificity of a diagnostic test
}

\description{
Sample size to estimate the sensitivity or specificity of a diagnostic test. 
}

\usage{
epi.ssdxsesp(test, type = "se", Py, epsilon, error = "relative", 
   nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{test}{scalar number, the prior estimate of diagnostic test performance (0 to 1).}
  \item{type}{character string. Options are \code{se} to estimate a sample size to determine diagnostic sensitivity and \code{sp} to estimate a sample size to determine diagnostic specificity.}
  \item{Py}{scalar number, an estimate of the prevalence of the outcome in the study population.}
  \item{epsilon}{scalar number, the maximum difference between the estimate and the unknown population value expressed in absolute or relative terms.}
  \item{error}{character string. Options are \code{absolute} for absolute error and \code{relative} for relative error.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar number, the level of confidence in the computed result.}
}

\value{
Returns an integer defining the required sample size.
}

\references{
Hajian-Tilaki K (2014). Sample size estimation in diagnostic test studies of biomedical informatics. Journal of Biomedical Informatics 48: 193 - 204. DOI: 10.1016/j.jbi.2014.02.013.
}

\note{
The sample size calculation method implemented in this function follows the approach described by Hajian-Tilaki (2014). 
}

\examples{
## EXAMPLE 1 (from Hajian-Tilaki 2014, p 195):
## A new diagnostic test has been developed and we'd like to conduct a study
## to determine its diagnostic sensitivity which we believe should be in the
## order of 0.80. How many subjects should be enrolled if the prevalence of 
## the disease outcome of interest is 0.10 and we'd like to be 95\% confident 
## that our estimate of sensitivity is within 0.07 of the true population 
## value?

epi.ssdxsesp(test = 0.80, type = "se", Py = 0.10, epsilon = 0.07, 
   error = "absolute", nfractional = FALSE, conf.level = 0.95)

## A total of 1255 subjects need to be enrolled to meet the requirements of the 
## study.
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
