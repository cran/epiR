\name{epi.ssdxsesp}

\alias{epi.ssdxsesp}

\title{
Sample size to estimate the sensitivity or specificity of a diagnostic test
}

\description{
Sample size to estimate the sensitivity or specificity of a diagnostic test. 
}

\usage{
epi.ssdxsesp(se, sp, Py, epsilon, error = "relative", 
   nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{se}{scalar number, the prior estimate of diagnostic test sensitivity (0 to 1).}
  \item{sp}{scalar number, the prior estimate of diagnostic test specificity (0 to 1).}
  \item{Py}{scalar number, an estimate of the prevalence of the outcome in the study population.}
  \item{epsilon}{scalar number, the maximum difference between the estimate and the unknown population value expressed in absolute or relative terms.}
  \item{error}{character string. Options are \code{absolute} for absolute error and \code{relative} for relative error.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar number, the level of confidence in the computed result.}
}

\value{
A data frame listing: \code{se}, the prior estimate of diagnostic senstivity as entered by the user, \code{sp}, the prior estimate of diagnostic specificity as entered by the user, \code{Py} the estimate of the prevalence of the outcome in the population as entered by the user, \code{se.n} the sample size to estimate diagnostic sensitivity, \code{sp.n} the sample size to estimate diagnostic specificity and \code{total.n} the total sample size for the study (which is the maximum of \code{se.n} and \code{sp.n}).
}

\references{
Buderer, NM (1996). Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. Academic Emergency Medicine, 3: 895 - 900. DOI: 10.1111/j.1553-2712.1996.tb03538.x

Hajian-Tilaki K (2014). Sample size estimation in diagnostic test studies of biomedical informatics. Journal of Biomedical Informatics 48: 193 - 204. DOI: 10.1016/j.jbi.2014.02.013.

van Hal S, Lahra M, Whiley D (2025) Determining targets for ceftriaxone resistance detection in Neisseria gonorrhoeae based on genotype-phenotype data: An observational study. Lancet Microbe: 101199. DOI: 10.1016/j.lanmic.2025.101199.

}

\note{
The sample size calculation method implemented in this function follows the approach described by Buderer (1996) and Hajian-Tilaki (2014). 
}

\examples{
## EXAMPLE 1 (from Hajian-Tilaki 2014, p 195):
## A new diagnostic test has been developed and we'd like to conduct a study
## to determine its diagnostic sensitivity which we believe should be in the
## order of 0.80. How many subjects should be enrolled if the prevalence of 
## the disease outcome of interest is 0.10 and we'd like to be 95\% confident 
## that our estimate of sensitivity is within 0.07 of the true population 
## value? For this example we want to estimate diagnostic sensitivity (only)
## so the estimated performance measure for specificity is set to NA.

epi.ssdxsesp(se = 0.80, sp = NA, Py = 0.10, epsilon = 0.07, 
   error = "absolute", nfractional = FALSE, conf.level = 0.95)

## A total of 1255 subjects need to be enrolled to meet the requirements of the 
## study.


## EXAMPLE 2 (from van Hal et al. 2025):
## van Hal et al (2025) determined targets for ceftriaxone resistance detection 
## in Neisseria gonorrhoeae based on genotype-phenotype. To estimate 
## diagnostic testing performance with high precision [sic], a sample size 
## of 8295 distinct MIC-geonotype pairings was required, to achieve a 
## sensitivity and specificity of 95%, based on a ceftriaxone resistance 
## prevalence of 0.22\% --- a value reflecting current rates in high-income 
## settings with robust surveillance systems.

## Code to reproduce the sample size calculations reported in the van Hal 
## et al. (2025) paper as follows:

epi.ssdxsesp(se = 0.95, sp = 0.95, Py = 0.0022, epsilon = 0.10, 
   error = "absolute", nfractional = FALSE, conf.level = 0.95)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
