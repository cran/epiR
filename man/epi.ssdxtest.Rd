\name{epi.ssdxtest}

\alias{epi.ssdxtest}

\title{
Sample size to validate a diagnostic test in the absence of a gold standard
}

\description{
Sample size to validate a diagnostic test in the absence of a gold standard. 
}

\usage{
epi.ssdxtest(pi, se, sp, epsilon.api, epsilon.ase, epsilon.asp, epsilon.asesp, 
    r = 1, nfractional = FALSE, verbose = FALSE, conf.level = 0.95)
}

\arguments{
  \item{pi}{vector of length two, the expected prevalence of the outcome of interest in the two populations (0 to 1), respectively.}
  \item{se}{vector of length two, the expected diagnostic sensitivity of the first and second test (0 to 1), respectively.}
  \item{sp}{vector of length two, the expected diagnostic specificity of the first and second test (0 to 1), respectively.}
  \item{epsilon.api}{vector of length two, the absolute error for the prevalence of the outcome of interest in the first and second study populations.}
  \item{epsilon.ase}{vector of length two, the absolute error for the sensitivity estimate of the first and second test.}
  \item{epsilon.asp}{vector of length two, the absolute error for the specificity estimate of the first and second test.}
  \item{epsilon.asesp}{vector of length two, the absolute error for the difference in the two sensitivity and specificity estimates.}
  \item{r}{scalar, the required number in the second population divided by the number in the first population.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{verbose}{logical, return sample size estimates for \code{se}, \code{sp}, and \code{pi}?}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
Hui and Walter (1980) describe an approach for estimating the sensitivity and specificity of a diagnostic test in the absence of a gold standard. Their method involves testing individuals from two populations with two conditionally independent diagnostic tests (neither of which is a gold standard). With such data, all six parameters of interest (two sensitivities, two specificities and two prevalences) can be estimated since there are six degrees of freedom available. The methodology in this function follows the sample size calculation methods described by Georgiadis et al (2005).

In their paper Georgiadis et al (2005) parameterise the uncertainty in the prevalence, sensitivity and specificity estimates in terms of the width of the confidence interval. For consistency with the other sample size calculation functions in \pkg{epiR} the amount of uncertainty in the prevalence, sensitivity and specificity is parameterised in absolute terms. Using this approach, if we set \code{se = c(0.80,0.90)} and \code{epsilon.ase = c(0.05,0.10)} the number of subjects to return an estimate of \code{se1} that is between \code{0.75} and \code{0.85} and \code{se2} that is between \code{0.80} and \code{1.0} will be returned.   
}

\value{
When \code{verbose = TRUE} a data frame listing the required sample size to estimate:
 
  \item{p1}{the prevalence of the outcome of interest in population 1.}
  \item{p2}{the prevalence of the outcome of interest in population 2.}
  \item{se1}{the sensitivity of the first diagnostic test.}
  \item{se2}{the sensitivity of the second diagnostic test.} 
  \item{sp1}{the specificity of the first diagnostic test.}
  \item{sp2}{the specificity of the second diagnostic test.} 
  \item{se1.se2}{the difference in the sensitivities of the two diagnostic tests.}
  \item{sp1.sp2}{the difference in the specificities of the two diagnostic tests.}

When \code{verbose = FALSE} a data frame listing the maximum of the sample size estimates listed when \code{verbose = TRUE}. 
}

\references{
Georgiadis M, Johnson W, Gardner I (2005) Sample size determination for estimation of the accuracy of two conditionally independent tests in the absence of a gold standard. Preventive Veterinary Medicine 71, 1 - 10. DOI: 10.1016/j.prevetmed.2005.04.004.

Hui SL, Walter SD (1980) Estimating the error rates of diagnostic tests. Biometrics 36, 167 - 171.

Nielsen SS, Gronbaek C, Agger JF, Houe H (2002) Maximum-likelihood estimation of sensitivity and specificity of ELISAs and faecal culture for diagnosis of paratuberculosis. Preventive Veterinary Medicine 53, 191 - 204. DOI: 10.1016/s0167-5877(01)00280-x.
}

\examples{
## EXAMPLE 1 (from Georgiadis et al. 2005, pp. 5):
## Nielsen et al. (2002) present data from the evaluation of a milk 
## antibody ELISA and faecal culture for the diagnosis of Mycobacterium avium 
## subsp. paratuberculosis infection in cattle. Because the ELISA detects 
## antibodies and culture is based on isolation of the bacterium in faeces 
## we can reasonably assume that the two tests are conditionally independent.

## How many cattle need to be sampled if we wanted to be 95\% certain that
## our estimate of diagnostic sensitivity and specificity of the two tests 
## is within 0.05 of the true population value assuming the number sampled 
## in the second population divided by the number sampled in the first 
## population is 0.817? The prevalence of Mycobacterium avium subsp. 
## paratuberculosis is thought to be 0.093 and 0.204, respectively. Assume 
## the sensitivity of the the ELISA and faecal culture is 0.349 and 0.534, 
## respectively. Assume the specificity of the ELISA and faecal culture is 
## 0.995 and 0.894, respectively.

epi.ssdxtest(pi = c(0.093,0.204), se = c(0.349,0.534), sp = c(0.995,0.894), 
   epsilon.api = c(0.05,0.05), epsilon.ase = c(0.05,0.05), 
   epsilon.asp = c(0.05,0.05), epsilon.asesp = c(0.05,0.05), 
    r = 0.817, nfractional = FALSE, verbose = FALSE, conf.level = 0.95)

## A total of 63,887 cattle need to be sampled (35,161 from population 1 and
## 28,726 from population 2) to meet the requirements of the study.

}

\keyword{univar}
