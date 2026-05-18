\name{epi.pooled}

\alias{epi.pooled}

\title{
Estimate herd test characteristics when pooled sampling is used
}

\description{
We may wish to designate a group of individuals (e.g., a herd) as being either diseased or non-diseased on the basis of pooled samples. This function estimates sensitivity and specificity of this testing regime at the group (or herd) level. 
}

\usage{
epi.pooled(Se, Sp, P, m, r, dilution = FALSE)
}

\arguments{
  \item{Se}{a vector of length one defining the diagnostic sensitivity of the test at the individual sample level.}
  \item{Sp}{a vector of length one defining the diagnostic specificity of the test at the individual sample level.}
  \item{P}{scalar, the expected true prevalence of disease in the group (herd).}
  \item{m}{scalar, defining the number of individual samples to make up a pooled sample.}
  \item{r}{scalar, defining the number of pooled samples per group (or herd).}
  \item{dilution}{logical. If \code{TRUE} the function calculates a value for pooled diagnostic sensitivity dependent on the expected number of positives in the pool. If \code{FALSE} pooled diagnostic sensitivity is assumed to be the same as individual sample sensitivity.}
}

\value{
A list containing the following: 
   \item{Se}{the diagnostic sensitivity of the test at the individual sample level as entered by the user.}
   \item{Sp}{the diagnostic specificity of the test at the individual sample level as entered by the user.}
   \item{PlSe}{the diagnostic sensitivity of the test when samples are pooled.}
   \item{PlSp}{the diagnostic specificity of the test when samples are pooled.}
  \item{HAPneg}{the apparent prevalence in a disease negative herd.}
  \item{HSe}{the estimated group (herd) level sensitivity.}
  \item{HSp}{the estimated group (herd) level specificity.}
}

\references{
Dohoo I, Martin W, Stryhn H (2014). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 126 - 128
.

Christensen J, Gardner IA (2000). Herd-level interpretation of test results for epidemiologic studies of animal diseases. Preventive Veterinary Medicine 45: 83 - 106.  
}

\examples{
## EXAMPLE 1:
## We want to test dairy herds for Johne's disease using faecal culture 
## which has a diagnostic sensitivity and specificity of 0.647 and 
## 0.981, respectively. Suppose faecal samples from five cows are pooled 
## and we collect six pooled samples per herd. What is the herd level 
## sensitivity and specificity using this approach assuming the true prevalence
## of disease is 0.12 and there are no dilution effects arising from pooling. 

epi.pooled(Se = 0.647, Sp = 0.981, P = 0.12, m = 5 , r = 6, dilution = FALSE)

## Herd level sensitivity is 0.927, herd level specificity is 0.562. 
## Sensitivity at the herd level is increased using the pooled sampling 
## approach. Herd level specificity is decreased.

## Re-calculate herd level sensitivity and specificity assuming there is a 
## dilution effect arising from pooling:

epi.pooled(Se = 0.647, Sp = 0.981, P = 0.12, m = 5 , r = 6, dilution = TRUE)

## Herd level sensitivity is 0.943, herd level specificity is 0.562. 
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
