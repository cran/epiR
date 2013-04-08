\name{epi.pooled}

\alias{epi.pooled}

\title{
Estimate herd test characteristics when pooled sampling is used
}

\description{
We may wish to designate a group of individuals (e.g. a herd) as being either diseased or non-diseased on the basis of pooled samples. This function estimates sensitivity and specificity of this testing regime at the group (or herd) level. 
}

\usage{
epi.pooled(se, sp, P, m, r)
}

\arguments{
  \item{se}{a vector of length one defining the sensitivity of the individual test used.}
  \item{sp}{a vector of length one defining the specificity of the individual test used.}
  \item{P}{scalar, defining the estimated true prevalence.}
  \item{m}{scalar, defining the number of individual samples to make up a pooled sample.}
  \item{r}{scalar, defining the number of pooled samples per group (or herd).}
}

\value{
A list containing the following: 
  \item{HAPneg}{the apparent prevalence in a disease negative herd.}
  \item{HSe}{the estimated group (herd) level sensitivity.}
  \item{HSp}{the estimated group (herd) level specificity.}
}

\references{
Dohoo I, Martin W, Stryhn H (2003). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 115 - 117
.

Christensen J, Gardner IA (2000). Herd-level interpretation of test results for epidemiologic studies of animal diseases. Preventive Veterinary Medicine 45: 83 - 106.  
}

\examples{
## We want to test dairy herds for Johne's disease using faecal culture 
## which has a sensitivity and specificity of 0.647 and 0.981, respectively. 
## Suppose we pool faecal samples from five cows together and use six pooled 
## samples per herd. What is the herd level sensitivity and specificity 
## based on this approach (assuming homogenous mixing)? 

epi.pooled(se = 0.647, sp = 0.981, P = 0.12, m = 5 , r = 6)

## Herd level sensitivity is 0.927, herd level specificity is 0.562. 
## Sensitivity at the herd level is increased using the pooled sampling 
## approach; herd level specificity is decreased. 
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
