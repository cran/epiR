\name{rsu.dxtest}

\alias{rsu.dxtest}

\title{
Sensitivity and specificity of diagnostic tests interpreted in series or parallel
}

\description{
Calculates the sensitivity and specificity of two diagnostic tests interpreted in series or parallel.
}

\usage{
rsu.dxtest(se, sp, interpretation = "series", covar = c(0,0))
}

\arguments{
\item{se}{a vector of length two defining the diagnostic sensitivity of the two tests.}
\item{sp}{a vector of length two defining the diagnostic specificity of the two tests.}
\item{interpretation}{a character string indicating how the test results should be interpreted. Options are \code{series} or \code{parallel}.}
\item{covar}{a vector of length two defining the covariance between test results for disease positive and disease negative groups. The first element of the vector is the covariance between test results for disease positive subjects. The second element of the vector is the covariance between test results for disease negative subjects. Use \code{covar = c(0,0)} (the default) if the tests are known to be independent.}
}

\value{
A list comprised of two elements:
 
\item{se}{scalar, the estimated diagnostic sensitivity.}
\item{sp}{scalar, the estimated diagnostic specificity.}
}

\references{
Dohoo I, Martin S, Stryhn H (2009). Veterinary Epidemiologic Research. AVC Inc Charlottetown, Prince Edward Island, Canada.

Martin S, Meek A, Willeberg P (1987). Veterinary Epidemiology Principles and Methods. Iowa State University Press Ames.
}

\note{
With \code{interpretation = "series"} a subject is declared test positive if both of the tests performed return a positive result. With \code{interpretation = "parallel"} a subject is declared test positive if one of the tests performed return a positive result. 

Intepreting test results in series increases diagnostic specificity. Interpreting test results in parallel increases diagnostic sensitivity. 
}


\examples{
## EXAMPLE 1:
## You would like to confirm the absence of disease in a study area. You 
## intend to use two tests: the first has a sensitivity and specificity of 
## 0.90 and 0.80, respectively. The second has a sensitivity and specificity 
## of 0.95 and 0.85, respectively. You need to make sure that an individual
## that returns a positive test really has disease, so the tests will be 
## interpreted in series (to improve specificity). 

## What is the diagnostic sensitivity and specificity of this testing 
## regime? 

rsu.dxtest(se = c(0.90,0.95), sp = c(0.80,0.85), 
   interpretation = "series", covar = c(0,0))
   
## Interpretation of these tests in series returns a diagnostic sensitivity
## of 0.855 and a diagnostic specificity of 0.97.

}
\keyword{methods}

