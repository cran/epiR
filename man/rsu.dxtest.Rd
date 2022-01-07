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
 
\item{independent}{a data frame listing sensitivity \code{se} and specificity \code{sp} assuming the tests are independent.}
\item{dependent}{a data frame listing sensitivity \code{se} and specificity \code{sp} calculated using the values of \code{covar}, as entered by the user.}

If \code{covar = c(0,0)} data frames \code{independent} and \code{dependent} will be the same. 

}

\references{
Dohoo I, Martin S, Stryhn H (2009). Veterinary Epidemiologic Research. AVC Inc Charlottetown, Prince Edward Island, Canada.

Gardner I, Stryhn H, Lind P, Collins M (2000). Conditional dependence between tests affects the diagnosis and surveillance of animal diseases. Preventive Veterinary Medicine 45: 107 - 122.

Martin S, Meek A, Willeberg P (1987). Veterinary Epidemiology Principles and Methods. Iowa State University Press Ames.
}

\note{
With \code{interpretation = "series"} a subject is declared test positive if both of the tests performed return a positive result. With \code{interpretation = "parallel"} a subject is declared test positive if one of the tests performed return a positive result. Intepreting test results in series increases diagnostic specificity. Interpreting test results in parallel increases diagnostic sensitivity.

How do I work out appropriate values for \code{covar}? Assume you have two diagnostic tests --- an indirect fluorescent antibody test (IFAT) and a polymerase chain reaction (PCR). The diagnostic sensitivity and specificity of the IFAT is 0.784 and 0.951, respectively. The diagnostic sensitivity and specificity of the PCR is 0.926 and 0.979, respectively. These tests are used on a group of individuals known to be disease positive and a group of individuals known to be disease negative. Results for the disease positive group are as follows:  

\tabular{lllll}{
           \tab IFAT      \tab            \tab            \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr
PCR        \tab Pos      	\tab Neg        \tab Total      \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr                   
Pos     	 \tab 134       \tab 29 	  	  \tab 163        \cr
Neg        \tab 4		      \tab 9          \tab 13         \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr 
Total      \tab 138     	\tab 38         \tab 176        \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr
}

Results for the disease negative group are as follows:  

\tabular{lllll}{
           \tab IFAT      \tab            \tab            \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr
PCR        \tab Pos      	\tab Neg        \tab Total      \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr                   
Pos     	 \tab 0         \tab 12 	  	  \tab 12         \cr
Neg        \tab 28	      \tab 534        \tab 562        \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr 
Total      \tab 28      	\tab 546        \tab 574        \cr
-----------\tab ----------\tab ---------- \tab ---------- \cr
}

The observed proportion of disease positive individuals with a positive test result to both tests as \code{p111}. For this example \code{p111 = 134 / 176 = 0.761}.

The observed proportion of disease negative individuals with a negative test result to both tests as \code{p000}. For this example \code{p000 = 534 / 574 = 0.930}.

Covariance for the disease positive group: \code{covar[1] = p111 - se[1] * se[2] = 0.761 - 0.784 * 0.926 = 0.035}.
 
Covariance for the disease negative group: \code{covar[2] = p000 - sp[1] * sp[2] = 0.930 - 0.951 * 0.979 = -0.001}.

The covariance for the disease positive group is small. The covariance for the disease negative group is negligible.
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
## of 0.855 and a diagnostic specificity of 0.970.


## EXAMPLE 2 (from Dohoo, Martin and Stryhn p 113):
## An IFAT and PCR are to be used to diagnose infectious salmon anaemia. 
## The diagnostic sensitivity and specificity of the IFAT is 0.784 and 0.951, 
## respectively. The diagnostic sensitivity and specificity of the PCR is 
## 0.926 and 0.979, respectively. It is known that the two tests are dependent, 
## with details of the covariance calculated above. What is the expected 
## sensitivity and specificity if the tests are to be interpreted in parallel?   

rsu.dxtest(se = c(0.784,0.926), sp = c(0.951,0.979), 
   interpretation = "parallel", covar = c(0.035,-0.001))

## Interpreting test results in parallel and accounting for the lack of 
## test indepdendence returns a diagnostic sensitivity of 0.949 and diagnostic 
## specificity of 0.930. 

}
\keyword{methods}

