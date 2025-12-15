\name{rsu.dxtest}

\alias{rsu.dxtest}

\title{
Sensitivity and specificity of diagnostic tests interpreted in series or parallel
}

\description{
Calculates the sensitivity and specificity of two or three diagnostic tests interpreted in series or parallel.
}

\usage{
rsu.dxtest(se, sp, covar.pos, covar.neg, tconf.int, method = "exact", 
   interpretation = "series", conf.int = 0.95, nsim = 999)
}

\arguments{
\item{se}{a matrix of diagnostic test sensitivities. For tests without confidence intervals for sensitivity this will be a 2 or 3-row, 1-column matrix with rows 1 to 3 representing diagnostic tests and column 1 listing the point estimate of diagnostic sensitivity. For tests with confidence intervals for sensitivity this will be a 2 or 3-row, 3-column matrix with rows 1 to 3 representing diagnostic tests and columns 1, 2 and listing the point estimate, the lower bound of the sensitivity confidence interval and the upper bound of the sensitivity confidence interval, respectively.}

\item{sp}{a matrix of diagnostic test specificities. For tests without confidence intervals for specificity this will be a 2 or 3-row, 1-column matrix with rows 1 to 3 representing diagnostic tests and column 1 listing the point estimate of diagnostic specificity. For tests with confidence intervals for specificity this will be a 2 or 3-row, 3-column matrix with rows 1 to 3 representing diagnostic tests and columns 1, 2 and listing the point estimate, the lower bound of the specificity confidence interval and the upper bound of the specificity confidence interval, respectively.}

\item{covar.pos}{a vector of length either two or four defining the covariance between test results for the disease positive group. For two diagnostic tests \code{covar.pos} is a vector of length one. For three diagnostic tests \code{covar.pos} is a vector of length four listing (in order) the covariance between test 1 and 2, test 1 and 3, test 2 and 3 and tests 1, 2 and 3. See the examples, below for details.}

\item{covar.neg}{a vector of length either two or four defining the covariance between test results for the disease negative group. For two diagnostic tests \code{covar.neg} is a vector of length one. For three diagnostic tests \code{covar.neg} is a vector of length four listing (in order) the covariance between test 1 and 2, test 1 and 3, test 2 and 3 and tests 1, 2 and 3. See the examples, below for details.}

\item{tconf.int}{scalar, the magnitude of the confidence interval for the listed diagnostic test performance measures (\code{se} and \code{sp}).}

\item{method}{a character string indicating the method used to calculate the confidence intervals for \code{se} and \code{sp}. Options are \code{method = "exact"}, \code{method = "wilson"}, \code{method = "agresti"}, \code{method = "clopper-pearson"} and \code{method = "jeffreys"}.}

\item{interpretation}{a character string indicating how the test results should be interpreted. Options are \code{series} or \code{parallel}.}

\item{conf.int}{magnitude of the returned confidence interval for sensitivity and specificity if tests are interpreted in series or parallel. Must be a single number between 0 and 1.}

\item{nsim}{scalar, defining the number of simulations to be used to estimate confidence interval for sensitivity and specificity if tests are interpreted in series or parallel.}
}

\value{
A list comprised of four elements:
 
\item{independent}{a data frame listing sensitivity \code{se} and specificity \code{sp} assuming the tests are independent.}
\item{dependent}{a data frame listing sensitivity \code{se} and specificity \code{sp} calculated using the values of \code{covar.pos} and \code{covar.neg}, as entered by the user.}
\item{covar.pos}{a data frame listing the covariance between test results for the disease positive group, as entered by the user.}
\item{covar.neg}{a data frame listing the covariance between test results for the disease negative group, as entered by the user.}

If \code{covar.pos = 0} and \code{covar.neg = 0} (for the two test situation) data frames \code{independent} and \code{dependent} will be identical. 

}

\references{
Dohoo I, Martin S, Stryhn H (2009). Veterinary Epidemiologic Research. AVC Inc Charlottetown, Prince Edward Island, Canada.

Gardner I, Stryhn H, Lind P, Collins M (2000). Conditional dependence between tests affects the diagnosis and surveillance of animal diseases. Preventive Veterinary Medicine 45: 107 - 122.

Jones G, Johnson W, Hanson T, Christensen R (2010). Identifiability of models for multiple diagnostic testing in the absence of a gold standard. Biometrics 66:  855 - 863. DOI: 10.1111/j.1541-0420.2009.01330.x.

Martin S, Meek A, Willeberg P (1987). Veterinary Epidemiology Principles and Methods. Iowa State University Press Ames.

McClure C, Hammell K, Stryhn H, Dohoo I, Hawkins L (2005). Application of surveillance data in evaluation of diagnostic tests for infectious salmon anemia. Diseases of Aquatic Organisms 63: 119 - 127. DOI: 10.3354/dao063119.

Toft N, Akerstedt J, Tharaldsen J, Hopp P (2007). Evaluation of three serological tests for diagnosis of Maedi-Visna virus infection using latent class analysis. Veterinary Microbiology 120: 77 - 86.
}

\note{
For \code{interpretation = "series"} a subject is declared test positive if both of the tests performed return a positive result. For \code{interpretation = "parallel"} a subject is declared test positive if one of the tests performed return a positive result. Intepreting test results in series increases diagnostic specificity. Interpreting test results in parallel increases diagnostic sensitivity.

For the two test situation \code{covar.pos} is a vector of length one quantifying the covariance between between the two tests for the disease positive group. For the three test situation \code{covar.pos} is a vector of length four listing (in order) the the covariance between test 1 and 2, test 1 and 3, test 2 and 3 and tests 1, 2 and 3 for the disease positive group: 

\tabular{ll}{
Element           \tab Test pair         \cr
------------------\tab ------------------\cr
1                 \tab Test 1 - Test 2   \cr
2                 \tab Test 1 - Test 3   \cr
3                 \tab Test 2 - Test 3   \cr
4                 \tab Test 1, 2, and 3  \cr
------------------\tab ------------------\cr
}

For the two test situation \code{covar.neg} is a vector of length one quantifying the covariance between between the two tests for the disease negative group. For the three test situation \code{covar.neg} is a vector of length four listing (in order) the the covariance between test 1 and 2, test 1 and 3, test 2 and 3 and tests 1, 2 and 3 for the disease negative group: 

\tabular{ll}{
Element           \tab Test pair         \cr
------------------\tab ------------------\cr
1                 \tab Test 1 - Test 2   \cr
2                 \tab Test 1 - Test 3   \cr
3                 \tab Test 2 - Test 3   \cr
4                 \tab Test 1, 2, and 3  \cr
------------------\tab ------------------\cr
}

How do I work out appropriate values for \code{covar.pos} and \code{covar.neg}? Assume you have two diagnostic tests --- an indirect fluorescent antibody test (IFAT) and a polymerase chain reaction (PCR). The diagnostic sensitivity and specificity of the IFAT is 0.784 and 0.951, respectively. The diagnostic sensitivity and specificity of the PCR is 0.926 and 0.979, respectively. These tests are used on a group of individuals known to be disease positive and a group of individuals known to be disease negative. Results for the disease positive group are as follows:  

\tabular{llll}{
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

\tabular{llll}{
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

Covariance for the disease positive group: \code{covar.pos = p111 - se[1] * se[2] = 0.761 - 0.784 * 0.926 = 0.035}.
 
Covariance for the disease negative group: \code{covar.neg = p000 - sp[1] * sp[2] = 0.930 - 0.951 * 0.979 = -0.001}.

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

se <- matrix(c(0.90,0.95), ncol = 1, byrow = TRUE); se
sp <- matrix(c(0.80,0.85), ncol = 1, byrow = TRUE); sp

rsu.dxtest(se = se, sp = sp, covar.pos = 0 , covar.neg = 0, 
   tconf.int = 0.95, method = "exact", interpretation = "series", 
   conf.int = 0.95, nsim = 999)

## Interpretation of these tests in series returns a diagnostic sensitivity
## of 0.855 and a diagnostic specificity of 0.970.


## EXAMPLE 2 (from Dohoo, Martin and Stryhn p 113):
## An IFAT and PCR are to be used to diagnose infectious salmon anaemia (IFA). 
## Counts of IFA positive salmon that tested positive to the IFAT and PCR and 
## counts of IFA negative salmon that tested positive to the IFAT and PCR are
## shown above. Calculate the sensitivity and specificity of the two tests
## and a 95\% confidence interval for sensitivity and specificity of the 
## two tests using the exact method.

test <- rep(c("ifat","pcr"), each = 2)
perf <- rep(c("se","sp"), times = 2)
num <- c(138,546,163,562)
den <- c(176,574,176,574)
dat.df <- data.frame(test, perf, num, den)

tmp <- epi.conf(dat = as.matrix(dat.df[,3:4]), ctype = "prevalence", 
   method = "exact", N = 1000, design = 1, conf.level = 0.95)
dat.df <- cbind(dat.df, tmp); dat.df

## The diagnostic sensitivity and specificity of the IFAT is 0.784 
## (95\% CI 0.716 to 0.842) and 0.951 (95\% CI 0.930 to 0.967), respectively.

## The diagnostic sensitivity and specificity of the PCR is 
## 0.926 (95\% CI 0.877 to 0.960) and 0.979 (95\% CI 0.964 to 0.989), 
## respectively. It is known that the two tests are dependent, 
## with positive and negative covariances calculated above. What is the 
## expected sensitivity and specificity if the tests are to interpreted 
## in parallel?   

se <- matrix(c(0.784,0.716,0.842,0.926,0.877,0.960), ncol = 3, 
   byrow = TRUE); se
sp <- matrix(c(0.951,0.930,0.967,0.979,0.964,0.989), ncol = 3, 
   byrow = TRUE); sp

rsu.dxtest(se = se, sp = sp, covar.pos = 0.035, covar.neg = -0.001, 
   tconf.int = 0.95, method = "exact", interpretation = "parallel", 
   conf.int = 0.95, nsim = 999)

## Interpreting test results in parallel and accounting for the lack of 
## test indepdendence returns a diagnostic sensitivity of 0.949 (95\% CI 0.938 
## to 0.957) and diagnostic specificity of 0.929 (95\% CI 0.906 to 0.947). 


## EXAMPLE 3:
## Three diagnostic tests for Brucella suis in dogs are available: the Rose
## Bengal test (RBT), complement fixation (CFT) and an c-ELISA. The diagnostic
## sensitivities of the three tests are 0.869, 0.873 and 0.897, respectively.
## The diagnostic specificities of the three tests are 0.895, 0.905, and 0.915, 
## respectively.

## For disease positive groups the covariance between test results are 
## as follows:

## RBT - CFT: 0.071
## RBT - c-ELISA: 0.018
## CFT - c-ELISA: 0.015
## RBT - CFT - c-ELISA: -0.001

covp012 <- 0.071
covp013 <- 0.018
covp023 <- 0.015
covp123 <- -0.001
covp <- c(covp012,covp013,covp023,covp123)

## For disease negative groups the covariance between test results are 
## as follows:

## RBT - CFT: 0.076
## RBT - c-ELISA: 0.062
## CFT - c-ELISA: 0.062
## RBT - CFT - c-ELISA: 0.049

covn012 <- 0.076
covn013 <- 0.062
covn023 <- 0.062
covn123 <- 0.049
covn <- c(covn012,covn013,covn023,covn123)

## What is the expected diagnostic sensitivity and specificity if all three 
## tests are run on an individual and interpreted in parallel?

se <- matrix(c(0.869,0.873,0.897), ncol = 1, byrow = TRUE); se
sp <- matrix(c(0.895,0.905,0.915), ncol = 1, byrow = TRUE); sp

rsu.dxtest(se = se, sp = sp, covar.pos = covp, covar.neg = covn, 
   tconf.int = 0.95, method = "exact", interpretation = "parallel", 
   conf.int = 0.95, nsim = 999)

## Interpreting the test results in parallel and accounting for depdendence 
## between the three tests returns a diagnostic sensitivity of 
## 0.991 and a diagnostic specificity of 0.873. If we didn't account for test
## dependence our estimate of diagnostic sensitivity and specificity would be 
## 0.998 and 0.741, respectively.

## What is the expected sensitivity and specificity if all three 
## tests are run on an individual and interpreted in series?

rsu.dxtest(se = se, sp = sp, covar.pos = covp, covar.neg = covn, 
   tconf.int = 0.95, method = "exact", interpretation = "series", 
   conf.int = 0.95, nsim = 999)
   
## Interpreting the test results in series and accounting for depdendence 
## between the three tests returns a diagnostic sensitivity of 
## 0.774 and a diagnostic specificity of 1.000.
}


\keyword{methods}

