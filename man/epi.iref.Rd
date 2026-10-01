\name{epi.iref}

\alias{epi.iref}

\title{Calculate diagnostic sensitivity and specificity using an imperfect reference test}

\description{
Calculate diagnostic sensitivity and specificity using an imperfect reference test.
}

\usage{
epi.iref(x, se.rs, sp.rs, method = "staquet", ci.method = "wilson", 
   conf.level = 0.95, warn = TRUE)
}

\arguments{
  \item{x}{a vector of length four, an object of class \code{table} or an object of class \code{grouped_df} from package \code{dplyr} containing the individual cell frequencies (see below).}  
  
  \item{se.rs}{scalar, the known diagnostic sensitivity of the reference test. Must be a single number between 0 and 1.} 
  
  \item{sp.rs}{scalar, the known diagnostic specificity of the reference test. Must be a single number between 0 and 1.}   
  
  \item{method}{character string indicating the method to use. Options are \code{brenner}, \code{gart_buck}, or \code{staquet}.}  
  
  \item{ci.method}{character string indicating the confidence interval calculation method to use. For \code{method == "gart_buck"} and \code{method = "staquet"} options are \code{wilson} and\code{delta}. For \code{method == "brenner"} options are \code{wilson} and \code{wilson.eff}.}  
  
  \item{conf.level}{magnitude of the confidence interval of the sensitivity and specificity estimates. Must be a single number between 0 and 1.}
  
  \item{warn}{logical. If \code{TRUE} warnings are issued in the event of input data irregularities.}
}

\details{The required 2 by 2 table format for argument \code{x} for \code{epi.iref} is shown below. Columns list counts for the reference test: counts of test positive study units in column 1 and test negative study units in column 2. Rows list counts for the index test (i.e., the test under investigation): counts of test positive study units in row 1 and counts of test negative study units in row 2. The labels \code{a}, \code{b}, \code{c} and \code{d} in the table below correspond to the order of study unit counts when argument \code{dat} is a vector.

\tabular{llll}{
------------ \tab ------------ \tab ------------ \tab ------------ \cr 
 		         \tab Reference +  \tab Reference -  \tab Total        \cr
------------ \tab ------------ \tab ------------ \tab ------------ \cr                  
Test +  	   \tab a		         \tab b		         \tab a + b        \cr
Test -       \tab c		         \tab d		         \tab c + d        \cr
------------ \tab ------------ \tab ------------ \tab ------------ \cr 
Total        \tab a + c        \tab b + d        \tab a + b + c + d\cr
------------ \tab ------------ \tab ------------ \tab ------------ \cr 
}
}

\value{
A list containing the following: 
  \item{uncorrected}{diagnostic sensitivity and specificity for the index test (and their confidence intervals), computed using the values provided in table \code{x}.}
  \item{uncorrected}{diagnostic sensitivity and specificity for the index test (and their confidence intervals) computed using the reference test as a quasi gold standard.}
  \item{prevalence}{apparent prevalence and true prevalence (and their confidence intervals) computed using the reference test as the quasi gold standard.}  
}

\references{
Brenner H (1996). Correcting for exposure misclassification using an alloyed gold standard. Epidemiology 7: 406 - 410. 

Gart J, Buck A (1966). Comparison of a screening test and a reference test in epidemiologic studies. II. A probabilistic model for the comparison of diagnostic tests. American Journal of Epidemiology 83: 593 - 602. DOI: 10.1093/oxfordjournals.aje.a120610.

Habibzadeh F (2023). On determining the sensitivity and specificity of a new diagnostic test through comparing its results against a non-gold-standard test. Biochemia Medica: Casopis Hrvatskoga Drustva Medicinskih Biokemicara 33, 010101. DOI: 10.11613/BM.2023.010101.

Staquet M, Rozencweig M, Lee Y, Muggia F (1981). Methodology for the assessment of new dichotomous diagnostic tests. Journal of Chronic Diseases 34, 599 - 610. DOI: 10.1016/0021-9681(81)90059-x.

Chikere C, Wilson K, Allen A, Vale L (2021). Comparative diagnostic accuracy studies with an imperfect reference standard --- a comparison of correction methods. BMC Medical Research Methodology 21: 67. DOI: 10.1186/s12874-021-01255-4.
}

\author{
Mark Stevenson (Melbourne Veterinary School, Faculty of Science, The University of Melbourne, Australia).
}

\note{
The diagnostic accuracy measures (sensitivity and specificity) of a new test can be estimated by comparing test results to a gold standard.

Brenner (1996), Gart and Buck (1996), and Staquet et al. (1981) have published methods to estimate the sensitivity and specificity of a binary response index test when the sensitivity and specificity of the imperfect reference standard are known and the index test and reference standard are conditionally independent.

Based on a methodological review by Chikere et al. (2021), under the assumption of conditional independence, the Staquet et al. correction method outperforms the Brenner correction method irrespective of the prevalence of disease and whether the performance of the reference standard is better or worse than the index test. When the prevalence of the disease is high (> 0.9) or low (< 0.1), the Staquet et al. correction method can produce illogical results (i.e., results outside [0,1]). When the two tests are dependent both methods fail to estimate the sensitivity and specificity of the index test particularly when the covariance terms between the index test and the reference standard is not close to zero.

This function will be of most use in the situation when investigators have compared the results of an index test to an imperfect reference standard and calculated sensitivity and specificity, labelling them as 'relative' sensitivity and specificity estimates. On the condition that: (1) the diagnostic sensitivity and specificity of the imperfect reference test is known (i.e., the test has been validated against a gold standard); and (2) the prevalence of disease is within the range of approximately 0.1 to 0.9 (see above) this function is useful to (quickly) show how how far relative diagnostic sensitivity and specificity estimates are from 'actual' sensitivity and specificity estimates. Bayesian latent class modelling is the recommended approach to determine the diagnostic sensitivity and specificity of a test in the absence of a gold standard.
}

\examples{
## EXAMPLE 1 (from Habibzadeh 2023):
## The results of a new (index) diagnostic test were compared with an existing 
## reference. Of 150 individuals that were positive to the reference test, 107
## were positive to the index test. Of 450 individuals that were negative 
## to the reference test, 104 were positive to the index test. The known 
## sensitivity of the reference test is 0.850. The known specificity of the 
## reference test is 0.900.

## What is the diagnostic sensitivity and specificity of the index test?

x <- c(107,104,43,346)
epi.iref(x = x, se.rs = 0.850, sp.rs = 0.900, method = "staquet", 
   ci.method = "wilson", conf.level = 0.95, warn = TRUE)

## Without correction (i.e., assuming the reference test) has perfect
## diagnostic sensitivity and specificity the sensitivity and specificity of 
## the index test is 0.71 (95\% CI 0.63 to 0.78) and 0.77 (95\% CI 0.73 to 
## 0.81), respectively.

## With correction (i.e., accounting for imperfect performance of the reference 
## test the sensitivity and specificity of the index test is 0.95 (95\% CI 0.91 
## to 0.98) and 0.80 (95\% CI 0.76 to 0.83), respectively. 
}

\keyword{univar}
