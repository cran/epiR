\name{epi.nomogram}

\alias{epi.nomogram}

\title{Post-test probability of an outcome given diagnostic test sensitivity and specificity}

\description{
Compute the post-test probability of an outcome given diagnostic test sensitivity and specificity.
}

\usage{
epi.nomogram(pretest.ppos, se, sp, lratio.pos = NA, lratio.neg = NA, 
   tconf.int = 0.95, method = "exact", verbose = FALSE, cred.int = 0.95, 
   nsim = 999)
}

\arguments{
  \item{pretest.ppos}{the pre-test probability of being outcome positive (0 - 1).}
  \item{se}{diagnostic test sensitivity (0 - 1). Either a vector of length one listing the point estimate of the diagnostic test sensitivity or a vector of length three listing (in order) the point estimate, the lower bound of the confidence interval and the upper bound of the confidence interval of the diagnostic test sensitivity.}
  \item{sp}{diagnostic test specificity (0 - 1). Either a vector of length one listing the point estimate of the diagnostic test specificity or a vector of length three listing (in order) the point estimate, the lower bound of the confidence interval and the upper bound of the confidence interval of the diagnostic test specificity.}
  \item{lratio.pos}{positive likelihood ratio. A vector of length one. Ignored if \code{se} and \code{sp} are not null.}
  \item{lratio.neg}{negative likelihood ratio. A vector of length one. Ignored if \code{se} and \code{sp} are not null.}
  \item{tconf.int}{scalar, the magnitude of the confidence interval for the listed diagnostic test performance measures (\code{se} and \code{sp}).}
  \item{method}{a character string indicating the method used to calculate the confidence intervals for \code{se} and \code{sp}. Options are \code{method = "exact"}, \code{method = "wilson"}, \code{method = "agresti"}, \code{method = "clopper-pearson"} and \code{method = "jeffreys"}.}
  \item{verbose}{logical, indicating whether detailed or summary results are to be returned.}
  \item{cred.int}{magnitude of the returned credible interval for post-test probability. Must be a single number between 0 and 1.}
  \item{nsim}{scalar, defining the number of simulations to be used for the post-test probability credible interval calculations.}
}

\details{
To return a credible interval for the post-test probability this function uses the objective Bayesian method described by Mossman and Berger (2001). In brief, the confidence limits entered by the user for \code{se} and \code{sp} and specification of the method used to estimate those confidence limits (i.e., \code{method = "exact"}, \code{method = "wilson"}, \code{method = "agresti"}, \code{method = "clopper-pearson"} or \code{method = "jeffreys"}) are used to return the number of study subjects in the study from which the diagnostic test accuracy results were obtained (e.g., the number of individuals outcome and test positive individuals and the number of outcome positive individuals outcome positive for diagnostic sensitivity). Random draws from binomial distributions using shape parameters based on the number of study subjects are then used to generate a distribution of outcome positive post-test probabilities and a distribution of outcome negative post-test probabilities. The setting for \code{credible.int} defines the magnitude of the reported credible interval for the post-test probabilities.  

}

\value{
A list containing the following:
  \item{pretest.ppos}{the pre-test probability of being outcome positive, as entered by the user.}
  \item{pretest.pneg}{the pre-test probability of being outcome negative, \code{1 - pretest.ppos}.}
  \item{postest.ppos}{the post-test probability of being outcome positive and the lower and upper bounds of the credible interval of the post-test probability of being outcome positive if confidence limits for \code{se} and \code{sp} have been provided.}
  \item{postest.pneg}{the post-test probability of being outcome negative and the lower and upper bounds of the credible interval of the post-test probability of being outcome negative if confidence limits for \code{se} and \code{sp} have been provided.}
    \item{lratio.pos}{the positive likelihood ratio.}
    \item{lratio.neg}{the negative likelihood ratio.}
}

\references{
Caraguel C, Vanderstichel R (2013). The two-step Fagan's nomogram: ad hoc interpretation of a diagnostic test result without calculation. Evidence Based Medicine 18: 125 - 128.

Haskins R, Osmotherly P, Tuyl F, Rivett D, (2014). Uncertainty in clinical prediction rules: The value of credible intervals. Journal of Orthopaedic and Sports Physical Therapy 44:  85 - 91. DOI: 10.2519/jospt.2014.4877.

Hunink M, Glasziou P (2001). Decision Making in Health and Medicine - Integrating Evidence and Values. Cambridge University Press, pp. 128 - 156.

Mossman D, Berger J (2001) Intervals for posttest probabilities: A comparison of 5 methods. Medical Decision Making 21:  498 - 507.
}

\examples{
## EXAMPLE 1:
## You are presented with a dog with lethargy, exercise intolerance, 
## weight gain and bilaterally symmetric truncal alopecia. You are 
## suspicious of hypothyroidism and take a blood sample to measure 
## basal serum thyroxine (T4).

## You believe that around 5\% of dogs presented to your clinic with 
## a signalment of general debility have hypothyroidism. The serum T4 
## has a sensitivity of 0.89 and specificity of 0.85 for diagnosing 
## hypothyroidism in the dog. The laboratory reports a serum T4 
## concentration of 22.0 nmol/L (reference range 19.0 to 58.0 nmol/L). 
## What is the post-test probability that this dog is hypothyroid?

epi.nomogram(pretest.ppos = 0.05, se = 0.890, sp = 0.850, 
   lratio.pos = NA, lratio.neg = NA, tconf.int = NA, method = "exact", 
   verbose = FALSE, cred.int = 0.95, nsim = 999)

## Given the reference range of 19.0 to 58.0 nmol/L and a test result of 
## 22.0 nmol/L we declare this a negative test result.

## If the pre-test probability of being hypothyroid positive is 0.05 and 
## the test is negative, the post-test probability of being hypothyroid 
## positive is 0.0068.


## EXAMPLE 2:
## A test with diagnostic sensitivity 0.980 (95\% CI 0.894 to 0.999) and
## diagnostic specificity 0.980 (95\% CI 0.893 to 0.999) has been developed 
## to detect feline hepatic lipidosis (FHL). After conducting a  thorough 
## clinical examination you are 75\% sure that your patient is FHL positive. 
## Using the assigned cutpoint for this test, your patient returns a positive 
## result. What is the post-test probability that this cat is FHL positive?

epi.nomogram(pretest.ppos = 0.75, se = c(0.980,0.894,0.999), 
   sp = c(0.980,0.893,0.999), lratio.pos = NA, lratio.neg = NA, 
   tconf.int = 0.95, method = "exact", verbose = FALSE, cred.int = 0.95, 
   nsim = 999)
   
## If the pre-test probability of being FHL positive is 0.75 and the test 
## is positive, the post-test probability of being FHL positive is 
## 0.99 (95% CrI 0.97 to 1.00). 

## If the pre-test probability of being FHL positive is 0.75 and the test 
## is negative, the post-test probability of being FHL positive is 
## 0.06 (95% CrI 0.0060 to 0.2000). If the cat returns a negative test result 
## there's up to a 20\% chance that it could still be FHL positive. False 
## negatives are a potential problem with this test. 

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
