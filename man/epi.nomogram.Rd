\name{epi.nomogram}

\alias{epi.nomogram}

\title{Post-test probability of disease given sensitivity and specificity of a test}

\description{
Compute the post-test probability of disease given sensitivity and specificity of a test.
}

\usage{
epi.nomogram(se, sp, lr, pre.pos, verbose = FALSE)
}

\arguments{
  \item{se}{test sensitivity (0 - 1).}
  \item{sp}{test specificity (0 - 1).}
  \item{lr}{a vector of length 2 listing the positive and negative likelihood ratio (respectively) of the test. Ignored if \code{se} and \code{sp} are not null.}
  \item{pre.pos}{the pre-test probability of the outcome.}
  \item{verbose}{logical, indicating whether detailed or summary results are to be returned.}
}

\value{
A list containing the following:
  \item{lr}{the likelihood ratio of a positive and negative test.}
  \item{prob}{the post-test probability of the outcome given a positive and negative test.}
}

\references{
Caraguel C, Vanderstichel R (2013). The two-step Fagan's nomogram: ad hoc interpretation of a diagnostic test result without calculation. Evidence Based Medicine 18: 125 - 128.

Hunink M, Glasziou P (2001). Decision Making in Health and Medicine - Integrating Evidence and Values. Cambridge University Press, pp. 128 - 156.
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

epi.nomogram(se = 0.89, sp = 0.85, lr = NA, pre.pos = 0.05, verbose = FALSE)

## Given a positive test result, the post-test probability of being 
## disease positive is 0.24.

## Given a negative test result, the post-test probability of being 
## disease negative is 0.0068.

## EXAMPLE 2:
## A dog is presented to you with severe pruritis. You suspect sarcoptic 
## mange and decide to take a skin scraping (LR+ 9000; LR- 0.1). The scrape 
## returns a negative result (no mites are seen). What is the post-test 
## probability that your patient has sarcoptic mange? You recall that you 
## diagnose around 3 cases of sarcoptic mange per year in a clinic that 
## sees approximately 2 -- 3 dogs per week presented with pruritic skin disease.

pre.pos <- 3 / (3 * 52)
epi.nomogram(se = NA, sp = NA, lr = c(9000, 0.1), pre.pos = pre.pos, 
   verbose = FALSE)

## If the skin scraping is negative the post-test probability that this dog 
## has sarcoptic mange is 0.002.

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
