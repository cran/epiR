\name{epi.nomogram}

\alias{epi.nomogram}

\title{Post-test probability of disease given sensitivity and specificity of a test}

\description{
Computes the post-test probability of disease given sensitivity and specificity of a test.
}

\usage{
epi.nomogram(se, sp, pre.pos, verbose = FALSE)
}

\arguments{
  \item{se}{test sensitivity (0 - 1).}
  \item{sp}{test specificity (0 - 1).}
  \item{pre.pos}{the pre-test probability of disease in the patient.}
  \item{verbose}{logical, indicating whether detailed or summary results are to be returned.}
}

\details{
}

\value{
A list containing the following:
  \item{likelihood.ratio}{the likelihood ratio of a positive and negative test.}
  \item{prob}{the post-test probability of disease given a positive and negative test.}
}

\references{
Hunink M, Glasziou P (2001). Decision Making in Health and Medicine - Integrating Evidence and Values. Cambridge University Press, pp. 128 - 156.
}

\author{
}

\note{
}

\seealso{
}

\examples{
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

epi.nomogram(se = 0.89, sp = 0.85, pre.pos = 0.05, verbose = FALSE)

## The post-test probability that this dog is hypothyroid is 24\%.

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
