\name{epi.fpos}

\alias{epi.fpos}

\title{
Calculate the expected number of false positive and false negative test results
}

\description{
Compute the number of true positives, false positives, true negatives and false negatives given the number of individuals tested, the design prevalence and diagnostic test sensitivity and specificity.
}

\usage{
epi.fpos(n, pstar, se.u, sp.u, conf.level = 0.95)
}

\arguments{
  \item{n}{scalar, integer, the number of surveillance units tested.}
  \item{pstar}{scalar, the expected design prevalence.}
  \item{se.u}{scalar, the surveillance unit sensitivity.}
  \item{sp.u}{scalar, the surveillance unit specificity.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\value{
A list containing the following:
 
   \item{test.pos}{the expected median number of test positives and the uncertainty in the number of test positives, as specified by \code{conf.level}.}
   \item{true.pos}{the expected median number of true positives and the uncertainty in the number of true positives, as specified by \code{conf.level}.}
   \item{false.pos}{the expected median number of false positives units and the uncertainty in the number of false positives, as specified by \code{conf.level}.}   
   
   \item{test.neg}{the expected median number of test negatives and the uncertainty in the number of test negative, as specified by \code{conf.level}.}
   \item{true.neg}{the expected median number of true negatives and the uncertainty in the number of true negative, as specified by \code{conf.level}.}
   \item{false.neg}{the expected median number of false negatives and the uncertainty in the number of false negative, as specified by \code{conf.level}.}  
}

\examples{
## EXAMPLE 1:
## A bulk milk tank test has been developed to detect bovine tuberculosis 
## (bTB) in dairy herds. The diagnostic sensitivity and specificity of the 
## test is 0.714 and 0.981, respectively.

## If there are 9,000 dairy herds in the population of interest and of that 
## group 15 are thought to be bTB positive how many false positive test 
## results can be expected if all herds are tested on a single occasion?

epi.fpos(n = 9000, pstar = 15 / 9000, se.u = 0.714, sp.u = 0.981, 
   conf.level = 0.95)$false.pos

## If all 9,000 herds are tested on a single occasion we can expect 171 
## (95\% CI 146 to 197) false positive.

## How many false negatives are expected using this test regime?

epi.fpos(n = 9000, pstar = 15 / 9000, se.u = 0.714, sp.u = 0.981, 
   conf.level = 0.95)$false.neg

## If all 9,000 herds are tested on a single occasion we can expect 4 
## (95\% CI 1 to 9) false negatives.
}

\keyword{univar}
