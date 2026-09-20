\name{epi.mbf}

\alias{epi.mbf}

\title{Minimum Bayes factor using the Goodman or Sellke et al. method.}

\description{
A function to calculate the minimum Bayes factor using the Goodman or Sellke et al. method.
}

\usage{
epi.mbf(H0, class = c("prior", "posterior"), p, method = c("goodman", "sellke"), 
   two.sided = TRUE)
}

\arguments{
  \item{H0}{scalar, prior or posterior estimate of the null hypothesis being true.}
  \item{class}{character, \code{prior} or \code{posterior}.}
  \item{p}{scalar, the p-value calculated for the comparison of interest.}
  \item{method}{character, which minimum Bayes factor calculation method to use. Options \code{goodman} or \code{sellke}. See details, below.}
  \item{two.sided}{logical, if \code{TRUE} the analysis assumes a two-sided test has been used.}
}

\value{
A list containing the following:
  \item{method}{the minimum Bayes factor method (\code{goodman} or \code{sellke}) as entered by the user.}
  
  \item{H0.prior}{the prior estimate that the null hypothesis is true. If \code{class} is set to \code{prior} this is the value of \code{H0} entered by the user. If \code{class} is set to \code{posterior} this is the minimum Bayes factor estimate of \code{H0}.}
  
  \item{H1.prior}{the prior estimate that the alternative hypothesis is true. If \code{class} is set to \code{prior} this is 1 minus the value of \code{H0} entered by the user. If \code{class} is \code{posterior} this is 1 - the minimum Bayes factor estimate of \code{H0}.}
  
  \item{H0.post}{the posterior estimate that the null hypothesis is true. If \code{class} is set to \code{prior} this is minimum Bayes factor estimate of \code{H0}. If \code{class} is set to \code{posterior} this is the value of \code{H0} entered by the user.}
  
  \item{H1.post}{the posterior estimate that the alternative  hypothesis is true. If \code{class} is set to \code{prior} this is 1 minus the minimum Bayes factor estimate of \code{H0}. If \code{class} is set to \code{posterior} this is 1 minus the value of \code{H0} entered by the user.}
  
  \item{interp}{a one sentence statement of the results.}
}

\details{
A Bayes factor is a number that tells you how much a piece of evidence should shift your belief between two competing explanations such as: 'this treatment works' versus 'ths treatment doesn't work'. 

The 'minimum' Bayes factor is the smallest possible shift in belief that a given piece of evidence (like a p-value from a study) could justify, even under the most generous assumptions favouring the surprising result. It is used to show that even seemingly strong evidence (such as a low p-value) often provides weaker support for a hypothesis than people intuitively assume, because it calculates the best-case scenario for that evidence and often still finds it underwhelming.

The Goodman (1999) method compares the null hypothesis to a single, specific alternative: The value of the parameter that best fits the data (the maximum likelihood estimate).

The Sellke, Bayarri and Berger (2001) method considers an entire class of reasonable prior distributions on the alternative hypothesis (essentially, any distribution that doesn't put special weight right at the null) and finds the minimum Bayes factor across that whole class.

The Sellke, Bayarri and Berger method is consistently more conservative (larger minimum Bayes factor estimates, i.e., weaker evidence against the null) than the Goodman method. This makes sense: the Goodman method is comparing against the single best-fitting alternative, which is the most favourable comparison possible. The Sellke-Bayarri-Berger method asks for the worst case across many reasonable alternatives, which can't be as extreme as cherry-picking the one alternative that happens to match the data perfectly.
}

\references{
Goodman S (1999) Toward evidence-based medical statistics. 1: The P value fallacy. Annals of Internal Medicine 130: 995 - 1004. DOI: 10.7326/0003-4819-130-12-199906150-00008.

Sellke T, Bayarri M, Berger J (2001) Calibration of p values for testing precise null hypotheses. American Statistician 55: 62 - 71.

}

\examples{
## EXAMPLE 1:
## A pharmaceutical company tests a new blood pressure medication against a 
## placebo in a clinical trial with 200 patients. At the start of the trial
## the researchers are moderately optimistic. Previous research on similar 
## compounds suggests the drug probably works, but there's real uncertainty. 
## They set their prior probability that H1 (the drug is effective) is true 
## at 0.60.

## The trial finds a statistically significant reduction in blood pressure 
## in the treatment group, with a p-value of 0.04.

epi.mbf(H0 = 0.40, class = "prior", p = 0.04, method = "goodman", 
   two.sided = TRUE)

## The calculated minimum Bayes factor is 0.121. With a prior estimate of the 
## probability that H1 is true of 0.6 the data provide at most a posterior 
## probability of 0.925 in favour of H1.

## Repeat these calculations using the Sellke, Bayarri and Berger method:

epi.mbf(H0 = 0.40, class = "prior", p = 0.04, method = "sellke", 
   two.sided = TRUE)
   
## The calculated minimum Bayes factor is 0.350. With a prior estimate of the 
## probability that H1 is true of 0.6 the data provide at most a posterior 
## probability of 0.811 in favour of H1.

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
