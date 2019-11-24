\name{epi.smd}

\alias{epi.smd}

\title{Fixed-effects meta-analysis of continuous outcomes using the standardised mean difference method
}

\description{
Computes the standardised mean difference and confidence intervals of the standardised mean difference for continuous outcome data.
}

\usage{
epi.smd(mean.trt, sd.trt, n.trt, mean.ctrl, sd.ctrl, n.ctrl, 
   names, method = "cohens", conf.level = 0.95)
}

\arguments{
  \item{mean.trt}{a vector, defining the mean outcome in the treatment group.}
  \item{sd.trt}{a vector, defining the standard deviation of the outcome in the treatment group.}
  \item{n.trt}{a vector, defining the number of subjects in the treatment group.}
  \item{mean.ctrl}{a vector, defining the mean outcome in the control group.}
  \item{sd.ctrl}{a vector, defining the standard deviation of the outcome in the control group.}
  \item{n.ctrl}{a vector, defining the number of subjects in the control group.}
  \item{names}{character string identifying each trial.}
  \item{method}{a character string indicating the method to be used. Options are \code{cohens} or \code{hedges} and \code{glass}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\value{
A list containing the following:
  \item{md}{standardised mean difference and its confidence interval computed for each trial.}
  \item{md.invar}{the inverse variance (fixed effects) summary standardised mean difference.}
  \item{md.dsl}{the DerSimonian and Laird (random effects) summary standardised mean difference.}
  \item{heterogeneity}{a vector containing \code{Q} the heterogeneity test statistic, \code{df} the degrees of freedom and its associated P-value.}
}
  
\references{
Deeks JJ, Altman DG, Bradburn MJ (2001). Statistical methods for examining heterogeneity and combining results from several studies in meta-analysis. In: Egger M, Davey Smith G, Altman D (eds). Systematic Review in Health Care Meta-Analysis in Context. British Medical Journal, London, pp. 290 - 291. 
}

\note{
The standardised mean difference method is used when trials assess the same outcome, but measure it in a variety of ways. For example: a set of trials might measure depression scores in psychiatric patients but use different methods to quantify depression. In this circumstance it is necessary to standardise the results of the trials to a uniform scale before they can be combined. The standardised mean difference method expresses the size of the treatment effect in each trial relative to the variability observed in that trial.
}

\seealso{
\code{\link{epi.dsl}, \link{epi.iv}, \link{epi.mh}}
}

\examples{
## EXAMPLE 1:
## A systematic review comparing assertive community treatment (ACT) for the 
## severely mentally ill was compared to standard care. A systematic review
## comparing ACT to standard care found three trials that assessed mental
## status after 12 months. All three trials used a different scoring system,
## so standardisation is required before they can be compared.

names <- c("Audini", "Morse", "Lehman")
mean.trt <- c(41.4, 0.95, -4.10)
mean.ctrl <- c(42.3, 0.89, -3.80)
sd.trt <- c(14, 0.76, 0.83)
sd.ctrl <- c(12.4, 0.65, 0.87)
n.trt <- c(30, 37, 67)
n.ctrl <- c(28, 35, 58)

epi.smd(mean.trt, sd.trt, n.trt, mean.ctrl, sd.ctrl, n.ctrl, 
   names, method = "cohens", conf.level = 0.95)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
