\name{epi.mh}

\alias{epi.mh}


\title{Fixed-effects meta-analysis of binary outcomes using the Mantel-Haenszel method}

\description{
Computes odds ratios (or risk ratios) and confidence intervals for binary outcome data. Computes the pooled odds ratio (or risk ratio) using the Mantel-Haenszel method. Performs a test of heterogeneity among trials. Performs a test for the overall difference between groups.
}

\usage{
epi.mh(ev.trt, n.trt, ev.ctrl, n.ctrl, names, method = "odds.ratio", 
   conf.level = 0.95)
}

\arguments{
  \item{ev.trt}{observed number of events in the treatment group.}
  \item{n.trt}{number in the treatment group.}
  \item{ev.ctrl}{observed number of events in the control group.}
  \item{n.ctrl}{number in the control group.}
  \item{names}{character string identifying each trial.}
  \item{method}{a character string indicating the method to be used. Options are \code{odds.ratio} or \code{risk.ratio}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Under the random-effects model, the assumption of a common treatment effect is relaxed, and the effect sizes are assumed to have a normal distribution with variance \code{tau.sq}. 
}

\value{
A list containing the following:
  \item{odds.ratio}{the names of the trials, the odds ratio for each trial, the standard error of the odds ratio for each trial, the lower and upper bounds of the confidence interval of the odds ratio for each trial, the pooled odds ratio, the standard error of the pooled odds ratio, and the lower and upper bounds of the confidence interval of the pooled odds ratio.}
  \item{risk.ratio}{the names of the trials, the risk ratio for each trial, the standard error of the risk ratio for each trial, the lower and upper bounds of the confidence interval of the risk ratio for each trial, the pooled risk ratio, the standard error of the pooled risk ratio, and the lower and upper bounds of the confidence interval of the pooled risk ratio.}
  \item{heterogeneity}{a vector containing \code{Q} the heterogeneity test statistic, \code{df} the degrees of freedom and its associated P-value.}
  \item{Hsq}{the relative excess of the  heterogeneity test statistic \code{Q} over the degrees of freedom \code{df}.}
  \item{Isq}{the percentage of total variation in study estimates that is due to heterogeneity rather than chance.} 
  \item{effect}{a vector containing \code{z} the test statistic for overall treatment effect and its associated P-value.}
}

\references{
Higgins JP, Thompson SG (2002). Quantifying heterogeneity in a meta-analysis. Statistics in Medicine 21: 1539 - 1558. 
}

\author{
}

\note{
Using this method, the raw weights are used to compute the pooled odds ratios and risk ratios.
}

\seealso{
\code{\link{epi.dsl}, \link{epi.iv}, \link{epi.smd}}
}

\examples{
data(epi.epidural)
epi.mh(ev.trt = epi.epidural$ev.trt, n.trt = epi.epidural$n.trt, 
   ev.ctrl = epi.epidural$ev.ctrl, n.ctrl = epi.epidural$n.ctrl, 
   names = as.character(epi.epidural$trial), method = "odds.ratio",
   conf.level = 0.05)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
