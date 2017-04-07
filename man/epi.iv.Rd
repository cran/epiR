\name{epi.iv}

\alias{epi.iv}

\title{Fixed-effect meta-analysis of binary outcomes using the inverse variance method
}

\description{
Computes individual study odds or risk ratios for binary outcome data. Computes the summary odds or risk ratio using the inverse variance method. Performs a test of heterogeneity among trials. Performs a test for the overall difference between groups (that is, after pooling the studies, do treated groups differ significantly from controls?).
}

\usage{
epi.iv(ev.trt, n.trt, ev.ctrl, n.ctrl, names, method = "odds.ratio", 
   alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
}

\arguments{
  \item{ev.trt}{observed number of events in the treatment group.}
  \item{n.trt}{number in the treatment group.}
  \item{ev.ctrl}{observed number of events in the control group.}
  \item{n.ctrl}{number in the control group.}
  \item{names}{character string identifying each trial.}
  \item{method}{a character string indicating the method to be used. Options are \code{odds.ratio} or \code{risk.ratio}.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{two.sided}, \code{greater} or \code{less}. }
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Using this method, the inverse variance weights are used to compute the pooled odds ratios and risk ratios. The inverse variance weights should be used to indicate the weight each trial contributes to the meta-analysis.

\code{alternative = "greater"} tests the hypothesis that the inverse variance summary measure of association is greater than 1. 

}

\value{
A list containing:
  \item{OR}{the odds ratio for each trial and the lower and upper bounds of the confidence interval of the odds ratio for each trial.}
  \item{RR}{the risk ratio for each trial and the lower and upper bounds of the confidence interval of the risk ratio for each trial.}
  \item{OR.summary}{the inverse variance summary odds ratio and the lower and upper bounds of the confidence interval of the inverse variance summary odds ratio.}
  \item{RR.summary}{the inverse variance summary risk ratio and the lower and upper bounds of the confidence interval of the inverse variance summary risk ratio.}
  \item{weights}{the raw and inverse variance weights assigned to each trial.}
  \item{heterogeneity}{a vector containing \code{Q} the heterogeneity test statistic, \code{df} the degrees of freedom and its associated P-value.}
  \item{Hsq}{the relative excess of the  heterogeneity test statistic \code{Q} over the degrees of freedom \code{df}.}
  \item{Isq}{the percentage of total variation in study estimates that is due to heterogeneity rather than chance.}  
  \item{effect}{a vector containing \code{z} the test statistic for overall treatment effect and its associated P-value.}
  }
  
\references{
Deeks JJ, Altman DG, Bradburn MJ (2001). Statistical methods for examining heterogeneity and combining results from several studies in meta-analysis. In: Egger M, Davey Smith G, Altman D (eds). Systematic Review in Health Care Meta-Analysis in Context. British Medical Journal, London, 2001, pp. 291 - 299.

Higgins JP, Thompson SG (2002). Quantifying heterogeneity in a meta-analysis. Statistics in Medicine 21: 1539 - 1558. 
}

\note{
The inverse variance method performs poorly when data are sparse, both in terms of event rates being low and trials being small. The Mantel-Haenszel method (\code{\link{epi.mh}}) is more robust when data are sparse.

Using this method, the inverse variance weights are used to compute the pooled odds ratios and risk ratios. 

The function checks each strata for cells with zero frequencies. If a zero frequency is found in any cell, 0.5 is added to all cells within the strata. 
}

\seealso{
\code{\link{epi.dsl}}, \code{\link{epi.mh}}, \code{\link{epi.smd}}
}

\examples{
data(epi.epidural)
epi.iv(ev.trt = epi.epidural$ev.trt, n.trt = epi.epidural$n.trt, 
   ev.ctrl = epi.epidural$ev.ctrl, n.ctrl = epi.epidural$n.ctrl, 
   names = as.character(epi.epidural$trial), method = "odds.ratio", 
   alternative = "two.sided", conf.level = 0.95)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
