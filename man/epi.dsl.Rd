\name{epi.dsl}

\alias{epi.dsl}

\title{Mixed-effects meta-analysis of binary outcomes using the DerSimonian and Laird method
}

\description{
Computes individual study odds or risk ratios for binary outcome data. Computes the summary odds or risk ratio using the DerSimonian and Laird method. Performs a test of heterogeneity among trials. Performs a test for the overall difference between groups (that is, after pooling the studies, do treated groups differ significantly from controls?).
}

\usage{
epi.dsl(ev.trt, n.trt, ev.ctrl, n.ctrl, names, method = "odds.ratio", 
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
}

\value{
A list containing the following:
  \item{OR}{the odds ratio for each trial, the standard error of the odds ratio for each trial, and the lower and upper bounds of the confidence interval of the odds ratio for each trial.}
  \item{RR}{the risk ratio for each trial, the standard error of the risk ratio for each trial, and the lower and upper bounds of the confidence interval of the risk ratio for each trial.}
  \item{OR.summary}{the DerSimonian and Laird summary odds ratio, the standard error of the DerSimonian and Laird summary odds ratio, the lower and upper bounds of the confidence interval of the DerSimonian and Laird summary odds ratio.}
  \item{RR.summary}{the DerSimonian and Laird summary risk ratio, the standard error of the DerSimonian and Laird summary risk ratio, the lower and upper bounds of the confidence interval of the DerSimonian and Laird summary risk ratio.}
  \item{weights}{the inverse variance and DerSimonian and Laird weights for each trial.}
  \item{heterogeneity}{a vector containing \code{Q} the heterogeneity test statistic, \code{df} the degrees of freedom and its associated P-value.}
  \item{Hsq}{the relative excess of the  heterogeneity test statistic \code{Q} over the degrees of freedom \code{df}.}
  \item{Isq}{the percentage of total variation in study estimates that is due to heterogeneity rather than chance.} 
  \item{tau.sq}{the variance of the treatment effect among trials.}
  \item{effect}{a vector containing \code{z} the test statistic for overall treatment effect and its associated P-value.}
}

\references{
Deeks JJ, Altman DG, Bradburn MJ (2001). Statistical methods for examining heterogeneity and combining results from several studies in meta-analysis. In: Egger M, Davey Smith G, Altman D (eds). Systematic Review in Health Care Meta-Analysis in Context. British Medical Journal, London, 2001, pp. 291 - 299.

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled Clinical Trials 7: 177 - 188. 

Higgins J, Thompson S (2002). Quantifying heterogeneity in a meta-analysis. Statistics in Medicine 21: 1539 - 1558. 
}

\author{
}

\note{
Under the random-effects model, the assumption of a common treatment effect is relaxed, and the effect sizes are assumed to have a normal distribution with variance \code{tau.sq}.

Using this method, the DerSimonian and Laird weights are used to compute the pooled odds ratio.

The function checks each strata for cells with zero frequencies. If a zero frequency is found in any cell, 0.5 is added to all cells within the strata.
}

\seealso{
\code{\link{epi.iv}}, \code{\link{epi.mh}}, \code{\link{epi.smd}}
}

\examples{
data(epi.epidural)
epi.dsl(ev.trt = epi.epidural$ev.trt, n.trt = epi.epidural$n.trt, 
   ev.ctrl = epi.epidural$ev.ctrl, n.ctrl = epi.epidural$n.ctrl, 
   names = as.character(epi.epidural$trial), method = "odds.ratio", 
   conf.level = 0.95)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
