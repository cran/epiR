\name{epi.empbayes}

\alias{epi.empbayes}

\title{Empirical Bayes estimates}

\description{
Computes empirical Bayes estimates of observed event counts using the method of moments.
}

\usage{
epi.empbayes(obs, pop)
}

\arguments{
  \item{obs}{a vector representing the observed event counts in each unit of interest.}
  \item{pop}{a vector representing the population count in each unit of interest.}
}

\details{
The gamma distribution is parameterised in terms of shape (\eqn{\alpha}) and scale (\eqn{\nu}) parameters. The mean of a given gamma distribution equals \eqn{\nu / \alpha}. The variance equals \eqn{\nu / \alpha^{2}}. The empirical Bayes estimate of event risk in each unit of interest equals \eqn{(obs + \nu) / (pop + \alpha)}.

This technique performs poorly when your data contains large numbers of zero event counts. In this situation a Bayesian approach for estimating \eqn{\alpha} and \eqn{\nu} would be advised.
}

\value{
A data frame with four elements: \code{gamma} the mean event risk across all units, \code{phi} the variance of event risk across all units, \code{alpha} the estimated shape parameter of the gamma distribution, and \code{nu} the estimated scale parameter of the gamma distribution.
}

\references{
Bailey TC,  Gatrell AC (1995). Interactive Spatial Data Analysis. Longman Scientific & Technical. London, pp. 303 - 308. 

Langford IH (1994). Using empirical Bayes estimates in the geographical analysis of disease risk. Area 26: 142 - 149. 

Meza J (2003). Empirical Bayes estimation smoothing of relative risks in disease mapping. Journal of Statistical Planning and Inference 112: 43 - 62.
}

\examples{
data(epi.SClip)
obs <- epi.SClip$cases; pop <- epi.SClip$population

est <- epi.empbayes(obs, pop)
empbayes.prop <- (obs + est[4]) / (pop + est[3])
raw.prop <- (obs) / (pop)
rank <- rank(raw.prop)
dat <- data.frame(rank, raw.prop, empbayes.prop)

plot(dat$rank, dat$raw.prop, type = "n", xlab = "Rank", ylab = "Risk")
points(dat$rank, dat$raw.prop, pch = 16, col = "red")
points(dat$rank, dat$empbayes.prop, pch = 16, col = "blue")
legend(x = "topleft", legend = c("Raw estimate", "Bayes adjusted estimate"), 
   col = c("red","blue"), pch = c(16,16), bty = "n")
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
