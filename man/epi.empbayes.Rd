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
  \item{obs}{a vector representing the observed disease counts in each region of interest.}
  \item{pop}{a vector representing the population count in each region of interest.}
}

\details{
The gamma distribution is sometimes parameterised in terms of shape and rate parameters. The rate parameter equals the inverse of the scale parameter. The mean of the distribution equals \eqn{ \delta / \alpha}. The variance of the distribution equals \eqn{ \delta / \alpha^{2}}. The empirical Bayes estimate of the proportion affected in each area equals \eqn{(obs + \delta) / (pop + \alpha)}.
}

\value{
A data frame with four elements: \code{gamma}: mean observed event count, \code{phi}: variance of observed event count, \code{alpha}: shape parameter of gamma distribution, and \code{delta}: scale parameter of gamma distribution.
}

\references{
Bailey TC,  Gatrell AC (1995). Interactive Spatial Data Analysis. Longman Scientific & Technical. London. 

Langford IH (1994). Using empirical Bayes estimates in the geographical analysis of disease risk. Area 26: 142 - 149. 
}

\author{
}

\note{
}

\seealso{
}

\examples{
data(epi.SClip)
obs <- epi.SClip$cases
pop <- epi.SClip$population

est <- epi.empbayes(obs, pop)
empbayes.prop <- (obs + est[4]) / (pop + est[3])
raw.prop <- (obs) / (pop)
rank <- rank(raw.prop)
dat <- as.data.frame(cbind(rank, raw.prop, empbayes.prop))

plot(dat$rank, dat$raw.prop, type = "n", xlab = "Rank", ylab = "Proportion")
points(dat$rank, dat$raw.prop, pch = 16, col = "red")
points(dat$rank, dat$empbayes.prop, pch = 16, col = "blue")
legend(5, 0.00025, legend = c("Raw estimate", "Bayes adjusted estimate"), 
   col = c("red","blue"), pch = c(16,16), bty = "n")
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
