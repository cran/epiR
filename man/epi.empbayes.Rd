\name{epi.empbayes}

\alias{epi.empbayes}

\title{Empirical Bayes estimates of observed event counts}

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
crude.p <- ((obs) / (pop)) * 100000
crude.r <- rank(crude.p)
ebay.p <- ((obs + est[4]) / (pop + est[3])) * 100000
dat <- data.frame(rank = c(crude.r, crude.r), 
   Method = c(rep("Crude", times = length(crude.r)), 
      rep("Empirical Bayes", times = length(crude.r))),
   est = c(crude.p, ebay.p)) 

## Scatter plot showing the crude and empirical Bayes adjusted lip cancer 
## incidence rates as a function of district rank for the crude lip 
## cancer incidence rates: 
                          
\dontrun{
library(ggplot2)

ggplot(dat, aes(x = rank, y = est, colour = Method)) + 
   geom_point() +
   scale_x_continuous(name = "District rank", 
      breaks = seq(from = 0, to = 60, by = 10), 
      labels = seq(from = 0, to = 60, by = 10), 
      limits = c(0,60)) +
  scale_y_continuous(limits = c(0,30), name = "Lip cancer incidence rates 
      (cases per 100,000 person years)") 
}
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
