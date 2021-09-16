\name{epi.sssupb}

\alias{epi.sssupb}

\title{
Sample size for a parallel superiority trial, binary outcome
}

\description{
Sample size for a parallel superiority trial, binary outcome. 
}

\usage{
epi.sssupb(treat, control, delta, n, r = 1, power, nfractional = FALSE, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
  \item{delta}{the equivalence limit, expressed as the absolute change in the outcome of interest that represents a clinically meaningful difference. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero.}
  \item{n}{scalar, the total number of study subjects in the trial.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{power}{scalar, the required study power.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{alpha}{scalar, defining the desired alpha level.}
}

\value{
A list containing the following: 
  \item{n.total}{the total number of study subjects required.}
  \item{n.treat}{the required number of study subject in the treatment group.}
  \item{n.control}{the required number of study subject in the control group.}
  \item{delta}{the equivalence limit, as entered by the user.}  
  \item{power}{the specified or calculated study power.}
}

\references{
Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, page 90.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Pocock SJ (1983). Clinical Trials: A Practical Approach. Wiley, New York.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). A proportion of subjects in the standard treatment group experience the outcome of interest \eqn{P_{s}} and a proportion of subjects in the new treatment group experience the outcome of interest \eqn{P_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero. 

For a superiority trial the null hypothesis is:

\eqn{H_{0}: P_{s} - P_{n} = 0} 
 
The alternative hypothesis is:

\eqn{H_{1}: P_{s} - P_{n} != 0}

When calculating the power of a study, the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

For a comparison of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.
}

\examples{
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 91):
## Suppose that a pharmaceutical company is interested in conducting a
## clinical trial to compare the efficacy of two antimicrobial agents 
## when administered orally once daily in the treatment of patients 
## with skin infections. In what follows, we consider the situation 
## where the intended trial is for testing superiority of the 
## test drug over the active control drug. For this purpose, the following 
## assumptions are made. First, sample size calculation will be performed 
## for achieving 80\% power at the 5\% level of significance.

## Assume the true mean cure rates of the treatment agents and the active 
## control are 85\% and 65\%, respectively. Assume the superiority
## margin is 5\%.

epi.sssupb(treat = 0.85, control = 0.65, delta = 0.05, n = NA, 
   r = 1, power = 0.80, nfractional = FALSE, alpha = 0.05)

## A total of 196 subjects need to be enrolled in the trial, 98 in the 
## treatment group and 98 in the control group.

}

\keyword{univar}
