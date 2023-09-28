\name{epi.sssupc}

\alias{epi.sssupc}

\title{
Sample size for a parallel superiority trial, continuous outcome
}

\description{
Sample size for a parallel superiority trial, continuous outcome. 
}

\usage{
epi.sssupc(treat, control, sigma, delta, n, power, r = 1, nfractional = FALSE, 
   alpha)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sigma}{the expected population standard deviation of the outcome of interest.}
  \item{delta}{the equivalence limit, expressed as the absolute change in the outcome of interest that represents a clinically meaningful difference. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero.}
  \item{n}{scalar, the total number of study subjects in the trial.}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
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
Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, page 61.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Pocock SJ (1983). Clinical Trials: A Practical Approach. Wiley, New York.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). In each group, the mean of the outcome of interest for subjects receiving the standard treatment is \eqn{N_{s}} and the mean of the outcome of interest for subjects receiving the new treatment is \eqn{N_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{N_{n}} and \eqn{N_{s}} as \eqn{\delta}. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero.

For a superiority trial the null hypothesis is:

\eqn{H_{0}: N_{s} - N_{n} = 0} 
 
The alternative hypothesis is:

\eqn{H_{1}: N_{s} - N_{n} != 0}

When calculating the power of a study, the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

For a comparison of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.}

\examples{
## EXAMPLE 1:
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease (CHD) using a parallel design. The primary 
## efficacy parameter is the concentration of high density lipoproteins
## (HDL). We consider the situation where the intended trial is to test 
## superiority of the test drug over the active control agent. Sample 
## size calculations are to be calculated to achieve 80\% power at the
## 5\% level of significance.

## In this example, we assume that if treatment results in a 5 unit 
## (i.e., delta = 5) increase in HDL it is declared to be superior to the
## active control. Assume the standard deviation of HDL is 10 units and 
## the HDL concentration in the treatment group is 20 units and the 
## HDL concentration in the control group is 20 units.

epi.sssupc(treat = 20, control = 20, sigma = 10, delta = 5, n = NA, 
   power = 0.80, r = 1, nfractional = FALSE, alpha = 0.05)

## A total of 100 subjects need to be enrolled in the trial, 50 in the 
## treatment group and 50 in the control group.
}

\keyword{univar}
