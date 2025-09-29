\name{epi.ssninfc}

\alias{epi.ssninfc}

\title{
Sample size for a non-inferiority trial, continuous outcome
}

\description{
Sample size for a non-inferiority trial, continuous outcome. 
}

\usage{
epi.ssninfc(treat, control, sigma, delta, n, power, r = 1, nfractional = FALSE, 
   alpha = 0.05)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sigma}{the expected population standard deviation of the outcome of interest.}
  \item{delta}{the equivalence limit, expressed as the absolute change in the outcome of interest that represents a clinically meaningful difference. For a non-inferiority trial the value entered for \code{delta} must be greater than or equal to zero.}
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
Blackwelder WC (1982). Proving the null hypothesis in clinical trials. Controlled Clinical Trials 3: 345 - 353. 

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2018). Sample Sizes for Clinical, Laboratory and Epidemiology Studies. Wiley Blackwell, New York.

Scott IA (2009). Non-inferiority trials: determining whether alternative treatments are good enough. Medical Journal of Australia 190: 326 - 330.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.

Zhong B (2009). How to calculate sample size in randomized controlled trial? Journal of Thoracic Disease 1: 51 - 54.
}

\author{
Many thanks to Aniko Szabo (Medical College of Wisconsin, Wisconsin USA) for improvements to the power calculations for this function and suggestions to improve the documentation.
}

\note{
The aim of a non-inferiority trial is show that a new treatment is not (much) inferior to a standard treatment. Showing non-inferiority can be of interest because: (a) it is often not ethically possible to do a placebo-controlled trial; (b) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is safer; (c) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is cheaper to produce or easier to administer; and (d) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints in clinical trial, but compliance will be better outside the clinical trial and hence efficacy better outside the trial.

Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). In each group, the mean of the outcome of interest for subjects receiving the standard treatment is \eqn{N_{s}} and the mean of the outcome of interest for subjects receiving the new treatment is \eqn{N_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{N_{n}} and \eqn{N_{s}} as \eqn{\delta}. For a non-inferiority trial the value entered for \code{delta} must be greater than or equal to zero.

For a non-inferiority trial the null hypothesis is:

\eqn{H_{0}: N_{n} - N_{s} \leq -\delta}
                                 
and the alternative hypothesis is:

\eqn{H_{1}: N_{n} - N_{s} > -\delta}

If the alternative hypothesis is accepted we conclude that \eqn{N_{n}} is not inferior to \eqn{N_{s}}. For a non-inferiority trial a one-sided test is used.

When calculating study power the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

Refer to the documentation for \code{\link{epi.ssequb}} for a comparison of the key features of superiority, equivalence and non-inferiority trials.}

\examples{
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 61 - 62):
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease using a parallel design. The primary efficacy 
## parameter is the LDL. In what follows, we will consider the situation
## where the intended trial is for testing non-inferiority of mean responses 
## in LDL. Assume that 80\% power is required at a 5\% level of significance.

## In this example we assume a -0.05 unit change in LDL is a clinically 
## meaningful difference. Assume the standard deviation of LDL is 0.10 units 
## and the LDL concentration in the treatment group is 0.20 units and the 
## LDL concentration in the control group is 0.20 units.

epi.ssninfc(treat = 0.20, control = 0.20, sigma = 0.10, delta = 0.05, n = NA, 
   power = 0.80, r = 1, nfractional = FALSE, alpha = 0.05)

## A total of 100 subjects need to be enrolled in the trial, 50 in the 
## treatment group and 50 in the control group.
}

\keyword{univar}
