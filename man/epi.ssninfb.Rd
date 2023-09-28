\name{epi.ssninfb}

\alias{epi.ssninfb}

\title{
Sample size for a non-inferiority trial, binary outcome
}

\description{
Sample size for a non-inferiority trial, binary outcome. 
}

\usage{
epi.ssninfb(treat, control, delta, n, power, r = 1, nfractional = FALSE, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
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

Scott IA (2009). Non-inferiority trials: determining whether alternative treatments are good enough. Medical Journal of Australia 190: 326 - 330.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.

Zhong B (2009). How to calculate sample size in randomized controlled trial? Journal of Thoracic Disease 1: 51 - 54.
}

\author{
Many thanks to Aniko Szabo (Medical College of Wisconsin, Wisconsin USA) for improvements to the power calculations for this function and suggestions to improve the documentation.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). A proportion of subjects in the standard treatment group experience the outcome of interest \eqn{P_{s}} and a proportion of subjects in the new treatment group experience the outcome of interest \eqn{P_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. For a non-inferiority trial the value entered for \code{delta} must be greater than or equal to zero.

For a non-inferiority trial the null hypothesis is:

\eqn{H_{0}: P_{s} - P_{n} \ge \delta} 
 
The alternative hypothesis is:

\eqn{H_{1}: P_{s} - P_{n} < \delta}

The aim of a non-inferiority trial is show that a new treatment is not (much) inferior to a standard treatment. Showing non-inferiority can be of interest because: (a) it is often not ethically possible to do a placebo-controlled trial; (b) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is safer; (c) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is cheaper to produce or easier to administer; and (d) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints in clinical trial, but compliance will be better outside the clinical trial and hence efficacy better outside the trial.

When calculating the power of a study, note that the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

For a comparison of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.
}

\examples{
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 90):
## A pharmaceutical company would like to conduct a clinical trial to
## compare the efficacy of two antimicrobial agents when administered orally 
## to patients with skin infections. Assume the true mean cure rate of the 
## treatment is 0.85 and the true mean cure rate of the control is 0.65. 
## We consider the proportion cured in the treatment group minus the proportion 
## cured in the control group (i.e., delta) of 0.10 or less to be of no clinical 
## significance.

## Assuming a one-sided test size of 5\% and a power of 80\% how many 
## subjects should be included in the trial?

epi.ssninfb(treat = 0.85, control = 0.65, delta = 0.10, n = NA, power = 0.80, 
   r = 1, nfractional = FALSE, alpha = 0.05)

## A total of 50 subjects need to be enrolled in the trial, 25 in the 
## treatment group and 25 in the control group.


## EXAMPLE 1 (cont.):
## Suppose only 40 subjects were enrolled in the trial, 20 in the treatment
## group and 20 in the control group. What is the estimated study power?

epi.ssninfb(treat = 0.85, control = 0.65, delta = 0.10, n = 40, power = NA, 
   r = 1, nfractional = FALSE, alpha = 0.05)

## With only 40 subjects the estimated study power is 0.73.


## EXAMPLE 2:
## Assume the true mean cure rate for a treatment group to be 0.40 and the true
## mean cure rate for a control group to be the same, 0.40. We consider a 
## difference of 0.10 in cured proportions (i.e., delta = 0.10) to be of no 
## clinical importance.

## Assuming a one-sided test size of 5\% and a power of 30\% how many 
## subjects should be included in the trial?

n <- epi.ssninfb(treat = 0.4, control = 0.4, delta = 0.10, n = NA, power = 0.3, 
   r = 1, nfractional = TRUE, alpha = 0.05)$n.total
n

## A total of 120 subjects need to be enrolled in the trial, 60 in the 
## treatment group and 60 in the control group.

## Re-run the function using n = 120 to confirm that power equals 0.30:

epi.ssninfb(treat = 0.4, control = 0.4, delta = 0.10, n = n, power = NA, 
   r = 1, nfractional = TRUE, alpha = 0.05)$power

## With 120 subjects the estimated study power is 0.30.
}

\keyword{univar}
