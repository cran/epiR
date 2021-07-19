\name{epi.ssninfb}

\alias{epi.ssninfb}

\title{
Sample size for a non-inferiority trial, binary outcome
}

\description{
Sample size for a non-inferiority trial, binary outcome. 
}

\usage{
epi.ssninfb(treat, control, delta, n, r = 1, power, nfractional = FALSE, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
  \item{delta}{the equivalence limit, expressed as the change in the outcome of interest that represents a clinically meaningful diffference.}
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
  \item{power}{the specified or calculated study power.}
}

\references{
Blackwelder WC (1982). Proving the null hypothesis in clinical trials. Controlled Clinical Trials 3: 345 - 353. 

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Scott IA (2009). Non-inferiority trials: determining whether alternative treatments are good enough. Medical Journal of Australia 190: 326 - 330.

Zhong B (2009). How to calculate sample size in randomized controlled trial? Journal of Thoracic Disease 1: 51 - 54.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (s) and a new treatment (n). In each group, a proportion of subjects respond to the treatment: \eqn{P_{s}} and \eqn{P_{n}}. 

With a superiority trial we specify the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. The null hypothesis is \eqn{H_{0}: P_{n} - P_{s} \leq \delta} and the alternative hypothesis is \eqn{H_{1}: P_{n} - P_{s} > \delta}.

An equivalence trial is used if want to prove that two treatments produce the same clinical outcomes. With an equivalence trial, we specify the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. The null hypothesis is \eqn{H_{0}: |P_{s} - P_{n}| \geq \delta} and the alternative hypothesis is \eqn{H_{1}: |P_{s} - P_{n}| < \delta}. In bioequivalence trials, a 90\% confidence interval is often used. The value of the maximum acceptable difference \eqn{\delta} is chosen so that a patient will not detect any change in effect when replacing the standard treatment with the new treatment.

With a non-inferiority trial, we specify the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. The null hypothesis is \eqn{H_{0}: P_{s} - P_{n} \geq \delta} and the alternative hypothesis is \eqn{H_{1}: P_{s} - P_{n} < \delta}. The aim of a non-inferiority trial is show that a new treatment is not (much) inferior to a standard treatment. Showing non-inferiority can be of interest because: (a) it is often not ethically possible to do a placebo-controlled trial, (b) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is safer, (c) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is cheaper to produce or easier to administer, (d) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints in clinical trial, but compliance will be better outside the clinical trial and hence efficacy better outside the trial.

For a summary of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.

When calculating the power of a study, note that the variable \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).
}

\examples{
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 90):
## A pharmaceutical company would like to conduct a clinical trial to
## compare the efficacy of two antimicrobial agents when administered orally 
## to patients with skin infections. Assume the true mean cure rate of the 
## treatment is 0.85 and the true mean cure rate of the control is 0.65. 
## We consider a difference of less than 0.10 in cure rate to be of no 
## clinical importance (i.e. delta = -0.10).

## Assuming a one-sided test size of 5\% and a power of 80\% how many 
## subjects should be included in the trial?

epi.ssninfb(treat = 0.85, control = 0.65, delta = -0.10, n = NA, r = 1, 
   power = 0.80, nfractional = FALSE, alpha = 0.05)

## A total of 50 subjects need to be enrolled in the trial, 25 in the 
## treatment group and 25 in the control group.


## EXAMPLE 1 (cont.):
## Suppose only 40 subjects were enrolled in the trial, 20 in the treatment
## group and 20 in the control group. What is the estimated study power?

epi.ssninfb(treat = 0.85, control = 0.65, delta = -0.10, n = 40, r = 1,
   power = NA, nfractional = FALSE, alpha = 0.05)

## With only 40 subjects the estimated study power is 0.73.

}

\keyword{univar}
