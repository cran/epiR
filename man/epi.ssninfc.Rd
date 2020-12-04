\name{epi.ssninfc}

\alias{epi.ssninfc}

\title{
Sample size for a non-inferiority trial, continuous outcome
}

\description{
Sample size for a non-inferiority trial, continuous outcome. 
}

\usage{
epi.ssninfc(treat, control, sd, delta, n, r = 1, power, nfractional = FALSE, 
   alpha)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sd}{the expected population standard deviation of the outcome of interest.}
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
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 61 - 62):
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease using a parallel design. The primary efficacy 
## parameter is the LDL. In what follows, we will consider the situation
## where the intended trial is for testing non-inferiority of mean responses 
## in LDL. Assume that 80\% power is required at a 5\% level of significance.

## In this example, we assume a 5% (i.e. delta = 0.05) change in LDL is 
## a clinically meaningful difference. Assume the standard deviation of 
## LDL is 0.10 and the LDL concentration in the treatment group is 0.20 
## units and the LDL concentration in the control group is 0.20 units.

epi.ssninfc(treat = 0.20, control = 0.20, sd = 0.10, delta = 0.05, n = NA, 
   r = 1, power = 0.80, nfractional = FALSE, alpha = 0.05)

## A total of 100 subjects need to be enrolled in the trial, 50 in the 
## treatment group and 50 in the control group.
}

\keyword{univar}
