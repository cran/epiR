\name{epi.noninfb}

\alias{epi.noninfb}

\title{
Estimate the sample size for a parallel non-inferiority trial, binary outcomes
}

\description{
Computes the sample size for a parallel non-inferiority trial with a binary outcome variable. 
}

\usage{
epi.noninfb(treat, control, delta, n, r = 1, power, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
  \item{delta}{the equivalence limit, expressed as a proportion.}
  \item{n}{scalar, the total number of study subjects in the trial.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{power}{scalar, the required study power.}
  \item{alpha}{scalar, defining the desired alpha level.}
}

\value{
A list containing one or more of the following: 
  \item{n.treat}{the required number of study subject in the treatment group.}
  \item{n.control}{the required number of study subject in the control group.}
  \item{n.total}{the total number of study subjects required.}
}

\references{
Blackwelder WC (1982). Proving the null hypothesis in clinical trials. Controlled Clinical Trials. 3: 345 - 353. 

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber. 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Scott IA (2009). Non-inferiority trials: determining whether alternative treatments are good enough. Medical Journal of Australia. 190: 326 - 330.

Zhong B (2009). How to calculate sample size in randomized controlled trial? Journal of Thoracic Disease. 1: 51 - 54.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (s) and a new treatment (n). In each group, a proportion of subjects respond to the treatment: Ps and Pn. 

With a superiority trial we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: Pn - Ps <= \code{delta} and the alternative hypothesis is H1: Pn - Ps > \code{delta}.

An equivalence trial is used if want to prove that two treatments produce the same clinical outcomes. With an equivalence trial, we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: |Ps - Pn| >= \code{delta} and the alternative hypothesis is H1: |Ps - Pn| < \code{delta}. In bioequivalence trials, a 90\% confidence interval is often used. The value of the maximum acceptable difference \code{delta} is chosen so that a patient will not detect any change in effect when replacing the standard treatment with the new treatment.

With a non-inferiority trial, we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: Ps - Pn >= \code{delta} and the alternative hypothesis is H1: Ps - Pn < \code{delta}. The aim of a non-inferiority trial is show that a new treatment is not (much) inferior to a standard treatment. Showing non-inferiority can be of interest because: (a) it is often not ethically possible to do a placebo-controlled trial, (b) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is safer, (c) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is cheaper to produce or easier to administer, (d) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints in clinical trial, but compliance will be better outside the clinical trial and hence efficacy better outside the trial.

For a summary of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.equivb}}.

When calculating the power of a study, note that the variable \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).
}

\examples{
## EXAMPLE 1:
## Suppose it is of interest to establish non-inferiority of a new treatment 
## as compared to the currently used standard treatment. A difference of less
## than 10\% is of no clinical importance. Thus, the non-inferiority margin 
## (delta) is set at 0.10. Assume the true cure rate for the new treatment
## is 0.85 and the control is 0.65. Assuming a one-sided test size of 2.5\% and 
## a power of 90\% how many subjects should be included in the trial?

epi.noninfb(treat = 0.85, control = 0.65, delta = 0.10, n = NA, r = 1,
   power = 0.80, alpha = 0.025)

## A total of 558 subjects need to be enrolled in the trial, 279 in the 
## treatment group and 279 in the control group.

## EXAMPLE 1 (cont.):
## Suppose only 400 subjects were enrolled in the trial, 200 in the treatment
## group and 200 in the control group. What is the estimated study power?

epi.noninfb(treat = 0.85, control = 0.65, delta = 0.10, n = 400, r = 1,
   power = NA, alpha = 0.025)

## With only 500 subjects the estimated study power is 0.66.

}

\keyword{univar}
