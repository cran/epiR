\name{epi.supc}

\alias{epi.supc}

\title{
Estimate the sample size for a parallel superiority trial, continuous outcomes
}

\description{
Computes the sample size for a parallel superiority trial with a continuous outcome variable. 
}

\usage{
epi.supc(treat, control, sd, delta, n, r = 1, power, alpha)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sd}{the expected population standard deviation of the outcome of interest.}
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
Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, page 61.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Pocock SJ (1983). Clinical Trials: A Practical Approach. Wiley, New York.
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
## EXAMPLE 1
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease (CHD) using a parallel design. The primary 
## efficacy parameter is the concentration of high density lipoproteins.
## (HDL). We consider the situation where the intended trial is to test 
## superiority of the test drug over the active control agent. Sample 
## size calculations are to be calculated to achieve 80\% power at the
## 5\% level of significance.

## In this example, we assume that if treatment results in a 5 unit 
## (i.e. delta = 5) increase in HDL it is declared to be superior to the
## active control. Assume the standard deviation of HDL is 10 units and 
## the HDL concentration in the treatment group is 20 units and the 
## HDL concentration in the control group is 20 units.

epi.supc(treat = 20, control = 20, sd = 10, delta = 5, n = NA, 
   r = 1, power = 0.80, alpha = 0.05)

## A total of 100 subjects need to be enrolled in the trial, 50 in the 
## treatment group and 50 in the control group.
}

\keyword{univar}
