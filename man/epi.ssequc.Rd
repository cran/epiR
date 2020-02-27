\name{epi.ssequc}

\alias{epi.ssequc}

\title{
Sample size for a parallel equivalence trial, continuous outcome
}

\description{
Sample size for a parallel equivalence trial, continuous outcome. 
}

\usage{
epi.ssequc(treat, control, sd, delta, n, r = 1, power, alpha)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sd}{the expected population standard deviation of the outcome of interest.}
  \item{delta}{the equivalence limit, expressed as the change in the outcome of interest that represents a clinically meaningful diffference.}
  \item{n}{scalar, the total number of study subjects in the trial.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{power}{scalar, the required study power.}
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
Bennett JE, Dismukes WE, Duma RJ, Medoff G, Sande MA, Gallis H, Leonard J, Fields BT, Bradshaw M, Haywood H, McGee Z, Cate TR, Cobbs CG, Warner JF and Alling DW (1979). A comparison of amphotericin B alone and combined flucytosine in the treatment of cryptococcal meningitis. New England Journal of Medicine 301: 126 - 131.

Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, pp. 91.

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (s) and a new treatment (n). In each group, a proportion of subjects respond to the treatment: Ps and Pn. 

With a superiority trial we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: Pn - Ps <= \code{delta} and the alternative hypothesis is H1: Pn - Ps > \code{delta}.

An equivalence trial is used if want to prove that two treatments produce the same clinical outcomes. With an equivalence trial, we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: |Ps - Pn| >= \code{delta} and the alternative hypothesis is H1: |Ps - Pn| < \code{delta}. In bioequivalence trials, a 90\% confidence interval is often used. The value of the maximum acceptable difference \code{delta} is chosen so that a patient will not detect any change in effect when replacing the standard treatment with the new treatment.

With a non-inferiority trial, we specify the maximum acceptable difference between Pn and Ps as \code{delta}. The null hypothesis is H0: Ps - Pn >= \code{delta} and the alternative hypothesis is H1: Ps - Pn < \code{delta}. The aim of a non-inferiority trial is show that a new treatment is not (much) inferior to a standard treatment. Showing non-inferiority can be of interest because: (a) it is often not ethically possible to do a placebo-controlled trial, (b) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is safer, (c) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints, but is cheaper to produce or easier to administer, (d) the new treatment is not expected to be better than the standard treatment on primary efficacy endpoints in clinical trial, but compliance will be better outside the clinical trial and hence efficacy better outside the trial.

For a summary of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.

When calculating the power of a study, note that the variable \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).
}

\examples{
## EXAMPLE 1 (from Machin, Campbell, Tan and Tan 2009 p. 113):
## It is anticipated that patients on a particular drug have a mean diastolic 
## blood pressure of 96 mmHg, as against 94 mmHg on an alternative. It is also 
## anticipated that the standard deviation of diastolic BP is approximately 
## 8 mmHg. If one wishes to confirm that the difference is likely to be less 
## than 5 mmHg, that is, one wishes to show equivalence, how many patients 
## are need to be enrolled in the trial? Assume 80\% power and 
## 95\% significance.

epi.ssequc(treat = 94, control = 96, sd = 8, delta = 5, n = NA, 
   r = 1, power = 0.80, alpha = 0.05)

## A total of 244 subjects need to be enrolled in the trial, 122 in the 
## treatment group and 122 in the control group.

## EXAMPLE 2 (from Chow S, Shao J, Wang H 2008, p. 64):
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease using a parallel design. The primary efficacy 
## parameter is the LDL. In what follows, we will consider the situation
## where the intended trial is for testing equivalence of mean responses 
## in LDL. Assume that 80\% power is required at a 5\% level of significance.

## In this example, we assume a 5 unit (i.e. delta = 5) change of LDL is 
## considered of clinically meaningful difference. Assume the standard 
## of LDL is 10 units and the LDL concentration in the treatment group is 20 
## units and the LDL concentration in the control group is 21 units.

epi.ssequc(treat = 20, control = 21, sd = 10, delta = 5, n = NA, 
   r = 1, power = 0.80, alpha = 0.05)

## A total of 216 subjects need to be enrolled in the trial, 108 in the 
## treatment group and 108 in the control group.

## EXAMPLE 2 (cont.):
## Suppose only 150 subjects were enrolled in the trial, 75 in the treatment
## group and 75 in the control group. What is the estimated study power?

epi.ssequc(treat = 0.20, control = 0.21, sd = 0.10, delta = 0.05, n = 150, 
   r = 1, power = NA, alpha = 0.05)

## With only 150 subjects the estimated study power is 0.58.
}

\keyword{univar}
