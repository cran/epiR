\name{epi.ssequc}

\alias{epi.ssequc}

\title{
Sample size for a parallel equivalence trial, continuous outcome
}

\description{
Sample size for a parallel equivalence trial, continuous outcome. 
}

\usage{
epi.ssequc(treat, control, sd, delta, n, r = 1, power, nfractional = FALSE, 
   alpha)
}

\arguments{
  \item{treat}{the expected mean of the outcome of interest in the treatment group.}
  \item{control}{the expected mean of the outcome of interest in the control group.}
  \item{sd}{the expected population standard deviation of the outcome of interest.}
  \item{delta}{the equivalence limit, expressed as the absolute change in the outcome of interest that represents a clinically meaningful difference. For an equivalence trial the value entered for \code{delta} must be greater than zero.}
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
Bennett JE, Dismukes WE, Duma RJ, Medoff G, Sande MA, Gallis H, Leonard J, Fields BT, Bradshaw M, Haywood H, McGee Z, Cate TR, Cobbs CG, Warner JF and Alling DW (1979). A comparison of amphotericin B alone and combined flucytosine in the treatment of cryptococcal meningitis. New England Journal of Medicine 301: 126 - 131.

Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, pp. 91.

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). In each group, the mean of the outcome of interest for subjects receiving the standard treatment is \eqn{N_{s}} and the mean of the outcome of interest for subjects receiving the new treatment is \eqn{N_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{N_{n}} and \eqn{N_{s}} as \eqn{\delta}. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero.

For an equivalence trial the null hypothesis is:

\eqn{H_{0}: |N_{s} - N_{n}| \ge \delta} 
 
The alternative hypothesis is:

\eqn{H_{1}: |N_{s} - N_{n}| < \delta}

An equivalence trial is used if want to prove that two treatments produce the same clinical outcomes. In bioequivalence trials, a 90\% confidence interval is often used. The value of the maximum acceptable difference \eqn{\delta} is chosen so that a patient will not detect any change in effect when replacing the standard treatment with the new treatment.

Note that when:

\deqn{sign(P_n - P_s - \delta) \neq sign(z_{1-\alpha} + z_{1-\beta})}

there is no solution for study power. For typical values of \eqn{\alpha} and \eqn{\beta} this would occur if \eqn{P_{n} - P_{s} - \delta < 0}. That is, when the targeted alternative is within the null hypothesis. The function issues a warning if these conditions are met.

When calculating the power of a study, the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

For a comparison of the key features of superiority, equivalence and non-inferiority trials, refer to the documentation for \code{\link{epi.ssequb}}.
}

\examples{
## EXAMPLE 1 (from Machin, Campbell, Tan and Tan 2009 p. 113):
## It is anticipated that patients on a particular drug have a mean diastolic 
## blood pressure of 96 mmHg, as against 94 mmHg on an alternative. It is also 
## anticipated that the standard deviation of diastolic BP is approximately 
## 8 mmHg. If one wishes to confirm that the difference is likely to be less 
## than 5 mmHg, that is, one wishes to show equivalence, how many patients 
## are needed to be enrolled in the trial? Assume 80\% power and 
## 95\% significance.

epi.ssequc(treat = 94, control = 96, sd = 8, delta = 5, n = NA, 
   r = 1, power = 0.80, nfractional = FALSE, alpha = 0.05)

## A total of 244 subjects need to be enrolled in the trial, 122 in the 
## treatment group and 122 in the control group.


## EXAMPLE 2 (from Chow S, Shao J, Wang H 2008, p. 64):
## A pharmaceutical company is interested in conducting a clinical trial
## to compare two cholesterol lowering agents for treatment of patients with
## congestive heart disease using a parallel design. The primary efficacy 
## parameter is the LDL. In what follows, we will consider the situation
## where the intended trial is for testing equivalence of mean responses 
## in LDL. Assume that 80\% power is required at a 5\% level of significance.

## In this example, we assume a 5 unit (i.e., delta = 5) change of LDL is 
## considered of clinically meaningful difference. Assume the standard 
## of LDL is 10 units and the LDL concentration in the treatment group is 20 
## units and the LDL concentration in the control group is 21 units.

epi.ssequc(treat = 20, control = 21, sd = 10, delta = 5, n = NA, 
   r = 1, power = 0.80, nfractional = FALSE, alpha = 0.05)

## A total of 216 subjects need to be enrolled in the trial, 108 in the 
## treatment group and 108 in the control group.


## EXAMPLE 2 (cont.):
## Suppose only 150 subjects were enrolled in the trial, 75 in the treatment
## group and 75 in the control group. What is the estimated study power?

epi.ssequc(treat = 20, control = 21, sd = 10, delta = 5, n = 150, 
   r = 1, power = NA, nfractional = FALSE, alpha = 0.05)

## With only 150 subjects enrolled the estimated study power is 0.58.
}

\keyword{univar}
