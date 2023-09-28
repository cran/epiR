\name{epi.ssequb}

\alias{epi.ssequb}

\title{
Sample size for a parallel equivalence or equality trial, binary outcome
}

\description{
Sample size for a parallel equivalence or equality trial, binary outcome. 
}

\usage{
epi.ssequb(treat, control, delta, n, power, r = 1, type = "equivalence", 
   nfractional = FALSE, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
  \item{delta}{the equivalence limit, expressed as the absolute change in the outcome of interest that represents a clinically meaningful difference. For an equivalence trial the value entered for \code{delta} must be greater than zero.}
  \item{n}{scalar, the total number of study subjects in the trial.}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{type}{a character string indicating the type of analysis to be run. Options are \code{equality} and \code{equivalence}. If \code{type = equality} is selected argument \code{delta} is ignored.}
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
Bennett J, Dismukes W, Duma R, Medoff G, Sande M, Gallis H, Leonard J, Fields B, Bradshaw M, Haywood H, McGee Z, Cate T, Cobbs C, Warner J, Alling D (1979). A comparison of amphotericin B alone and combined with flucytosine in the treatment of cryptoccal meningitis. New England Journal of Medicine 301, 126 - 131. DOI: 10.1056/NEJM197907193010303.

Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, pp. 91.

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.
}

\note{
A summary of the key features of superiority, equivalence and non-superiority trial comparisons is shown in the following table (adapted from Campbell et al., 2018 [page 170] and Wang et al., 2017):

\tabular{lllll}{
--------------- \tab ------------------------- \tab ------------------------ \tab ---------------------- \tab ---------------------- \cr
Test for        \tab Null hypothesis           \tab Alt hypothesis           \tab Type I error           \tab Type II error          \cr
--------------- \tab ------------------------- \tab ------------------------ \tab ---------------------- \tab ---------------------- \cr
Superiority     \tab H0:  Ps - Pn  == 0        \tab H1:  Ps - Pn  != 0       \tab 2 sided, 0.050         \tab 1 sided, 0.10 or 0.20  \cr
Equivalence     \tab H0: |Ps - Pn| >= delta    \tab H1: |Ps - Pn| < delta    \tab 2 sided, 0.050         \tab 2 sided, 0.10 or 0.20  \cr
Non-inferiority \tab H0:  Ps - Pn  >= delta    \tab H1:  Ps - Pn  < delta    \tab 1 sided, 0.050         \tab 1 sided, 0.10 or 0.20  \cr
--------------- \tab ------------------------- \tab ------------------------ \tab ---------------------- \tab ---------------------- \cr
}

With a superiority trial the aim is to estimate \eqn{Ps - Pn} with a view to claiming a difference between groups.

With an equivalence trial the aim is not to estimate \eqn{Ps - Pn} but to judge if it is within the margins defined by \code{delta}.

With a non-inferiority trial the aim is not to estimate \eqn{Ps - Pn} but to judge if it is within the margins defined by \code{delta}.

Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). A proportion of subjects in the standard treatment group experience the outcome of interest \eqn{P_{s}} and a proportion of subjects in the new treatment group experience the outcome of interest \eqn{P_{n}}. We specify the absolute value of the maximum acceptable difference between \eqn{P_{n}} and \eqn{P_{s}} as \eqn{\delta}. 

In an equality trial we test the hypothesis that the means of two (treatment, control) groups are different. For an equality trial the null hypothesis is:

\eqn{H_{0}: |P_{s} - P_{n}| = 0} 
 
and the alternative hypothesis is:

\eqn{H_{1}: |P_{s} - P_{n}| \neq 0}

In an equivalence trial we test the hypothesis that the means of two (treatment, control) groups are not too different, where clinical criteria are used to define 'not too different' (\code{delta}). For an equivalence trial the null hypothesis is:

\eqn{H_{0}: |P_{s} - P_{n}| \ge \delta} 
 
The alternative hypothesis is:

\eqn{H_{1}: |P_{s} - P_{n}| < \delta}

An equality trial is used when you want to confirm that the numeric value of the outcome for treatment and control groups differ. An equivalence trial is used when you want to prove that the clinical outcome for treatment and control groups are the same. For a superiority trial the value entered for \code{delta} must be greater than or equal to zero. Note that when:

\deqn{sign(P_n - P_s - \delta) \neq sign(z_{1-\alpha} + z_{1-\beta})}

there is no solution for study power. For typical values of \eqn{\alpha} and \eqn{\beta} this would occur if \eqn{P_{n} - P_{s} - \delta < 0}. That is, when the targeted alternative is within the null hypothesis. The function issues a warning if these conditions are met.

When calculating the power of a study, the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).
}

\examples{
## EXAMPLE 1 (from Machin, Campbell, Tan and Tan 2009 p. 113):
## Bennett, Dismukes, Duma et al. (1979) designed a clinical trial to test 
## whether combination chemotherapy for a shorter period would be at least 
## as good as conventional therapy for patients with cryptococcal meningitis. 
## They recruited 39 patients to each treatment arm and wished to conclude 
## that a difference of less than 20\% in response rate between the treatments
## would indicate equivalence. Assuming a one-sided test size of 10\% and a 
## power of 80\% what would be a realistic sample size if the trial were 
## to be repeated?

epi.ssequb(treat = 0.50, control = 0.50, delta = 0.20, n = NA, power = 0.80, 
   r = 1, type = "equivalence", nfractional = FALSE, alpha = 0.10)

## A total of 166 subjects need to be enrolled in the trial, 83 in the 
## treatment group and 83 in the control group.
}

\keyword{univar}
