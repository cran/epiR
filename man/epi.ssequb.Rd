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
   nfractional = FALSE, alpha = 0.05)
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
Bennett J, Dismukes W, Duma R, Medoff G, Sande M, Gallis H, Leonard J, Fields B, Bradshaw M, Haywood H, McGee Z, Cate T, Cobbs C, Warner J, Alling D (1979). A comparison of amphotericin B alone and combined with flucytosine in the treatment of cryptoccal meningitis. New England Journal of Medicine 301, 126 - 131. DOI: 10.1056/nejm197907193010303.

Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, pp. 91.

Ewald B (2013). Making sense of equivalence and non-inferiority trials. Australian Prescriber 36: 170 - 173.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Julious SA (2009). Estimating Samples Sizes in Clinical Trials. CRC, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2009). Sample Size Tables for Clinical Studies. Wiley Blackwell, New York.

Machin D, Campbell MJ, Tan SB, Tan SH (2018). Sample Sizes for Clinical, Laboratory and Epidemiology Studies. Wiley Blackwell, New York.

Wang B, Wang H, Tu X, Feng C (2017). Comparisons of superiority, non-inferiority, and equivalence trials. Shanghai Archives of Psychiatry 29, 385 - 388. DOI: 10.11919/j.issn.1002-0829.217163.
}

\note{
Consider a clinical trial comparing two groups, a standard treatment (\eqn{s}) and a new treatment (\eqn{n}). A proportion of subjects in the standard treatment group experience the outcome of interest \eqn{P_{s}} and a proportion of subjects in the new treatment group experience the outcome of interest \eqn{P_{n}}.

For an equivalence trial the null hypothesis is:

\eqn{H_{0}: P_{n} - P_{s} \leq -\delta}

\eqn{H_{0}: P_{n} - P_{s} \geq +\delta}
                                 
and the alternative hypothesis is:

\eqn{H_{1}: P_{n} - P_{s} > -\delta}

\eqn{H_{1}: P_{n} - P_{s} < +\delta}

If the alternative hypothesis is accepted we conclude that \eqn{P_{n}} is therapeutically equivalent to \eqn{P_{s}}. For an equivalence trial a two-sided test is used.

For an equality trial the null hypothesis is:

\eqn{H_{0}: P_{n} - P_{s} = 0}

and the alternative hypothesis is:

\eqn{H_{1}: P_{n} - P_{s} \neq 0}

If the alternative hypothesis is accepted we conclude that \eqn{P_{n}} differs from \eqn{P_{s}}. For an equality trial a two-sided test is used.

Note that when:

\deqn{sign(P_n - P_s - \delta) \neq sign(z_{1-\alpha} + z_{1-\beta})}

there is no solution for study power. For typical values of \eqn{\alpha} and \eqn{\beta} this would occur if \eqn{P_{n} - P_{s} - \delta < 0}. That is, when the targeted alternative is within the null hypothesis. The function issues a warning if these conditions are met.

When calculating study power the argument \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).

A summary of the key features of superiority, equivalence and non-superiority trial comparisons is shown in the following table (adapted from Machin et al., 2018 [page 170] and Wang et al., 2017):

\tabular{llll}{
-------------------------------- \tab -------------------------------- \tab -------------------------------- \tab -------------------------------- \cr
Test for                         \tab Null hypothesis \tab Alternative hypothesis \tab Test \cr
-------------------------------- \tab -------------------------------- \tab -------------------------------- \tab -------------------------------- \cr
Superiority                      \tab \eqn{H_{0}: P_{n} - P_{s} \leq \delta} \tab \eqn{H_{1}: P_{n} - P_{s} > \delta} \tab 1-sided or 2-sided \cr
Equivalence                      \tab \eqn{H_{0}: P_{n} - P_{s} \leq -\delta} \tab \eqn{H_{1}: P_{n} - P_{s} > -\delta} \tab Two 1-sided tests \cr
                                 \tab \eqn{H_{0}: P_{n} - P_{s} \geq +\delta} \tab \eqn{H_{1}: P_{n} - P_{s} < +\delta} \tab \cr
Equality                         \tab \eqn{H_{0}: P_{n} - P_{s} == 0} \tab \eqn{H_{1}: P_{n} - P_{s} \neq 0} \tab 2-sided \cr
Non-inferiority                  \tab \eqn{H_{0}: P_{n} - P_{s} \leq -\delta} \tab \eqn{H_{1}: P_{n} - P_{s} > -\delta} \tab 1-sided \cr
-------------------------------- \tab -------------------------------- \tab -------------------------------- \tab -------------------------------- \cr
}
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
