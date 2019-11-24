\name{epi.ssequb}

\alias{epi.ssequb}

\title{
Sample size for a parallel equivalence trial, binary outcome
}

\description{
Sample size for a parallel equivalence trial, binary outcome. 
}

\usage{
epi.ssequb(treat, control, delta, n, r = 1, power, alpha)
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
A list containing the following: 
  \item{n.total}{the total number of study subjects required.}
  \item{n.treat}{the required number of study subject in the treatment group.}
  \item{n.control}{the required number of study subject in the control group.}
  \item{power}{the specified or calculated study power.}
}

\references{
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

To summarise (adapted from Machin et al. 2009, page 105):

\tabular{lllll}{
--------------- \tab ------------------------ \tab ----------------------- \tab --------------------- \tab --------------------- \cr
Test for        \tab Null hypothesis          \tab Alt hypothesis          \tab Type I                \tab Type II               \cr
--------------- \tab ------------------------ \tab ----------------------- \tab --------------------- \tab --------------------- \cr
Superiority     \tab H0:  Pn - Ps  <= delta   \tab H1:  Pn - Ps  > delta   \tab 2 sided, 5.0\%        \tab 1 sided, 10 or 20\%   \cr
Equivalence     \tab H0: |Pn - Ps| >= delta   \tab H1: |Pn - Ps| < delta   \tab 1 sided, 5.0\%        \tab 1 sided, 20\%         \cr
Non-inferiority \tab H0:  Pn - Ps  >= delta   \tab H1:  Pn - Ps  < delta   \tab 1 sided, 2.5\%        \tab 1 sided, 10 or 20\%   \cr
--------------- \tab ------------------------ \tab ----------------------- \tab --------------------- \tab --------------------- \cr
}

Superiority trial: H1 is that the new treatment is better than the standard treatment.\cr
Equivalence trial: H1 is that the new treatment is not too different from the standard treatment.\cr
Non-inferiority trial: 	H1 is that the new treatment is not much worse than the standard treatment.

When calculating the power of a study, note that the variable \code{n} refers to the total study size (that is, the number of subjects in the treatment group plus the number in the control group).
}

\examples{
## EXAMPLE 1 (from Machin, Campbell, Tan and Tan 2009 p. 113):
## Bennett, Dismukes, Duma et al. (1979) designed a clinical trial to test 
## whether combination chemotherapy for a shorter period would be at least 
## as good as conventional therapy for patients with cryptococcal meningitis. 
## They recruited 39 patients to each treatment arm and wished to conclude 
## that a difference of less than 20\% in response rate between the treatments
## would indicate equivalence. Assuming a one-sided test size of 10\%, a 
## power of 80\% and an overall response rate of 50\%, what would be a 
## realistic sample size if the trial were to be repeated?

epi.ssequb(treat = 0.50, control = 0.50, delta = 0.20, n = NA, r = 1, 
   power = 0.80, alpha = 0.10)

## A total of 166 subjects need to be enrolled in the trial, 83 in the 
## treatment group and 83 in the control group.
}

\keyword{univar}
