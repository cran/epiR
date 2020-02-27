\name{epi.sssupb}

\alias{epi.sssupb}

\title{
Sample size for a parallel superiority trial, binary outcome
}

\description{
Sample size for a parallel superiority trial, binary outcome. 
}

\usage{
epi.sssupb(treat, control, delta, n, r = 1, power, alpha)
}

\arguments{
  \item{treat}{the expected proportion of successes in the treatment group.}
  \item{control}{the expected proportion of successes in the control group.}
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
Chow S, Shao J, Wang H (2008). Sample Size Calculations in Clinical Research. Chapman & Hall/CRC Biostatistics Series, page 90.

Julious SA (2004). Sample sizes for clinical trials with normal data. Statistics in Medicine 23: 1921 - 1986.

Pocock SJ (1983). Clinical Trials: A Practical Approach. Wiley, New York.
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
## EXAMPLE 1 (from Chow S, Shao J, Wang H 2008, p. 91):
## Suppose that a pharmaceutical company is interested in conducting a
## clinical trial to compare the efficacy of two antimicrobial agents 
## when administered orally once daily in the treatment of patients 
## with skin infections. In what follows, we consider the situation 
## where the intended trial is for testing superiority of the 
## test drug over the active control drug. For this purpose, the following 
## assumptions are made. First, sample size calculation will be performed 
## for achieving 80\% power at the 5\% level of significance.

## Assume the true mean cure rates of the treatment agents and the active 
## control are 85\% and 65\%, respectively. Assume the superiority
## margin is 5\%.

epi.sssupb(treat = 0.85, control = 0.65, delta = 0.05, n = NA, 
   r = 1, power = 0.80, alpha = 0.05)

## A total of 196 subjects need to be enrolled in the trial, 98 in the 
## treatment group and 98 in the control group.

}

\keyword{univar}
