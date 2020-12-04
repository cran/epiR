\name{rsu.sssep.rsfreecalc}

\alias{rsu.sssep.rsfreecalc}

\title{
Sample size to achieve a desired surveillance system sensitivity to detect disease at a specified design prevalence assuming representative sampling, imperfect unit sensitivity and specificity
}

\description{
Calculates the sample size to achieve a desired surveillance system sensitivity to detect disease at a specified design prevalence assuming representative sampling, imperfect unit sensitivity and specificity .}

\usage{
rsu.sssep.rsfreecalc(N, pstar, mse.p, msp.p, se.u, sp.u, method = "hypergeometric", 
   max.ss = 32000)
}

\arguments{
  \item{N}{scalar, integer representing the total number of subjects eligible to be sampled.}
  \item{pstar}{scalar, numeric, representing the design prevalence, the hypothetical outcome prevalence to be detected. See details, below.}
  \item{mse.p}{scalar, numeric (0 to 1) representing the desired population level sensitivity. See details, below.}
  \item{msp.p}{scalar, numeric (0 to 1) representing the desired population level specificity. See details, below.}
  \item{se.u}{scalar (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}
  \item{sp.u}{scalar, numeric (0 to 1) representing the specificity of the diagnostic test at the surveillance unit level.}
  \item{method}{a character string indicating the calculation method to use. Options are \code{binomial} or \code{hypergeometric}.}
  \item{max.ss}{scalar, integer defining the maximum upper limit for required sample size.}
}

\details{
Type I error is the probabilty of rejecting the null hypothesis when in reality it is true. In disease freedom studies this is the situation where you declare a population as disease negative when, in fact, it is actually disease positive. Type I error equals \code{1 - SeP}. 

Type II error is the probabilty of accepting the null hypothesis when in reality it is false. In disease freedom studies this is the situation where you declare a population as disease positive when, in fact, it is actually disease negative. Type II error equals \code{1 - SpP}.

Argument \code{pstar} can be expressed as either a proportion or integer. Where the input value for \code{pstar} is between 0 and 1 the function interprets \code{pstar} as a prevalence. Where the input value for \code{pstar} is an integer greater than 1 the function interprets \code{pstar} as the number of outcome-positive individuals in the population of individuals at risk. A value for design prevalence is then calculated as \code{pstar / N}.     
}

\value{
A list comprised of two data frames: \code{summary} and \code{details}. Data frame \code{summary} lists:
 
\item{n}{the minimum number of individuals to be sampled.}
\item{N}{the total number of individuals eligible to be sampled.}
\item{c}{the cut-point number of positives to achieve the specified surveillance system (population-level) sensitivity and specificity.}
\item{pstar}{the design prevalence.}
\item{p1}{the probability that the population has the outcome of interest at the specified design prevalence.}
\item{se.p}{the calculated population level sensitivity.}
\item{sp.p}{the calculated population level specificity.}
  
Data frame \code{details} lists: 

\item{n}{the minimum number of individuals to be sampled.}
\item{se.p}{the calculated population level sensitivity.}
\item{sp.p}{the calculated population level specificity.}
}

\references{
Cameron A, Baldock C (1998a). A new probability formula for surveys to substantiate freedom from disease. Preventive Veterinary Medicine 34: 1 - 17.

Cameron A, Baldock C (1998b). Two-stage sampling in surveys to substantiate freedom from disease. Preventive Veterinary Medicine 34: 19 - 30.

Cameron A (1999). Survey Toolbox for Livestock Diseases --- A practical manual and software package for active surveillance of livestock diseases in developing countries. Australian Centre for International Agricultural Research, Canberra, Australia.
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## brucellosis in dairy herds using a bulk milk tank test assuming a design 
## prevalence of 0.05. Assume the total number of dairy herds in your study 
## area is 5000 and the bulk milk tank test to be used has a diagnostic 
## sensitivity of 0.95 and a specificity of 1.00. How many herds need to be 
## sampled to be 95\% certain that the prevalence of brucellosis in dairy herds 
## is less than the design prevalence if less than a specified number of 
## tests return a positive result?

rsu.sssep.rsfreecalc(N = 5000, pstar = 0.05, mse.p = 0.95, msp.p = 0.95, 
   se.u = 0.95, sp.u = 0.98, method = "hypergeometric", max.ss = 32000)$summary

## A system sensitivity of 95\% is achieved with a total sample size of 194 
## herds, assuming a cut-point of 7 or more positive herds are required to 
## return a positive survey result.
}

\keyword{univar}

