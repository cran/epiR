\name{epi.realrisk}

\alias{epi.realrisk}

\title{An R version of the Winton Centre's RealRisk calculator}

\description{
A function to return absolute risks from odds, incidence risk and hazard ratios.
}

\usage{
epi.realrisk(measure, measure.type, rE0)
}

\arguments{
  \item{measure}{scalar, the numeric estimate of the specified measure of association.}
  \item{measure.type}{character string, detailing the measure of association type. Options are \code{odds.ratio}, \code{risk.ratio}, \code{hazard.ratio} and \code{pct.change} (percentage change).}
  \item{rE0}{scalar, the prevalence of the outcome in the exposure negative ('baseline') group. Must be a number between 0 and 1.}
}

\value{
A data frame with four elements: \code{measure}, the numeric estimate of the measure of association (as entered by the user), \code{measure.type} the measure of association type (as entered by the user), \code{rE0} the prevalence of the outcome in the exposure negative ('baseline') group and \code{rE1} the prevalence of the outcome in the exposure positive ('comparison') group.}

\details{
For numerous worked examples and analytical details see \url{https://realrisk.wintoncentre.uk/}.    
}

\references{
Horby P, Lim W, Emberson J, Mafham M, Bell J, Linsell L, Staplin N, Brightling C, Ustianowski A, Elmahi E, Prudon B, Green C, Felton T, Chadwick D, Rege K, Fegan C, Chappell L, Faust S, Jaki T, Jeffery K, Montgomery A, Rowan K, Juszczak E, Baillie J, Haynes R, Landray M (2021). Dexamethasone in Hospitalized Patients with COVID-19. New England Journal of Medicine 384, 693 - 704.

Parekh T, Pemmasani S, Desai R (2020). Risk of stroke with E-cigarette and combustible cigarette use in young adults. American Journal of Preventive Medicine 58, 446 - 452.

The Winton Centre for Risk and Evidence Communication, The University of Cambridge. URL: \url{https://realrisk.wintoncentre.uk/}.
 
}

\examples{
## EXAMPLE 1:
## A cross-sectional study was conducted to determine if 18 to 44 year olds 
## who routinely used e-cigarettes were at greater risk of stroke compared 
## with individuals of the same age who did not smoke (Parekh et al., 2021). 
## The risk of stroke among those 18 to 44 years of age was estimated to 
## be 0.6\% per year. 

## Current dual use of e-cigarettes and combustible cigarettes was associated 
## with a 2.91 times higher odds of stroke compared with non-smokers 
## (adjusted OR = 2.91, 95\% CI 1.62 to 5.25). Report these findings in 
## absolute terms.

epi.realrisk(measure = 2.91, measure.type = "odds.ratio", rE0 = 0.006)

## For those who do not smoke we would expect around 1 in every 100 to suffer
## from stroke. For those who smoke we would expect around 2 in every 100
## to suffer from stroke.


## EXAMPLE 2:
## In a controlled, open-label trial comparing a range of possible treatments 
## in patients who were hospitalised with COVID-19 (Horby et al., 2021), 
## patients were randomly assigned to receive oral or intravenous 
## dexamethasone for up to 10 days or to receive usual care alone. The 
## primary outcome for this study was 28-day mortality. The risk of 28-day 
## mortality among patients receiving mechanical ventilation and usual care 
## was 0.414.

## The risk of death within 28 days of admission among mechanically 
## ventilated COVID-19 patients receiving dexamethasone was 0.64 
## (95\% CI 0.51 to 0.82) times that of mechanically ventilated COVID-19 
## patients receiving usual care. Report these findings in absolute terms.

epi.realrisk(measure = c(0.51,0.64,0.82), measure.type = "risk.ratio", 
   rE0 = 0.414)

## Out of 100 UK patients receiving mechanical ventilation and 
## usual care for COVID-19, we would expect around 41 to die after 28 days. 
## Out of 100 UK patients receiving mechanical ventilation and dexamethasone 
## for COVID-19, we would expect around 26 (95\% CI 21 to 34) to die after 
## 28 days.
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
