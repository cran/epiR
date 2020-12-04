## ---- echo = FALSE, message = FALSE-------------------------------------------

# If you want to create a PDF document paste the following after line 9 above:
#   pdf_document:
#     toc: true
#     highlight: tango
#     number_sections: no
#     latex_engine: xelatex    
# header-includes: 
#    - \usepackage{fontspec}

knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----ssrs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------
library(pander)
panderOptions('table.split.table', Inf)

set.caption("Functions to estimate sample size using representative population sampling data.")

ssrs.tab <- " 
Sampling        | Outcome               | RSurveillance              | epiR
Representative  | Prob disease freedom  | `n.pfree`                  | `rsu.sspfree.rs`
Representative  | SSe                   | `n.freedom`                | `rsu.sssep.rs`
Two stage representative | SSe         | `n.2stage`                 | `rsu.sssep.rs2st`
Representative  | SSe                   | `n.freecalc`               | `rsu.sssep.rsfreecalc`
Pooled representative    | SSe                   | `n.pooled`                 | `rsu.sssep.rspool`"

ssrs.df <- read.delim(textConnection(ssrs.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(ssrs.df) <- unname(as.list(ssrs.df[1,])) # put headers on
ssrs.df <- ssrs.df[-1,] # remove first row
row.names(ssrs.df) <- NULL
pander(ssrs.df, style = 'rmarkdown')

## ----seprs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----

set.caption("Functions to estimate surveillance system sensitivity (SSe) using representative population sampling data.")

seprs.tab <- " 
Sampling                  | Outcome      | RSurveillance        | epiR
Representative            | SSe          | `sep.binom`          | `rsu.sep.rs`
Representative            | SSe          | `sep.hypergeo`       | `rsu.sep.rs`
Representative            | SSe          | `sep`                | `rsu.sep.rs`
Two stage representative  | SSe          | `sep.sys`            | `rsu.sep.rs2st`
Representative            | SSe          | `sse.combined`       | `rsu.sep.rsmult`
Representative            | SSe          | `sep.freecalc`       | `rsu.sep.rsfreecalc`
Representative            | SSe          | `sep.binom.imperfect`| `rsu.sep.rsfreecalc`
Pooled representative     | SSe          | `sep.pooled`         | `rsu.sep.rspool`
Representative            | SSe          | `sep.var.se`         | `rsu.sep.rsvarse`
Representative            | SSp          | `spp`                | `rsu.spp.rs`
Representative            | SSp          | `sph.hp`             | `rsu.spp.rs`"

seprs.df <- read.delim(textConnection(seprs.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(seprs.df) <- unname(as.list(seprs.df[1,])) # put headers on
seprs.df <- seprs.df[-1,] # remove first row
row.names(seprs.df) <- NULL
pander(seprs.df, style = 'rmarkdown')

## ----pfreers.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----

set.caption("Functions to estimate the probability of disease freedom using representative population sampling data.")

pfreers.tab <- " 
Sampling          | Outcome                             |  RSurveillance      | epiR
Representative    | Prob disease of freedom             | `pfree.1`           | `rsu.pfree.rs`
Representative    | Prob disease of freedom             | `pfree.calc`        | `rsu.pfree.rs`
Representative    | Equilibrium prob of disease freedom | `pfree.equ`         | `rsu.pfree.equ`"

pfreers.df <- read.delim(textConnection(pfreers.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(pfreers.df) <- unname(as.list(pfreers.df[1,])) # put headers on
pfreers.df <- pfreers.df[-1,] # remove first row
row.names(pfreers.df) <- NULL
pander(pfreers.df, style = 'rmarkdown')

## ----ssrb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------

set.caption("Functions to estimate sample size using risk based sampling data.")

ssrb.tab <- " 
Sampling       | Outcome | RSurveillance        | epiR
Risk-based     | SSe     | `n.rb`               | `rsu.sssep.rbsrg`
Risk-based     | SSe     | `n.rb.varse`         | `rsu.sssep.rbmrg`
Risk-based     | SSe     | `n.rb.2stage.1`      | `rsu.sssep.rb2st1rf`
Risk-based     | SSe     | `n.rb.2stage.2`      | `rsu.sssep.rb2st2rf`"

ssrb.df <- read.delim(textConnection(ssrb.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(ssrb.df) <- unname(as.list(ssrb.df[1,])) # put headers on
ssrb.df <- ssrb.df[-1,] # remove first row
row.names(ssrb.df) <- NULL
pander(ssrb.df, style = 'rmarkdown')

## ----seprb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----

set.caption("Functions to estimate surveillance system sensitivity (SSe) using risk based sampling data.")

seprb.tab <- " 
Sampling       | Outcome | RSurveillance           | epiR
Risk-based     | SSe     | `sep.rb.bin.varse`      | `rsu.sep.rb`
Risk-based     | SSe     | `sep.rb.bin`            | `rsu.sep.rb1rf`
Risk-based     | SSe     | `sep.rb.hypergeo`       | `rsu.sep.rb1rf`
Risk-based     | SSe     | `sep.rb2.bin`           | `rsu.sep.rb2rf`
Risk-based     | SSe     | `sep.rb2.hypergeo`      | `rsu.sep.rb2rf`
Risk-based     | SSe     | `sep.rb.hypergeo.varse` | `rsu.sep.rbvarse`
Risk-based     | SSe     | `sse.rb2stage`          | `rsu.sep.rb2stage`"
seprb.df <- read.delim(textConnection(seprb.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(seprb.df) <- unname(as.list(seprb.df[1,])) # put headers on
seprb.df <- seprb.df[-1,] # remove first row
row.names(seprb.df) <- NULL
pander(seprb.df, style = 'rmarkdown')

## ----sepcen.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----

set.caption("Functions to estimate surveillance system sensitivity (SSe) using census data.")

sepcen.tab <- " 
Sampling       | Outcome | RSurveillance           | epiR
Risk-based     | SSe     | `sep.exact`             | `rsu.sep.cens`"
sepcen.df <- read.delim(textConnection(sepcen.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(sepcen.df) <- unname(as.list(sepcen.df[1,])) # put headers on
sepcen.df <- sepcen.df[-1,] # remove first row
row.names(sepcen.df) <- NULL
pander(sepcen.df, style = 'rmarkdown')

## ----seppas.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----

set.caption("Functions to estimate surveillance system sensitivity (SSe) using passively collected surveillance data.")

sepcen.tab <- " 
Sampling       | Outcome | RSurveillance           | epiR
Risk-based     | SSe     | `sep.passive`           | `rsu.sep.pass`"
seppas.df <- read.delim(textConnection(sepcen.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(seppas.df) <- unname(as.list(seppas.df[1,])) # put headers on
seppas.df <- seppas.df[-1,] # remove first row
row.names(seppas.df) <- NULL
pander(seppas.df, style = 'rmarkdown')

## ----misc.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------

set.caption("Miscellaneous functions.")

misc.tab <- " 
Details                             | RSurveillance           | epiR
Adjusted risk                       | `adj.risk`              | `rsu.adjrisk`
Adjusted risk                       | `adj.risk.sim`          | `rsu.adjrisk`
Series test interpretation, Se      | `se.series`             | `rsu.dxtest`
Parallel test interpretation, Se    | `se.parallel`           | `rsu.dxtest`

Series test interpretation, Sp      | `sp.series`             | `rsu.dxtest`
Parallel test interpretation, Sp    | `sp.parallel`           | `rsu.dxtest`

Effective probability of infection  | `epi.calc`              | `rsu.epinf`
Design prevalence back calculation  | `pstar.calc`            | `rsu.pstar`
Prob disease is less than design prevalence |                 | `rsu.sep`"
misc.df <- read.delim(textConnection(misc.tab), header = FALSE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)

names(misc.df) <- unname(as.list(misc.df[1,])) # put headers on
misc.df <- misc.df[-1,] # remove first row
row.names(misc.df) <- NULL
pander(misc.df, style = 'rmarkdown')

