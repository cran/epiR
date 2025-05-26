## ----echo = FALSE, message = FALSE--------------------------------------------
library(dplyr); library(flextable); library(knitr); library(officer);  library(tidyr)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----echo = FALSE, results = 'asis'-------------------------------------------
twobytwo.df <- data.frame("exp" = c("Exposure+","Exposure-","Total"), "dpos" = c("a","c","a + c"), "dneg" = c("b","c","b + c"), "total" = c("a + b","c + d","a + b + c + d"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total"),
  h1 = c("", "Outcome+", "Outcome-", "Total"), stringsAsFactors = FALSE)

# Create table:
border_h = fp_border(color = "black", width = 2)

flextable(twobytwo.df) %>%
  width(j = 1, width = 2.00) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.00) %>%
  width(j = 4, width = 4.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%

  set_caption("A 2 x 2 contingency table.")

## ----echo = FALSE, results = 'asis'-------------------------------------------
irr.df <- data.frame("exp" = c("Exposure+","Exposure-","Total"), "dpos" = c("a","c","a + c"), "dneg" = c("b","c","b + c"), "total" = c("a + b","c + d","a + b + c + d"), risk = c("RE+ = a \U00F7 (a + b)", "RE- = c \U00F7 (c + d)", "RT = (a + c) \U00F7 (a + b + c + d)"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total","risk"),
  h1 = c("", "Outcome+", "Outcome-", "Total", "Risk"), stringsAsFactors = FALSE)

# Create table:
border_h = fp_border(color = "black", width = 2)

flextable(irr.df) %>%
  width(j = 1, width = 2.00) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.00) %>%
  width(j = 4, width = 2.00) %>%
  width(j = 5, width = 4.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%

  set_caption("A 2 x 2 table with incidence risks calculated for the exposed, the unexposed and the total study population.")


## ----risk_ratio, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "The incidence risk ratio."----
knitr::include_graphics("risk_ratio.png")

## ----echo = FALSE, results = 'asis'-------------------------------------------
orcohort.df <- data.frame("exp" = c("Exposure+","Exposure-","Total"), "dpos" = c("a","c","a + c"), "dneg" = c("b","d","b + d"), "total" = c("a + b","c + d","a + b + c + d"), odds = c("OE+ = a \U00F7 b","OE- = c \U00F7 d", "OT = (a + c) \U00F7 (b + d)"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total","odds"),
  h1 = c("", "Outcome+", "Outcome-", "Total", "Odds"), stringsAsFactors = FALSE)

# Create table:
border_h = fp_border(color = "black", width = 2)

flextable(orcohort.df) %>%
  width(j = 1, width = 2.00) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.00) %>%
  width(j = 4, width = 2.00) %>%
  width(j = 5, width = 4.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%

  set_caption("A 2 x 2 table with the outcome odds calculated for the exposed, unexposed and the total study population.")


## ----echo = FALSE, results = 'asis'-------------------------------------------
orcc.df <- data.frame("exp" = c("Outcome+ (case)","Outcome- (control)","Total"), "dpos" = c("a","c","a + c"), "dneg" = c("b","d","b + d"), "total" = c("a + b","c + d","a + b + c + d"), odds = c("OD+ = a \U00F7 b","OD- = c \U00F7 d", "OT = (a + c) \U00F7 (b + d)"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total","odds"),
  h1 = c("", "Exposure+", "Exposure-", "Total", "Odds"), stringsAsFactors = FALSE)

# Create table:
border_h = fp_border(color = "black", width = 2)

flextable(orcc.df) %>%
  width(j = 1, width = 2.00) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.00) %>%
  width(j = 4, width = 2.00) %>%
  width(j = 5, width = 4.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%

  set_caption("A 2 x 2 table with the exposure odds calculated for cases, controls and the total study population.")


## ----attributable_risk, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "The attributable risk in the exposed."----
knitr::include_graphics("attributable_risk.png")

## ----attributable_fraction, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "The attributable fraction in the exposed."----
knitr::include_graphics("attributable_fraction.png")

## ----population_attributable_risk, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "The attributable risk in the population."----
knitr::include_graphics("population_attributable_risk.png")

## ----population_attributable_fraction, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "The attributable fraction in the population."----
knitr::include_graphics("population_attributable_fraction.png")

## -----------------------------------------------------------------------------
dat.v01 <- c(13,2163,5,3349); dat.v01

# View the data in the usual 2 by 2 table format:
matrix(dat.v01, nrow = 2, byrow = TRUE)

## ----message = FALSE----------------------------------------------------------
library(epiR)

epi.2by2(dat = dat.v01, method = "cross.sectional", elab = "Dry food", olab = "FLUTD", digits = 2, conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")

## ----message = FALSE----------------------------------------------------------
epi.2by2(dat = dat.v01, method = "cross.sectional", elab = "Dry food", olab = "FLUTD", digits = 2, conf.level = 0.95, units = 100, interpret = TRUE, outcome = "as.columns")

## ----message = FALSE----------------------------------------------------------
library(MASS)

# Load and view the data:
dat.df02 <- birthwt; head(dat.df02)

## ----message = FALSE----------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02

## ----message = FALSE----------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02
dat.tab02 <- dat.tab02[2:1,2:1]; dat.tab02

## ----message = FALSE----------------------------------------------------------
# Variables low, smoke and race as factors. Put an 'f' in front of the variable names to remind you that they're factors:

dat.df02$flow <- factor(dat.df02$low, levels = c(1,0))
dat.df02$fsmoke <- factor(dat.df02$smoke, levels = c(1,0))
dat.df02$frace <- factor(dat.df02$race, levels = c(1,2,3))

dat.tab02 <- table(dat.df02$fsmoke, dat.df02$flow, dnn = c("Smoke", "Low BW")); dat.tab02

## ----message = FALSE----------------------------------------------------------
dat.epi02 <- epi.2by2(dat = dat.tab02, method = "cohort.count", elab = "Smoke", olab = "Low BW", digits = 2, conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi02

## ----message = FALSE----------------------------------------------------------
names(dat.epi02$massoc.detail)

## ----message = FALSE----------------------------------------------------------
dat.epi02$massoc.detail$OR.strata.wald
# Wald confidence intervals: 2.0 (95% CI 1.1 to 3.8)

dat.epi02$massoc.detail$OR.strata.score
# Score confidence intervals: 2.0 (95% CI 1.1 to 3.8)


## ----message = FALSE----------------------------------------------------------
library(dplyr); library(tidyr)

dat.df03 <- birthwt; head(dat.df03)

# Here we set the factor levels and tabulate the data in a single call using pipe operators:
dat.tab03 <- dat.df03 %>%
  mutate(flow = factor(low, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(fsmoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) %>%
  group_by(fsmoke, flow) %>%
  summarise(n = n()) 

# View the data:
dat.tab03

# View the data in conventional 2 by 2 table format:
pivot_wider(dat.tab03, id_cols = c(fsmoke), 
   names_from = flow, values_from = n)

## ----message = FALSE----------------------------------------------------------
dat.epi03 <- epi.2by2(dat = dat.tab03, method = "cohort.count", elab = "Smoke", olab = "Low BW", digits = 2, conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi03

## ----message = FALSE----------------------------------------------------------
dat.df04 <- birthwt; head(dat.df04)

dat.df04$flow <- factor(dat.df04$low, levels = c(1,0))

dat.df04$fage <- ifelse(dat.df04$age > 23, 0,1)
dat.df04$fage <- factor(dat.df04$fage, levels = c(1,0))
dat.df04$fsmoke <- factor(dat.df04$smoke, levels = c(1,0))

# Race is coded 1 = white, 2 = black and 3 = other. Set white as the reference level:
dat.df04$frace <- ifelse(dat.df04$race == 1, 0, 1)
dat.df04$frace <- factor(dat.df04$frace, levels = c(1,0))

# Empty vectors to collect results:
rfactor <- ref <- or.p <- or.l <- or.u <- c() 

# The candidate risk factors are in columns 12 to 14 of data frame dat.df04:
for(i in 12:14){
  tdat.tab04 <- table(dat.df04[,i], dat.df04$flow)
  tdat.epi04 <- epi.2by2(dat = tdat.tab04, method = "cohort.count", 
   digits = 2, conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
  
  trfactor <- as.character(names(dat.df04)[i]) 
  rfactor <- c(rfactor, trfactor) 
  
  tref <- as.character(paste("Reference: ", trfactor, " - ", levels(dat.df04[,i])[2], sep = ""))
  ref <- c(ref, tref)
  
  tor.p <- as.numeric(tdat.epi04$massoc.detail$OR.strata.wald[1])
  or.p <- c(or.p, tor.p)
  
  tor.l <- as.numeric(tdat.epi04$massoc.detail$OR.strata.wald[2])
  or.l <- c(or.l, tor.l)
  
  tor.u <- as.numeric(tdat.epi04$massoc.detail$OR.strata.wald[3])
  or.u <- c(or.u, tor.u)
}

gdat.df04 <- data.frame(ybrk = 1:3, ylab = rfactor, ref, or.p, or.l, or.u)
gdat.df04

## ----odds_ratios, echo = TRUE, message = FALSE, fig.align = "center", out.width = "80%", fig.cap = "Risk factors for low birth weight babies. Error bar plot showing the point estimate of the odds ratio and its 95% confidence interval for maternal age, smoking and race."----
library(ggplot2); library(scales)

xbrk <- seq(from = -2, to = 2, by = 1)
xlab <- 2^xbrk

ggplot(data = gdat.df04, aes(x = log2(or.p), y = ybrk)) +
  theme_bw() +
  geom_point() + 
  geom_errorbarh(aes(xmin = log2(or.l), xmax = log2(or.u), height = 0.2)) + 
  scale_x_continuous(breaks = xbrk, labels = xlab, limits = range(xbrk), 
   name = "Odds ratio") + 
  scale_y_continuous(breaks = gdat.df04$ybrk, labels = gdat.df04$ylab, name = "Risk factor") + 
  geom_vline(xintercept = log2(1), linetype = "dashed") + 
  annotate("text", x = log2(0.25), y = gdat.df04$ybrk, label = gdat.df04$ref, hjust = 0, size = 3) +
  coord_fixed(ratio = 0.75 / 1) + 
  theme(axis.title.y = element_text(vjust = 0))

## ----message = FALSE----------------------------------------------------------
dat.df05 <- birthwt; head(dat.df05)

dat.df05$flow <- factor(dat.df05$low, levels = c(1,0))

dat.df05$fsmoke <- factor(dat.df05$smoke, levels = c(1,0))
dat.df05$frace <- factor(dat.df05$race, levels = c(1,2,3))

dat.tab05 <- table(dat.df05$fsmoke, dat.df05$flow, dat.df05$frace, 
   dnn = c("Smoke", "Low BW", "Race")); dat.tab05

## ----message = FALSE----------------------------------------------------------
dat.epi05 <- epi.2by2(dat = dat.tab05, method = "cohort.count", elab = "smoke",
   olab = "Low BW", digits = 2, conf.level = 0.95, units = 100, 
   interpret = FALSE, outcome = "as.columns")
dat.epi05

## ----message = FALSE----------------------------------------------------------
dat.df06 <- birthwt

dat.tab06 <- dat.df06 %>%
  mutate(flow = factor(low, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(fsmoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(frace = factor(race)) %>%
  group_by(frace, fsmoke, flow) %>%
  summarise(n = n()) 
dat.tab06

# View the data in conventional 2 by 2 table format:
pivot_wider(dat.tab06, id_cols = c(frace, fsmoke), 
   names_from = flow, values_from = n)

## ----message = FALSE----------------------------------------------------------
dat.epi06 <- epi.2by2(dat = dat.tab06, method = "cohort.count", 
   elab = "Smoke", olab = "Low BW", digits = 2, conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")

dat.epi06

## ----mantel_haenszel, echo = TRUE, message = FALSE, fig.align = "center", out.width = "80%", fig.cap = "Risk factors for low birth weight babies. Error bar plot showing the odds of having a low birth weight baby for smokers of maternal race categories 1, 2 and 3 and the Mantel-Haenszel odds of having a low birth weight baby for smokers, adjusted for maternal race."----

xbrk <- seq(from = -5, to = 5, by = 1)
xlab <- 2^xbrk

nstrata <- dat.epi06$n.strata
ybrk <- c(1:nstrata, max(nstrata) + 1)
ylab <- c("M-H", paste("Strata ", 1:nstrata, sep = ""))

or.p <- c(dat.epi06$massoc.detail$OR.mh$est, 
   dat.epi06$massoc.detail$OR.strata.cfield$est)
or.l <- c(dat.epi06$massoc.detail$OR.mh$lower, 
   dat.epi06$massoc.detail$OR.strata.cfield$lower)
or.u <- c(dat.epi06$massoc.detail$OR.mh$upper, 
   dat.epi06$massoc.detail$OR.strata.cfield$upper)
gdat.df06 <- data.frame(ybrk, ylab, or.p, or.l, or.u)

ggplot(data = gdat.df06, aes(x = log2(or.p), y = ybrk)) +
  theme_bw() +
  geom_point() + 
  geom_errorbarh(aes(xmin = log2(or.u), xmax = log2(or.l), height = 0.2)) + 
  scale_x_continuous(breaks = xbrk, labels = xlab, limits = range(xbrk), 
   name = "Odds ratio") + 
  scale_y_continuous(breaks = gdat.df06$ybrk, labels = gdat.df06$ylab, name = "Risk factor") + 
  geom_vline(xintercept = log2(1), linetype = "dashed") + 
  annotate("text", x = log2(0.03125), y = gdat.df06$yat, label = gdat.df06$ref, hjust = 0, size = 3) +
  coord_fixed(ratio = 0.75 / 1) + 
  theme(axis.title.y = element_text(vjust = 0))


