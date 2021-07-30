## ---- echo = FALSE, message = FALSE-------------------------------------------
library(knitr); library(kableExtra)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----echo = FALSE, results = 'asis'-------------------------------------------
twobytwo <- data.frame("Dis pos" = c("a","c","a+c"), "Dis neg" = c("b","c","b+c"), "Total" = c("a+b","c+d","a+b+c+d"))
colnames(twobytwo) <- c("Dis pos","Dis pos","Total")
rownames(twobytwo) <- c("Exp pos","Exp neg","Total")

kbl(twobytwo, caption = "A 2 by 2 table.") %>%
   column_spec(1, bold = FALSE, width = "5em") %>%
   column_spec(2, bold = FALSE, width = "5em") %>%
   column_spec(3, bold = FALSE, width = "5em")
   # row_spec(row = 1, bold = TRUE)

## ----echo = FALSE, results = 'asis'-------------------------------------------
irr <- data.frame("Dis pos" = c("a","c","a+c"), "Dis neg" = c("b","c","b+c"), "Total" = c("a+b","c+d","a+b+c+d"),"Risk" = c("RE+ = a/(a+b)","RE- = c/(c+d)", "RT = (a+c)/(a+b+c+d)"))

colnames(irr) <- c("Dis pos","Dis pos","Total","Risk")
rownames(irr) <- c("Exp pos","Exp neg","Total")

kbl(irr, caption = "A 2 by 2 table with incidence risks calculated for the exposed, the unexposed and the total study population.") %>%
   column_spec(1, bold = FALSE, width = "5em") %>%
   column_spec(2, bold = FALSE, width = "5em") %>%
   column_spec(3, bold = FALSE, width = "5em") %>%
   column_spec(4, bold = FALSE, width = "10em")
   # row_spec(row = 1, bold = TRUE)

## ----echo = FALSE, results = 'asis'-------------------------------------------
or.cohort <- data.frame("Dis pos" = c("a","c","a+c"), "Dis neg" = c("b","d","b+d"), "Total" = c("a+b","c+d","a+b+c+d"),"Odds" = c("OE+ = a/b","OE- = c/d", "OT = (a+c)/(b+d)"))

colnames(or.cohort) <- c("Dis pos","Dis pos","Total","Odds")
rownames(or.cohort) <- c("Exp pos","Exp neg","Total")

kbl(or.cohort, caption = "A 2 by 2 table with the odds of disease calculated for the exposed, the unexposed and the total study population.") %>%
   column_spec(1, bold = FALSE, width = "5em") %>%
   column_spec(2, bold = FALSE, width = "5em") %>%
   column_spec(3, bold = FALSE, width = "5em") %>%
   column_spec(4, bold = FALSE, width = "10em")
   # row_spec(row = 1, bold = TRUE)

## ----echo = FALSE, results = 'asis'-------------------------------------------
or.cc <- data.frame("Case" = c("a","c","a+c","OD+ = a/c"), "Control" = c("b","d","b+d","OD- = b/d"), "Total" = c("a+b","c+d","a+b+c+d","OT = (a+b)/(c+d)"))

colnames(or.cc) <- c("Case","Control","Total")
rownames(or.cc) <- c("Exp pos","Exp neg","Total","Odds")

kbl(or.cc, caption = "A 2 by 2 table with the odds of exposure calculated for cases, controls and the total study population.") %>%
   column_spec(1, bold = FALSE, width = "5em") %>%
   column_spec(2, bold = FALSE, width = "5em") %>%
   column_spec(3, bold = FALSE, width = "5em") %>%
   column_spec(4, bold = FALSE, width = "10em")
   # row_spec(row = 1, bold = TRUE)

## -----------------------------------------------------------------------------
dat.v01 <- c(13,2163,5,3349); dat.v01

# View the data in the usual 2 by 2 table format:
matrix(dat.v01, nrow = 2, byrow = TRUE)

## -----------------------------------------------------------------------------
library(epiR)

epi.2by2(dat = dat.v01, method = "cross.sectional", conf.level = 0.95, units = 100, 
   interpret = FALSE, outcome = "as.columns")

## -----------------------------------------------------------------------------
epi.2by2(dat = dat.v01, method = "cross.sectional", conf.level = 0.95, units = 100, 
   interpret = TRUE, outcome = "as.columns")

## -----------------------------------------------------------------------------
library(MASS)

# Load and view the data:
dat.df02 <- birthwt; head(dat.df02)

## -----------------------------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02

## -----------------------------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02
dat.tab02 <- dat.tab02[2:1,2:1]; dat.tab02

## -----------------------------------------------------------------------------
dat.df02$low <- factor(dat.df02$low, levels = c(1,0))
dat.df02$smoke <- factor(dat.df02$smoke, levels = c(1,0))
dat.df02$race <- factor(dat.df02$race, levels = c(1,2,3))

dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02

## -----------------------------------------------------------------------------
dat.epi02 <- epi.2by2(dat = dat.tab02, method = "cohort.count", conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi02

## -----------------------------------------------------------------------------
library(tidyverse)

dat.df03 <- birthwt; head(dat.df03)

# Here we set the factor levels and tabulate the data in a single call using pipe operators:
dat.tab03 <- dat.df03 %>%
  mutate(low = factor(low, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(smoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) %>%
  group_by(smoke, low) %>%
  summarise(n = n()) 

# View the data:
dat.tab03

## View the data in conventional 2 by 2 table format:
pivot_wider(dat.tab03, id_cols = c(smoke), 
   names_from = low, values_from = n)

## -----------------------------------------------------------------------------
dat.epi03 <- epi.2by2(dat = dat.tab03, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi03

## -----------------------------------------------------------------------------
dat.df04 <- birthwt; head(dat.df04)

dat.tab04 <- table(dat.df04$smoke, dat.df04$low, dat.df04$race, 
   dnn = c("Smoke", "Low BW", "Race")); dat.tab04

## -----------------------------------------------------------------------------
dat.epi04 <- epi.2by2(dat = dat.tab04, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi04

## -----------------------------------------------------------------------------
dat.df05 <- birthwt; head(dat.df05)

dat.tab05 <- dat.df05 %>%
  mutate(low = factor(low, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(smoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) %>%
  mutate(race = factor(race)) %>%
  group_by(race, smoke, low) %>%
  summarise(n = n()) 
dat.tab05

## View the data in conventional 2 by 2 table format:
pivot_wider(dat.tab05, id_cols = c(race, smoke), 
   names_from = low, values_from = n)

## -----------------------------------------------------------------------------
dat.epi05 <- epi.2by2(dat = dat.tab05, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi05

## -----------------------------------------------------------------------------
library(ggplot2); library(scales)

nstrata <- 1:length(unique(dat.tab05$race))
strata.lab <- paste("Strata ", nstrata, sep = "")
y.at <- c(nstrata, max(nstrata) + 1)
y.lab <- c("M-H", strata.lab)
x.at <- c(0.25,0.5,1,2,4,8,16,32)

or.p <- c(dat.epi05$massoc.detail$OR.mh$est, 
   dat.epi05$massoc.detail$OR.strata.cfield$est)
or.l <- c(dat.epi05$massoc.detail$OR.mh$lower, 
   dat.epi05$massoc.detail$OR.strata.cfield$lower)
or.u <- c(dat.epi05$massoc.detail$OR.mh$upper, 
   dat.epi05$massoc.detail$OR.strata.cfield$upper)
gdat.df05 <- data.frame(y.at, y.lab, or.p, or.l, or.u)

ggplot(data = gdat.df05, aes(x = or.p, y = y.at)) +
   geom_point() + 
   geom_errorbarh(aes(xmax = or.l, xmin = or.u, height = 0.2)) + 
   labs(x = "Odds ratio", y = "Strata") + 
   scale_x_continuous(trans = log2_trans(), breaks = x.at, 
      limits = c(0.25,32)) + 
   scale_y_continuous(breaks = y.at, labels = y.lab) + 
   geom_vline(xintercept = 1, lwd = 1) + 
   coord_fixed(ratio = 0.75 / 1) + 
   theme(axis.title.y = element_text(vjust = 0))

