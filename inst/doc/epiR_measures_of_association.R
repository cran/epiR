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

## ---- message = FALSE---------------------------------------------------------
library(epiR)

epi.2by2(dat = dat.v01, method = "cross.sectional", conf.level = 0.95, units = 100, 
   interpret = FALSE, outcome = "as.columns")

## ---- message = FALSE---------------------------------------------------------
epi.2by2(dat = dat.v01, method = "cross.sectional", conf.level = 0.95, units = 100, 
   interpret = TRUE, outcome = "as.columns")

## ---- message = FALSE---------------------------------------------------------
library(MASS)

# Load and view the data:
dat.df02 <- birthwt; head(dat.df02)

## ---- message = FALSE---------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02

## ---- message = FALSE---------------------------------------------------------
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW")); dat.tab02
dat.tab02 <- dat.tab02[2:1,2:1]; dat.tab02

## ---- message = FALSE---------------------------------------------------------
# Variables low, smoke and race as factors. Put an 'f' in front of the variable names to remind you that they're factors:

dat.df02$flow <- factor(dat.df02$low, levels = c(1,0))
dat.df02$fsmoke <- factor(dat.df02$smoke, levels = c(1,0))
dat.df02$frace <- factor(dat.df02$race, levels = c(1,2,3))

dat.tab02 <- table(dat.df02$fsmoke, dat.df02$flow, dnn = c("Smoke", "Low BW")); dat.tab02

## ---- message = FALSE---------------------------------------------------------
dat.epi02 <- epi.2by2(dat = dat.tab02, method = "cohort.count", conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi02

## ---- message = FALSE---------------------------------------------------------
names(dat.epi02$massoc.detail)

## ---- message = FALSE---------------------------------------------------------
dat.epi02$massoc.detail$OR.strata.wald
# Wald confidence intervals: 2.02 (95% CI 1.08 to 3.78)

dat.epi02$massoc.detail$OR.strata.score
# Score confidence intervals: 2.02 (95% CI 1.08 to 3.77)


## ---- message = FALSE---------------------------------------------------------
library(tidyverse)

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

## ---- message = FALSE---------------------------------------------------------
dat.epi03 <- epi.2by2(dat = dat.tab03, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi03

## ---- message = FALSE---------------------------------------------------------
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
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
  
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

gdat.df04 <- data.frame(yat = 1:3, ylab = rfactor, ref, or.p, or.l, or.u)
gdat.df04

## ---- message = FALSE---------------------------------------------------------
library(ggplot2); library(scales)

x.at <- c(0.25,0.5,1,2,4,8,16,32)

ggplot(data = gdat.df04, aes(x = or.p, y = yat)) +
  theme_bw() +
  geom_point() + 
  geom_errorbarh(aes(xmax = or.l, xmin = or.u, height = 0.2)) + 
  scale_x_continuous(trans = log2_trans(), breaks = x.at, limits = c(0.25,8), 
   name = "Odds ratio") + 
  scale_y_continuous(breaks = gdat.df04$yat, labels = gdat.df04$ylab, 
   name = "Risk factor") + 
  geom_vline(xintercept = 1, lwd = 1) + 
  annotate("text", x = 0.25, y = gdat.df04$yat, label = gdat.df04$ref, hjust = 0, size = 3) +
  coord_fixed(ratio = 0.75 / 1) + 
  theme(axis.title.y = element_text(vjust = 0))

## ---- message = FALSE---------------------------------------------------------
dat.df05 <- birthwt; head(dat.df05)

dat.df05$flow <- factor(dat.df05$low, levels = c(1,0))

dat.df05$fsmoke <- factor(dat.df05$smoke, levels = c(1,0))
dat.df05$frace <- factor(dat.df05$race, levels = c(1,2,3))

dat.tab05 <- table(dat.df05$fsmoke, dat.df05$flow, dat.df05$frace, 
   dnn = c("Smoke", "Low BW", "Race")); dat.tab05

## ---- message = FALSE---------------------------------------------------------
dat.epi05 <- epi.2by2(dat = dat.tab05, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi05

## ---- message = FALSE---------------------------------------------------------
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

## ---- message = FALSE---------------------------------------------------------
dat.epi06 <- epi.2by2(dat = dat.tab06, method = "cohort.count", 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi06

## ---- message = FALSE---------------------------------------------------------
nstrata <- 1:length(unique(dat.tab06$frace))
strata.lab <- paste("Strata ", nstrata, sep = "")
y.at <- c(nstrata, max(nstrata) + 1)
y.lab <- c("M-H", strata.lab)
x.at <- c(0.25,0.5,1,2,4,8,16,32)

or.p <- c(dat.epi06$massoc.detail$OR.mh$est, 
   dat.epi06$massoc.detail$OR.strata.cfield$est)
or.l <- c(dat.epi06$massoc.detail$OR.mh$lower, 
   dat.epi05$massoc.detail$OR.strata.cfield$lower)
or.u <- c(dat.epi06$massoc.detail$OR.mh$upper, 
   dat.epi05$massoc.detail$OR.strata.cfield$upper)
gdat.df06 <- data.frame(y.at, y.lab, or.p, or.l, or.u)

ggplot(data = gdat.df06, aes(x = or.p, y = y.at)) +
  theme_bw() + 
  geom_point() + 
  geom_errorbarh(aes(xmax = or.l, xmin = or.u, height = 0.2)) + 
  labs(x = "Odds ratio", y = "Strata") + 
  scale_x_continuous(trans = log2_trans(), breaks = x.at, 
   limits = c(0.25,32)) + 
  scale_y_continuous(breaks = y.at, labels = y.lab) + 
  geom_vline(xintercept = 1, lwd = 1) + 
  coord_fixed(ratio = 0.75 / 1) + 
  theme(axis.title.y = element_text(vjust = 0))

