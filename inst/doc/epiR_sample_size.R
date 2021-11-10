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

## ----message = FALSE----------------------------------------------------------
library(epiR)
epi.sssimpleestb(N = 1E+06, Py = 0.15, epsilon = 0.20, error = "relative", se = 1, sp = 1, nfractional = FALSE, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.sscohortt(irexp1 = 50/1000, irexp0 = 70/1000, FT = 5, n = NA, power = 0.80, r = 1, 
   design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95)$n.total

## ----message = FALSE----------------------------------------------------------
epi.sscc(OR = 2.0, p1 = NA, p0 = 0.30, n = NA, power = 0.80, 
   r = 1, phi.coef = 0, design = 1, sided.test = 2, conf.level = 0.95, 
   method = "unmatched", nfractional = FALSE, fleiss = FALSE)$n.total

## ----message = FALSE----------------------------------------------------------
epi.ssninfb(treat = 0.85, control = 0.65, delta = 0.10, n = NA, 
   r = 1, power = 0.80, nfractional = FALSE, alpha = 0.05)$n.total

## ----message = FALSE----------------------------------------------------------
epi.ssclus1estb(b = 75, Py = 0.46, epsilon = 0.10, error = "relative", rho = 0.20, conf.level = 0.95)$n.psu

## ----message = FALSE----------------------------------------------------------
epi.ssclus1estb(b = c(75,35), Py = 0.46, epsilon = 0.10, error = "relative", rho = 0.20, conf.level = 0.95)$n.psu

## ----message = FALSE----------------------------------------------------------
# From first principles:
n.crude <- epi.sssimpleestb(N = 1E+06, Py = 0.20, epsilon = 0.05 / 0.20, 
   error = "relative", se = 1, sp = 1, nfractional = FALSE, conf.level = 0.95)
n.crude

# A total of 246 subjects need to be enrolled into the study. Calculate the design effect:
rho <- 0.02; b <- 20
D <- rho * (b - 1) + 1; D
# The design effect is 1.38. Our crude sample size estimate needs to be increased by a factor of 1.38.

n.adj <- ceiling(n.crude * D)
n.adj
# After accounting for lack of independence in the data a total of 340 subjects need to be enrolled into the study. How many clusters are required?

ceiling(n.adj / b)

# Do all of the above using epi.ssclus2estb:
epi.ssclus2estb(b = 20, Py = 0.20, epsilon = 0.05 / 0.20, error = "relative", 
   rho = 0.02, nfractional = FALSE, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
n.crude <- epi.sssimpleestb(N = 1E+06, Py = 0.15, epsilon = 0.20,
   error = "relative", se = 1, sp = 1, nfractional = FALSE, conf.level = 0.95)
n.crude

rho <- 0.09; b <- 10
D <- rho * (b - 1) + 1; D

n.adj <- ceiling(n.crude * D)
n.adj

# Similar to the example above, we can do all of these calculations using epi.ssclus2estb:
epi.ssclus2estb(b = 10, Py = 0.15, epsilon = 0.20, error = "relative", 
   rho = 0.09, nfractional = FALSE, conf.level = 0.95)

