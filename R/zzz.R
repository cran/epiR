.First.lib <- function(libname, pkgname)
{
   ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version")
   ver <- as.character(ver)	
   message("Package epiR ", ver, " is loaded")
   message("Type help(epi.about) for summary information")
   message("\n")
}