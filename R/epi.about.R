"epi.about" <- function()
{
  cat("\n")
  cat("-----------------------------------------------------------\n")
  ver <- packageDescription("epiR", lib.loc = NULL, fields = "Version")
  cat(paste("epiR version", ver))
  cat("\n")
  cat("An R package for the analysis of epidemiological data.")
  cat("\n")
  cat("See http://epicentre.massey.ac.nz/ for details.")
  cat("\n")
  cat("-----------------------------------------------------------\n")
  invisible()
}
