"epi.about" <- function()
{
  cat("\n")
  cat("-----------------------------------------------------------\n")
  ver <- package.description("epi", lib = NULL, field="Version")
  cat(paste("epi version", ver,  "is now loaded"))
  cat("\n")
  cat("-----------------------------------------------------------\n")
  invisible()
}
