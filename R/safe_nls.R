# Function to perform safe NLS fitting
safe_nls <- function(x, KmGuess=200) {
  # Create a unique starting guess for Vmax
  VmaxGuess <- max(x$v0.mass.norm)
  
  # Safely fit the plots by NLS
  mmFit <- tryCatch(
    nls(formula=mmForm, data=x, start=list(Vmax=VmaxGuess, Km=KmGuess)),
    error=function(err) NA)
  
  mmFit
}