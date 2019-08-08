## Pulls out coefficients of MM mod

safe_coef <- function(nls_obj, param, est.or.err) {
  
  # Check reasonable object. One day I'll rewrite this with nls2
  if(!(class(nls_obj) == "nls" | class(nls_obj) == "nls2")) {
    if(is.na(nls_obj)) {
      result <- NA
      return(result)
    } else {
      stop("The argument is neither of class nls or nls2, nor is it NA. I don't know what to do.")
    }
  }
  
  sum <- summary(nls_obj) 
  #browser()
  
  result <- sum$parameters[param, est.or.err]
  
  result
}
