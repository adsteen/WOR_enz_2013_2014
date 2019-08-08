predict_mm <- function(nls_obj, s=0:400) {
  
  if(!(class(nls_obj) == "nls" | class(nls_obj) == "nls2")) {
    if(is.na(nls_obj)) {
      preds <- rep(NA, length(s))
    } else {
      stop("The argument is neither of class nls or nls2, nor is it NA. I don't know what to do.")
    }
  } else{
    
    newdata <- data.frame(conc=s)
    
    preds <- predict(nls_obj, newdata)
    
    # # S is vector, Vmax & Km are constants
    # if(length(Vmax)!=1 | !is.numeric(Vmax)) {
    #   stop("Vmax should be a numeric vector of length 1")
    # }
    # if(length(Km)!=1 | !is.numeric(Km)) {
    #   stop("Vmax should be a numeric vector of length 1")
    # }
  }

  # Calculate v0 based on Michaelis-Menten curve
  data.frame(s=s, preds)
}

