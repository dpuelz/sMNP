ymatrix.sMNP <- function(data, base=NULL, extra=FALSE, verbose=verbose) { 
  ## checking and formatting Y
  Y <- model.response(data)
 # standard Multinomial Probit model        
    Y <- as.factor(Y)
    lev <- levels(Y)
    counts <- table(Y)
    if (any(counts == 0)) {
      warning(paste("group(s)", paste(lev[counts == 0], collapse = " "), "are empty"))
      Y <- factor(Y, levels  = lev[counts > 0])
      lev <- lev[counts > 0]
    }
    p <- length(lev)
    Y <- as.matrix(unclass(Y)) - 1
  
  if(extra)
    return(list(Y=Y, lev=lev, p=p))
  else
    return(Y)
}
