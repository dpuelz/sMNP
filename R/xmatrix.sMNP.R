xmatrix.sMNP <- function(formula, data = parent.frame(), choiceX=NULL,
                        cXnames=NULL, base=NULL, n.dim, lev,
                        verbose=FALSE, extra=FALSE) {
  call <- match.call()
  mf <- match.call(expand = FALSE)
  mf$choiceX <- mf$cXnames <- mf$n.dim <- mf$lev <-
    mf$verbose <- mf$extra <- NULL  
  
  ## get variables
  mf[[1]] <- as.name("model.frame.default")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  X <- model.matrix.default(Terms, mf)
  xvars <- as.character(attr(Terms, "variables"))[-1]
  if ((yvar <- attr(Terms, "response")) > 0)
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  p <- n.dim
  n.obs <- nrow(X)
  n.cov <- ncol(X)
  
  ## expanding X
  allvnames <- Xnew <- NULL
  if (ncol(X) > 0) {
    Xcnames <- colnames(X)
    for (i in 1:n.cov) {
      Xv <- X[, Xcnames[i]]
      Xtmp <- varnames <- NULL
      for (j in 1:n.dim) {
        allvnames <- c(allvnames, paste(Xcnames[i], ":", lev[j], sep=""))
        for (k in 1:n.dim)
          varnames <- c(varnames, paste(Xcnames[i], ":", lev[j], sep=""))
        tmp <- matrix(0, nrow = n.obs, ncol = n.dim)
        tmp[, j] <- Xv
        Xtmp <- cbind(Xtmp, tmp)
      }
      colnames(Xtmp) <- varnames
      Xnew <- cbind(Xnew, Xtmp)
    }
  }
  
  ## checking and adding choice-specific variables
  if (!is.null(choiceX)) {
    cX <- eval(choiceX, data)
    cXn <- unique(names(cX))
    if (sum(is.na(pmatch(cXn, lev))) > 0)
      stop(paste("Error: Invalid input for `choiceX.'\n Some variables do not exist."))
#    if (is.na(pmatch(base, cXn)))
      xbase <- NULL
#    else
#      xbase <- as.matrix(cX[[base]])
    if (length(cXn) < n.dim)
      stop(paste("Error: Invalid input for `choiceX.'\n You must specify the choice-specific varaibles at least for all non-base categories."))
    if (!is.null(xbase) && length(cXn) != p)
      stop(paste("Error: Invalid input for `choiceX.'\n You must specify the choice-specific variables at least for all non-base categories."))
#    if(!is.null(xbase) && verbose)
#      cat("The choice-specific variables of the base category are subtracted from the #corresponding variables of the non-base categories.\n\n")
    for (i in 1:length(cXnames)) 
      for (j in 1:n.dim) {
        if (length(cXnames) != ncol(as.matrix(cX[[lev[j]]])))
            stop(paste("Error: The number of variables in `choiceX' and `cXnames' does not match."))  
        tmp <- matrix(as.matrix(cX[[lev[j]]])[,i], ncol=1)
        colnames(tmp) <- paste(cXnames[i], ":", lev[j], sep="") 
        Xnew <- cbind(Xnew, tmp)
      }
  }
  
  
  if(extra)
    return(list(X=Xnew, coefnames=c(allvnames, cXnames), nChoiceSpec = length(cX), nVarByResp = n.cov))
  else
    return(Xnew)
}
