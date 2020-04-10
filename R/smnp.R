smnp <- function(formula, data = parent.frame(), choiceX = NULL,
                cXnames = NULL, base = NULL, latent = FALSE,
                invcdf = FALSE, n.draws = 5000, p.var = "Inf", 
                p.df = n.dim+1, p.scale = 1, coef.start = 0, b0 = NULL, 
                cov.start = 1, burnin = 0, thin = 0, corrected = 1, verbose = FALSE) {   
  call <- match.call()
  mf <- match.call(expand = FALSE)
  mf$choiceX <- mf$cXnames <- mf$n.draws <- mf$latent <-
    mf$p.var <- mf$p.df <- mf$p.scale <- mf$coef.start <- mf$invcdf <-
      mf$cov.start <- mf$verbose <- mf$burnin <- mf$thin <- mf$corrected <- NULL   
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)


  
  ## obtaining Y
  tmp <- ymatrix.sMNP(mf, extra=TRUE, verbose=verbose)
  Y <- tmp$Y
  lev <- tmp$lev
#  base <- tmp$base
  p <- tmp$p
  n.dim <- p  
  if (p < 3)
    stop("The number of alternatives should be at least 3.")
  if(verbose) 
    cat("The total number of alternatives is ", p, ".\n\n", sep="") 
    
    if(! is.null(b0)){
    	b0 <- sample(1:p, size = 1) - 1
    	}
  
  ### obtaining X
  tmp <- xmatrix.sMNP(formula, data=eval.parent(data),
                     choiceX=call$choiceX, cXnames=cXnames, 
                     n.dim=n.dim, lev=lev, 
                     verbose=verbose, extra=TRUE)
  X <- tmp$X
  coefnames <- tmp$coefnames
  n.cov <- ncol(X) / n.dim
  
  ## listwise deletion for X
  na.ind <- apply(is.na(X), 1, sum)
  if (ncol(Y) == 1)
    na.ind <- na.ind + is.na(Y)
  Y <- Y[na.ind==0,]
  X <- X[na.ind==0,]
  n.obs <- nrow(X)
  
  if (verbose) {
    cat("The dimension of beta is ", n.cov, ".\n\n", sep="")
    cat("The number of observations is ", n.obs, ".\n\n", sep="")
    if (sum(na.ind>0)>0) {
      if (sum(na.ind>0)==1)
        cat("The observation ", (1:length(na.ind))[na.ind>0], " is dropped due to missing values.\n\n", sep="")
      else {
        cat("The following ", sum(na.ind>0), " observations are dropped due to missing values:\n", sep="")
        cat((1:length(na.ind))[na.ind>0], "\n\n")
      }
    }
  } 
  
  ## checking the prior for beta
  p.imp <- FALSE 
  if (p.var == Inf) {
    p.imp <- TRUE
    p.prec <- diag(0, n.cov-tmp$nChoiceSpec/n.dim)
    if (verbose)
      cat("Improper prior will be used for beta.\n\n")
  }
  else if (is.matrix(p.var)) {
    if (ncol(p.var) != n.cov || nrow(p.var) != n.cov)
      stop("The dimension of `p.var' should be ", n.cov, " x ", n.cov, sep="")
    if (sum(sign(eigen(p.var)$values) < 1) > 0)
      stop("`p.var' must be positive definite.")
    p.prec <- solve(p.var)
  }
  else {
    p.var <- diag(p.var, n.cov)
    p.prec <- solve(p.var)
  }
  p.mean <- rep(0, n.cov)

  ## checking prior for Sigma
  p.df <- eval(p.df)
  if (length(p.df) > 1)
    stop("`p.df' must be a positive integer.")
  if (p.df < n.dim)
    stop(paste("`p.df' must be at least ", n.dim, ".", sep=""))
  if (abs(as.integer(p.df) - p.df) > 0)
    stop("`p.df' must be a positive integer.")
  if (!is.matrix(p.scale))  
    p.scale <- diag(1 + 1/(n.dim-1), n.dim-1) - 1/(n.dim - 1)
  if (ncol(p.scale) != n.dim-1 || nrow(p.scale) != n.dim-1)
    stop("`p.scale' must be ", n.dim-1, " x ", n.dim-1, sep="")
  if (sum(sign(eigen(p.scale)$values) < 1) > 0)
    stop("`p.scale' must be positive definite.")
    
    print(p.scale)

  Signames <- NULL
  for(j in 1:n.dim)
    for(k in 1:n.dim)
      if (j<=k)
        Signames <- c(Signames, paste(lev[j+1],
                                      ":", lev[k+1], sep="")) 

  ## checking starting values
  if (length(coef.start) == 1)
    coef.start <- rep(coef.start, n.cov)
  else if (length(coef.start) != n.cov)
    stop(paste("The dimenstion of `coef.start' must be  ",
               n.cov, ".", sep=""))
  if (!is.matrix(cov.start)) {
    cov.start <- diag(1 + 1/(n.dim-1), n.dim-1) - 1/(n.dim - 1)
  }
  else if (ncol(cov.start) != n.dim || nrow(cov.start) != n.dim)
    stop("The dimension of `cov.start' must be ", n.dim, " x ", n.dim, sep="")
  else if (sum(sign(eigen(cov.start)$values) < 1) > 0)
    stop("`cov.start' must be a positive definite matrix.")
  
  ## checking thinnig and burnin intervals
  if (burnin < 0)
    stop("`burnin' should be a non-negative integer.") 
  if (thin < 0)
    stop("`thin' should be a non-negative integer.")
  keep <- thin + 1
  
  ## running the algorithm
  if (latent)
    n.par <- n.cov + n.dim*(n.dim-1)/2 + n.dim*n.obs + 1
  else
    n.par <- n.cov + n.dim*(n.dim-1)/2 + 1
  if(verbose)
    cat("Starting Gibbs sampler corrected (puelz)...\n")
  # recoding NA into -1
  Y[is.na(Y)] <- -1
  
  W <- matrix(rnorm(p * n.obs), nc = p)
  wMaxs <- apply(W, 1, max)
  whichMax <- apply(W, 1, which.max)
#return(list(W=W, wMaxs = wMaxs, whichMax = whichMax, Y=Y))  
  for(i in 1:n.obs){
  	hld <- W[i,Y[i]+1]
  	W[i, Y[i]+1] <- wMaxs[i]
  	W[i, whichMax[i]] <- hld
  }
  
  W <- W - rowMeans(W)
#cat(tmp$nChoiceSpec/n.dim, " nChoiceSpec\n")


  ### mean subtract covariates that vary by choice:
  for(i in 1:length(tmp$nChoiceSpec/n.dim)){
  	nCol <- ncol(X)
  	indices <- (nCol - n.dim + 1):nCol - (i-1) * n.dim
  	tmp2 <- X[,indices]
  	tmp2 <- tmp2 - rowMeans(tmp2)
  	X[,indices] <- tmp2
  	
  }
  
 #print(head(X))
#print(head(W))

#return(list(Y=Y,X=X)) 

  
  # this function is located in the src directory within MNP.c
  param <- .C("cMNPgibbs", as.integer(n.dim), 
              as.integer(tmp$nChoiceSpec/n.dim), as.integer(tmp$nVarByResp), 
              as.integer(n.obs), as.integer(n.draws),
              as.double(p.prec), as.integer(p.df),
              as.double(p.scale), as.double(X), as.double(W), as.integer(b0), 
              as.integer(Y), as.double(coef.start), as.double(cov.start), 
              as.integer(p.imp), as.integer(invcdf),
              as.integer(burnin), as.integer(keep), 
              as.integer(verbose), as.integer(latent), as.integer(corrected),
              pdStore = double(n.par*floor((n.draws-burnin)/keep)),
              PACKAGE="sMNP")$pdStore 
  param <- matrix(param, ncol = n.par,
                  nrow = floor((n.draws-burnin)/keep), byrow=TRUE)
  if (latent) {
    W <- array(as.vector(t(param[,(n.par-n.dim*n.obs+1):n.par])),
               dim = c(n.dim, n.obs, floor((n.draws-burnin)/keep)),
               dimnames = list(lev[-1], rownames(Y), NULL))
    param <- param[,1:(n.par-n.dim*n.obs)]
    }
  else
    W <- NULL
  #colnames(param) <- c(coefnames, Signames)
    
  ##recoding -1 back into NA
  Y[Y==-1] <- NA

  ## returning the object
  res <- list(param = param, x = X, y = Y, w = W, call = call, alt = lev,
              n.alt = p, invcdf = invcdf, 
              p.mean = if(p.imp) NULL else p.mean, p.var = p.var, 
              p.df = p.df, p.scale = p.scale, burnin = burnin, thin = thin) 
  class(res) <- "sMNP"
  return(res)
}
  


