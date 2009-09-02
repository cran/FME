## -----------------------------------------------------------------------------
## Collinearity indices
## -----------------------------------------------------------------------------

collin <- function(sensfun, parset = NULL, N = NULL, which = NULL) {

  if (is.null(colnames(sensfun)))colnames(sensfun) <- 1:ncol(sensfun)

  if (!is.null(which)) {
    nx  <- attr(sensfun, "nx")
    var <- attr(sensfun, "var")
    TYP <- attr(sensfun, "Type")

    if (! is.numeric(which)) {
      ln <- length(which)
      Select <- which (var %in% which)
      if(length(Select) != ln)
        stop("not all variables in 'which' are in 'sensfun'")
    } else {
       Select <- which
       if (max(Select) > nx)
         stop("index in 'which' too large")
    }
    ii <- NULL

    if (TYP == 1)
     for (i in Select)   ii <- c(ii, ((i-1)*nx):(i*nx))
    else
      for (i in Select)  ii <- c(ii, (nx[i]+1):nx[i+1])
    sensfun <- sensfun[ii,]
  }

  if (colnames(sensfun)[1]=="x" && colnames(sensfun)[2] == "var")
    Sens <- sensfun[,-(1:2)]
  else Sens <- sensfun

  npar <- ncol(Sens)
  L2   <- sqrt(colSums(Sens*Sens))

  iNa <- 0
  ## Check for non-identifiable parameters
  if (any(L2 == 0 ) ) {
    iNa <- which(L2 == 0)
    warning (paste("Sensitivity of parameter", colnames(Sens)[iNa], "is 0! "))
  }
  if (npar > 14 & is.null(parset) & is.null(N))
    warning ("will reduce collinearity estimates: too many combinations")

  normSens <- t(t(Sens) / L2)

  Collin <- NULL
  ## internal function to generate collinearities
  ## for a given set of parameter combinations

  collFun <- function(cc) {
    Collset <- NULL
    for (i in 1:nrow(cc)) {
      ii    <- cc[i,]
      S     <- normSens[,ii]
      Nident <- (iNa != 0 & iNa %in% ii)
      if (Nident) {
        id <- Inf
      } else {
        id  <- 1/sqrt(min(eigen(t(S) %*% S)$value))
      }
      psub     <- rep(0, npar)
      psub[ii] <- 1
      n        <- ncol(cc)
      Collset <- rbind(Collset, c(psub, n, id))
    }
    return(rbind(Collin, Collset))
  }
  if (is.null(parset)) {

    combin <- function(n, v) { # combinations of n elements from a vector p (length (p) > = n)
      if (n == 1)
        matrix(data = v, ncol = 1)
      else if (n >= length(v))
        matrix(data = v, nrow = 1)
      else
        rbind(cbind(v[1], combin(n-1, v[-1])), combin(n, v[-1]))
    }

    pset   <- 1:npar
    if (is.null(N)) nset <- 2:npar else nset <- N
    for (n in nset) {
      numcomb <- choose(npar, n)
      if (numcomb < 5000) {
        cc  <- combin(n, pset)
        Collin <- collFun(cc)
      }
    }
  } else {
    if (! is.vector(parset))
      stop("'parset' should be a vector")
    if (is.character(parset)) {
      ln <- length(parset)
      pnames<-colnames(Sens)
      parset <- which(pnames %in% parset)
      if (length(parset) != ln)
        stop ("Not all parameters in parset known")
    }

    parset <- matrix(data = parset, nrow = 1)
    Collin <- collFun(parset)
  }

  Collin <- as.data.frame(Collin)

  class(Collin) <- c("collin", "data.frame")
  names(Collin) <- c(colnames(Sens), "N", "collinearity")

  return(Collin)
}

## -----------------------------------------------------------------------------
## S3 methods of collin
## -----------------------------------------------------------------------------

plot.collin <- function(x, ...) {

  dots <- list(...)
  dots$ylab <- if(is.null(dots$ylab)) "Collinearity index" else dots$ylab
  dots$xlab <- if(is.null(dots$xlab)) "Number of parameters" else dots$xlab
  dots$main <- if(is.null(dots$main)) "Collinearity" else dots$main

  nc <- ncol(x)
  do.call("stripchart", c(alist(x[,nc] ~ x[,nc-1], method="stack",
    vertical = TRUE), dots))
}


## -----------------------------------------------------------------------------

print.collin <- function (x, ...)
  print(format(as.data.frame(unclass(x)), digits = 2, scientific = FALSE, ...))
