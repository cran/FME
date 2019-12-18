\name{FME-package}
\alias{FME-package}
\alias{FME}
\docType{package}
\title{
  A Flexible Modelling Environment for Inverse Modelling, Sensitivity,
  Identifiability, Monte Carlo Analysis.
}
\description{
  R-package FME contains functions to run complex applications of models
  that produce output as a function of input parameters.

  Although it was created to be used with models consisting of ordinary
  differential equations (ODE), partial differential equations (PDE) or
  differential algebraic equations (DAE), it can work with other models.

  It contains:
  \itemize{
    \item Functions to allow fitting of the model to data.

      Function \code{modCost} estimates the (weighted) residuals between
      model output and data, variable and model costs.

      Function \code{modFit} uses the output of \code{modCost} to find
      the best-fit parameters. It provides a wrapper around \R's
      built-in minimisation routines (\code{optim}, \code{nlm},
      \code{nlminb}) and \code{nls.lm} from package \code{minpack.lm}.

      Package FME also includes an implementation of the pseudo-random
      search algorithm (function \code{\link{pseudoOptim}}).

    \item Function \code{sensFun} estimates the sensitivity functions of
      selected output variables as a function of model parameters. This
      is the basis of uni-variate, bi-variate and multi-variate
      sensitivity analysis.

    \item Function \code{collin} uses as input the sensitivity functions
      and estimates the "collinearity" index for all possible parameter
      sets.  This multivariate sensitivity estimate measures approximate
      linear dependence and is useful to derive which parameter sets are
      identifiable given the data set.

    \item Function \code{sensRange} produces 'envelopes' around the
      sensitivity variables, consisting of a time series or a
      1-dimensional set, as a function of the sensitivity parameters.
      It produces "envelopes" around the variables.

    \item Function \code{modCRL} calculates the values of single
      variables as a function of the sensitivity parameters. This
      function can be used to run simple "what-if" scenarios

    \item Function \code{modMCMC} runs a Markov chain Monte Carlo
      (Bayesian analysis). It implements the delayed rejection -
      adaptive Metropolis (DRAM) algorithm.

    \item FME also contains functions to generate multiple parameter
      values arranged according to a grid (\code{Grid}) multinormal
      (\code{Norm}) or uniform (\code{Unif}) design, and a latin
      hypercube sampling (\code{Latinhyper}) function }
}
\details{

  bug corrections:
  
  \itemize{

  \item version 1.3.6, sensFun:  corrected calculation of L2 norm 
    (now consistent with help page),
  
  \item version 1.3,  modCost: minlogp was not correctly estimated 
    if more than one observed variable (used the wrong sd).
  }
}
\author{
  Karline Soetaert

  Thomas Petzoldt
}
\references{
  Soetaert, K. and Petzoldt, T. 2010.  Inverse Modelling, Sensitivity and
  Monte Carlo Analysis in R Using Package FME.  Journal of Statistical
  Software 33(3) 1--28. \url{https://doi.org/10.18637/jss.v033.i03}
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(modCost)
example(sensFun)
example(modMCMC)
example(modCRL)


## open the directory with documents
browseURL(paste(system.file(package = "FME"), "/doc", sep = ""))

## open the directory with examples
browseURL(paste(system.file(package = "FME"), "/doc/examples", sep = ""))

## the vignettes
vignette("FME")
vignette("FMEdyna")
vignette("FMEsteady")
vignette("FMEother")
vignette("FMEmcmc")

edit(vignette("FME"))
edit(vignette("FMEdyna"))
edit(vignette("FMEsteady"))
edit(vignette("FMEother"))
edit(vignette("FMEmcmc"))

}
}

\keyword{package}