\name{sensRange}
\alias{sensRange}
\alias{summary.sensRange}
\alias{plot.sensRange}
\alias{plot.summary.sensRange}
\title{
  Sensitivity Ranges of a Timeseries or 1-D Variables
}

\description{
  Given a model consisting of differential equations, estimates
  the global effect of certain (sensitivity) parameters on a time series or on
  1-D spatial series of selected sensitivity variables.

  This is done by drawing parameter values according to some predefined
  distribution, running the model with each of these parameter combinations,
  and calculating the values of the selected output variables at each output
  interval.

  This function thus produces 'envelopes' around the sensitivity variables.
}
\usage{
sensRange(func,parms=NULL,sensvar=NULL,dist="unif",
  parInput=NULL, parRange=NULL, parMean=NULL, parCovar=NULL,
  map=1, num=100,  ...)
  
\method{summary}{sensRange}(object,  \dots)

\method{plot}{summary.sensRange}(x, xyswap =FALSE,
              which=NULL, legpos="topleft",
              col=c(grey(0.8), grey(0.7)),
              quant =FALSE, ask = NULL, \dots)

\method{plot}{sensRange}(x, xyswap =FALSE,
              which=NULL, ask=NULL, \dots)
}

\arguments{
  \item{func }{an R-function that has as first argument \code{parms} and that
    returns a matrix or data.frame with the values of the output variables
    (columns) at certain output intervals (rows), and -optionally- a mapping
    variable (by default the first column).
  }
  \item{parms }{parameters passed to \code{func}; should be either a vector,
    or a list with named elements.
    If \code{NULL}, then the first element of \code{parInput} is taken.
  }
  \item{sensvar }{the output variables for which the sensitivity needs to be
    estimated. Either \code{NULL}, the default, which selects all variables,
    or a vector with variable \code{names} (which should be present
    in the matrix returned by \code{func}), or a vector with \code{indices}
    to variables as present in the output matrix (note that the column of this
    matrix with the mapping variable should not be selected).
  }
  \item{dist }{the distribution according to which the parameters should be
    generated, one of \code{"unif"} (uniformly random samples),
    \code{"norm"}, (normally distributed random samples),
    \code{"latin"} (latin hypercube distribution),
    \code{"grid"} (parameters arranged on a grid).
    The input parameters for the distribution are specified by \code{parRange}
    (min,max), except for the normally distributed parameters, in which case
    the distribution is specified by the parameter means \code{parMean} and
    the variance-covariance matrix, \code{parCovar}. Note that, if the
    distribution is \code{"norm"} and \code{parRange} is given, then a
    truncated distribution will be generated.
    (This is useful to prevent for instance that certain parameters become
    negative). Ignored if \code{parInput} is specified.
  }
  \item{parRange }{the range (min, max) of the sensitivity parameters, a matrix
    or (preferred) a data.frame with one row for each parameter, and two
    columns with the minimum (1st) and maximum (2nd) value.
    The rownames of \code{parRange} should be parameter names that are known
    in argument \code{parms}. Ignored if \code{parInput} is specified.
  }
  \item{parInput }{a matrix with dimension (*, npar) with the values of the
    sensitivity parameters.
  }
  \item{parMean }{only when \code{dist} is \code{"norm"}: the mean value of
    each parameter. Ignored if \code{parInput} is specified.
  }
  \item{parCovar }{only when \code{dist} is \code{"norm"}: the parameter's
    variance-covariance matrix.
  }
  \item{num }{the number of times the model has to be run. Set large enough.
    If \code{parInput} is specified, then \code{num} parameters are selected
    randomly (from the rows of \code{parInput}.
  }
  \item{map }{the column number with the (independent) mapping variable in the
    output matrix returned by \code{func}. For dynamic models solved by
    integration, this will be the (first) column with \code{time}. For 1-D
    spatial output, this column will be some distance variable.
    Set to NULL if there is no mapping variable. Mapping variables should not
    be selected for estimating sensitivity ranges; they are used for plotting.
  }
  \item{object }{an object of class \code{sensRange}.
  }
  \item{x }{an object of class \code{sensRange}.
  }
  \item{legpos }{position of the legend; set to \code{NULL} to avoid plotting
  a legend.
  }
  \item{xyswap }{if \code{TRUE}, then x-and y-values are swapped and the y-axis
    is from top to bottom. Useful for drawing vertical profiles.
  }
  \item{which }{the name or the index to the variables that should be plotted.
    Default = all variables.
  }
  \item{col }{the two colors of the polygons that should be plotted.
  }
  \item{quant }{if \code{TRUE}, then the median surrounded by the quantiles
    q25-q75 and q95-q95 are plotted, else the min-max and mean+- sd are plotted.
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.
  }
  \item{... }{additional arguments passed to \code{func} or to the methods.
  }
}
\value{
  a data.frame of type \code{sensRange} containing the parameter set and
  the corresponding values of the sensitivity output variables.
  
  The list returned by \code{sensRange} has a method for the generic functions
   \code{\link{summary}},\code{\link{plot}} and \code{plot.summary}
  - see note.
  
}

\details{
  Models solved by integration (i.e. by using one of \code{'ode', 'ode.1D',
  'ode.band', 'ode.2D'}), have the output already in a form usable by
  \code{sensRange}.

  The output for models solved by a steady-state solver (i.e. one of
  \code{'steady', 'steady.1D','steady.band', 'steady.2D'}, needs to be
  rearranged - see examples.
}
\note{
  The following \emph{methods} are included:
  \itemize{
    \item \link{summary}, estimates summary statistics for the sensitivity
      variables, a data.frame with as many rows as there are mapping variables
      (or rows in the matrix returned by \code{func})
      and the following columns: \code{x}, the mapping value, \code{Mean},
      the mean, \code{sd}, the standard deviation, \code{Min}, the minimal
      value, \code{Max}, the maximal value, \code{q25},  \code{q50}, \code{q75},
      the 25th, 50 and 75\% quantile
  
    \item \link{plot}, produces a "matplot" of the \code{sensRange} output,
      one plot for each sensitivity variable and with the mapping variable on
      the x-axis.

      Each variable will be plotted in a separate figure, and
      the figures aligned in a rectangular grid, unless par \code{mfrow} is
      passed as an argument.

    \item summary.plot, produces a plot of the summary of the \code{sensRange}
      output, one plot for each sensitivity variable and with the ranges and
      mean +- standard deviation or the quantiles as coloured polygons.
      
      Each variable will be plotted in a separate figure, and
      the figures aligned in a rectangular grid, unless par \code{mfrow} is
      passed as an argument.

  }

  The data.frame of type \code{sensRange} has several attributes, which
  remain hidden, and which are generally not of practical use
  (they are needed for the S3 methods).

  There is one exception, i.e. if parameter values are imposed via
  argument \code{parInput}, and these parameters are generated by a Markov chain
  (\code{\link{modMCMC}}). If the number of draws, \code{num}, is less
  than the number of rows in \code{parInput}, then \code{num} random draws
  will be taken. Attribute, "pset" then contains the index to the parameters that
  have been selected.
  The \code{sensRange} method only represents the distribution of the model
  response variables as a function of the parameter values. But an additional
  source of noise is due to the \emph{model error}, as represented by the
  sampled values of sigma in the Markov chain.
  In order to represent also this source of error, gaussian noise should be
  added to each sensitivity output variables, with a standard deviation that
  corresponds to the original parameter draw - see vignette "FMEother".

}
\examples{

## =======================================================================
## Bacterial growth model from Soetaert and Herman, 2009
## =======================================================================

pars <- list(gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

solveBact <- function(pars)
{
 derivs <- function(t,state,pars)     # returns rate of change
  {
  with (as.list(c(state,pars)), {

  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact

  return(list(c(dBact,dSub)))
                              })
  }

 state   <- c(Bact=0.1,Sub = 100)
 tout    <- seq(0,50,by=0.5)
 # ode solves the model by integration...
 return(as.data.frame(ode(y=state,times=tout,func=derivs,parms=pars)))
}

out <- solveBact(pars)

mf  <-par(mfrow=c(2,2))

plot(out$time,out$Bact,main="Bacteria",,
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)

# the sensitivity parameters
parRanges <- data.frame(min=c(0.4,0.4,0.),max=c(0.6,0.6,0.02))
rownames(parRanges)<- c("gmax","eff","rB")
parRanges

tout    <- 0:50
# sensitivity to rB; equally-spaced parameters ("grid")
SensR <- sensRange(func=solveBact,parms=pars,dist="grid",
                   sensvar="Bact",parRange=parRanges[3,],num=50)

Sens  <-summary(SensR)
plot(Sens,legpos="topleft",xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to rB",mfrow=NULL)

# sensitivity to all; latin hypercube
Sens2 <- summary(sensRange(func=solveBact,parms=pars,dist="latin",
                   sensvar=c("Bact","Sub"),parRange=parRanges,num=50))

# Plot all variables; plot mean +- sd, min max
plot(Sens2,xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to gmax,eff,rB",mfrow=NULL)

par(mfrow=mf )

# Select one variable for plotting; plot the quantiles
plot(Sens2,xlab="time, hour",ylab="molC/m3", which="Bact", quant = TRUE)

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

