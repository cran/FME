\name{sensFun}
\alias{sensFun}
\alias{summary.sensFun}
\alias{plot.sensFun}
\alias{plot.summary.sensFun}
\alias{pairs.sensFun}

\title{
  Local Sensitivity Analysis
}

\description{
  Given a model consisting of differential equations, estimates the
  local effect of certain parameters on selected sensitivity variables
  by calculating a matrix of so-called sensitivity functions.  In this
  matrix the (i,j)-th element contains

  \deqn{\frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta \Theta _j}
  {\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i}

  and where \eqn{y_i}{y_i} is an output variable (at a certain time instance),
  \eqn{\Theta_j}{par_j} is a parameter, and
  \eqn{\Delta y_i}{varscale_i} is the scaling of variable \eqn{y_i}{y_i},
  \eqn{\Delta \Theta_j}{parscale_j} is the scaling of parameter
  \eqn{\Theta_j}{par_j}.
}
\usage{
sensFun(func, parms, sensvar = NULL, senspar = names(parms),
  varscale = NULL, parscale = NULL, tiny = 1e-8, map = 1, ...)

\method{summary}{sensFun}(object, vars = FALSE, \dots)

\method{pairs}{sensFun}(x, which = NULL, \dots)

\method{plot}{sensFun}(x, which = NULL, legpos="topleft", ask = NULL, \dots)

\method{plot}{summary.sensFun}(x, which = 1:nrow(x), \dots)

}

\arguments{
  \item{func }{an R-function that has as first argument \code{parms} and
    that returns a matrix or data.frame with the values of the output
    variables (columns) at certain output intervals (rows), and
    -- optionally -- a mapping variable (by default the first column).
  }
  \item{parms }{parameters passed to \code{func}; should be either a vector,
    or a list with named elements.
    If \code{NULL}, then the first element of \code{parInput} is taken.
  }
  \item{sensvar }{the output variables for which the sensitivity needs
    to be estimated. Either \code{NULL}, the default, which selects all
    variables, or a vector with variable \code{names} (which should be
    present in the matrix returned by \code{func}), or a vector with
    \code{indices} to variables as present in the output matrix (note
    that the column of this matrix with the mapping variable should not
    be selected).
  }
  \item{senspar }{the parameters whose sensitivity needs to be
    estimated, the default=all parameters. Either a vector with
    parameter \emph{names}, or a vector with \emph{indices} to positions
    of parameters in \code{parms}.
  }
  \item{varscale }{the scaling (weighing) factor for sensitivity
    variables, \code{NA} indicates that the variable value is used.
  }
  \item{parscale }{the scaling (weighing) factor for sensitivity
    parameters, \code{NA} indicates that the parameter value is used.
  }
  \item{tiny }{the perturbation, or numerical difference, factor, see
    details.
  }
  \item{map }{the column number with the (independent) mapping variable
    in the output matrix returned by \code{func}. For dynamic models
    solved by integration, this will be the (first) column with
    \code{time}. For 1-D spatial output, this column will be some
    distance variable.  Set to NULL if there is no mapping
    variable. Mapping variables should not be selected for estimating
    sensitivity functions; they are used for plotting.
  }
  \item{... }{additional arguments passed to \code{func} or to the
    methods.
  }
  \item{object }{an object of class \code{sensFun}.
  }
  \item{x }{an object of class \code{sensFun}.
  }
  \item{vars }{if FALSE: summaries per parameter are returned; if
    \code{TRUE}, summaries per parameter and per variable are returned.
  }
  \item{which }{the name or the index to the variables that should be
    plotted.  Default = all variables.
  }
  \item{legpos }{position of the legend; set to \code{NULL} to avoid
    plotting a legend.
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask = ...)} and
    \code{\link{dev.interactive}}.
  }
}
\value{
  
  a data.frame of class \code{sensFun} containing the sensitivity
  functions this is one row for each sensitivity variable at each
  independent (time or position) value and the following columns:
  
  \code{x}, the value of the independent (mapping) variable, usually
  time (solver= "ode.."), or distance (solver= "steady.1D")

  \code{var}, the name of the observed variable,

  \code{...}, a number of columns, one for each sensitivity parameter

  The data.frame returned by \code{sensFun} has methods for the generic
  functions  \code{\link{summary}}, \code{\link{plot}},
  \code{\link{pairs}} -- see note.
}

\details{

  There are essentially two ways in which to use function \code{sensFun}.

  \itemize{
    \item When \code{func} returns a matrix or data frame with output
      values, \code{sensFun} can be used for sensitivity analysis,
      estimating the \emph{impact} of parameters on output variables.
       
    \item When \code{func} returns an instance of class \code{modCost}
      (as returned by a call to function \code{\link{modCost}}), then
      \code{sensFun} can be used for \emph{parameter
      identifiability}. In this case the results from \code{sensFun} are
      used as input to function \link{collin}. See the help file for
      \code{collin}.
  }

  For each sensitivity parameter, the number of sensitivity functions
  estimated is: length(sensvar) * length(mapping variable), i.e. one for
  each element returned by \code{func} (except the mapping variable).

  The sensitivity functions are estimated numerically. This means that
  each parameter value \eqn{\Theta_j}{par_j} is perturbed as
  \eqn{\max{(tiny,\Theta_j \cdot (1+tiny))}}{max(tiny,par_j)*(1+tiny)}
}

\note{
  The following methods are provided:
  
  \itemize{
    \item \emph{summary}. Produces summary statistics of the sensitivity
      functions, a \code{data.frame} with: one row for each parameter
      and the following columns:
  
      \itemize{
        \item L1: the L1-norm \eqn{\frac{1}{n} \cdot \sum{|S_{ij}|}}{sum(abs(Sij))/n},

        \item L2: the L2-norm \eqn{\cdot \sqrt{\frac{1}{n}\sum{S_{ij}
          \cdot S_{ij}}}}{sqrt(sum(Sij^2)/n)},

        \item Mean: the mean of the sensitivity functions,

        \item Min: the minimal value of the sensitivity functions,

        \item Max: the maximal value of the sensitivity functions.
      }
    \item \emph{var }the summary of the variables sensitivity functions, a
      data.frame with the same columns as \code{model}
      and one row for each parameter + variable combination.
      This is only outputted if the variable names are effectively known

    \item \emph{plot } plots the sensitivity functions for each parameter;
      each parameter has its own color.
    
      By default, the sensitivity functions for all variables are plotted
      in one figure, unless \code{which} gives a selection of variables;
      in that case, each variable will be plotted in a separate figure, and
      the figures aligned in a rectangular grid, unless par \code{mfrow} is
      passed as an argument.

    \item \emph{pairs } produces a pairs plot of the sensitivity results;
      per parameter.

      By default, the sensitivity functions for all variables are plotted
      in one figure, unless \code{which} gives a selection of variables.
      
      Overrides the default \code{gap = 0}, \code{upper.panel = NA}, and
      \code{diag.panel}.
  }
}

\examples{
## =======================================================================
## Bacterial growth model as in Soetaert and Herman, 2009
## =======================================================================
pars <- list(gmax = 0.5, eff = 0.5,
              ks = 0.5, rB = 0.01, dB = 0.01)

solveBact <- function(pars) {
  derivs <- function(t, state, pars) { # returns rate of change
    with (as.list(c(state, pars)), {
      dBact <-  gmax * eff * Sub/(Sub + ks) * Bact - dB * Bact - rB * Bact
      dSub  <- -gmax       * Sub/(Sub + ks) * Bact + dB * Bact
      return(list(c(dBact, dSub)))
    })
  }
  state   <- c(Bact = 0.1, Sub = 100)
  tout    <- seq(0, 50, by = 0.5)
  ## ode solves the model by integration ...
  return(as.data.frame(ode(y = state, times = tout, func = derivs,
    parms = pars)))
}

out <- solveBact(pars)

plot(out$time, out$Bact, ylim = range(c(out$Bact, out$Sub)),
     xlab = "time, hour", ylab = "molC/m3", type = "l", lwd = 2)
lines(out$time, out$Sub,lty = 2,lwd = 2)
lines(out$time, out$Sub + out$Bact)

legend("topright", c("Bacteria", "Glucose", "TOC"),
       lty = c(1, 2, 1), lwd = c(2, 2, 1))

## sensitivity functions
SnsBact <- sensFun(func = solveBact, parms = pars,
                   sensvar = "Bact", varscale = 1)
head(SnsBact)
plot(SnsBact)
summary(SnsBact)
plot(summary(SnsBact))

SF <- sensFun(func = solveBact, parms = pars,
             sensvar = c("Bact", "Sub"), varscale = 1)
head(SF)
tail(SF)

summary(SF, var = TRUE)

plot(SF)
plot(SF, which = c("Sub","Bact"))
pm <- par(mfrow = c(1,3))
plot(SF, which = c("Sub", "Bact"), mfrow = NULL)
plot(SF, mfrow = NULL)
par(mfrow = pm)

## Bivariate sensitivity
pairs(SF)  # same color
pairs(SF, which = "Bact", col = "green", pch = 15)
pairs(SF, which = c("Bact", "Sub"), col = c("green", "blue"))
mtext(outer = TRUE, side = 3, line = -2,
      "Sensitivity functions", cex = 1.5)

## pairwise correlation
cor(SnsBact[,-(1:2)])

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\references{
  Soetaert, K. and Herman, P. M. J., 2009. A Practical Guide to Ecological
  Modelling -- Using R as a Simulation Platform. Springer, 390 pp.

  Brun, R., Reichert, P. and Kunsch, H.R., 2001.  Practical
  Identificability Analysis of Large Environmental Simulation Models.
  Water Resour. Res. 37(4): 1015--1030.

  Soetaert, K. and Petzoldt, T., 2010.  Inverse Modelling, Sensitivity
  and Monte Carlo Analysis in R Using Package FME.  Journal of
  Statistical Software 33(3) 1--28.
  \url{http://www.jstatsoft.org/v33/i03}
}

\keyword{utilities}

