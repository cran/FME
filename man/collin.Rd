\name{collin}
\alias{collin}
\alias{plot.collin}
\alias{print.collin}

\title{
  Estimates the Collinearity of Parameter Sets
}
\description{
  Based on the sensitivity functions of model variables to a selection
  of parameters, calculates the "identifiability" of sets of parameter.

  The sensitivity functions are a matrix whose (i,j)-th element contains
  \deqn{\frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta
  \Theta _j} {\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i} and where
  \eqn{y_i}{y_i} is an output variable, at a certain (time) instance, i,
  \eqn{\Delta y_i}{varscale_i} is the scaling of variable
  \eqn{y_i}{y_i}, \eqn{\Delta \Theta_j}{parscale_j} is the scaling of
  parameter \eqn{\Theta_j}{par_j}.

  Function \code{collin} estimates the collinearity, or identifiability of all
  parameter sets or of one parameter set.

  As a rule of thumb, a collinearity value less than about 20 is
  "identifiable".
}

\usage{
collin(sensfun, parset = NULL, N = NULL, which = NULL, maxcomb = 5000)

\method{print}{collin}(x, ...)

\method{plot}{collin}(x, ...)
}
\arguments{
  \item{sensfun }{model sensitivity functions as estimated by \code{SensFun}.
  }
  \item{parset }{one selected parameter combination, a vector with their names
    or with the indices to the parameters.
  }
  \item{N }{the number of parameters in the set; if \code{NULL} then all
    combinations will be tried. Ignored if \code{parset} is not \code{NULL}.
  }
  \item{which }{the name or the index to the observed variables that should be
    used. Default = all observed variables.
  }
  \item{maxcomb }{the maximal number of combinations that can be tested.
    If too large, this may produce a huge output. The number of combinations of 
    n parameters out of a total of p parameters is \code{choose(p, n)}.
  }
  \item{x }{an object of class \code{collin}.
  }
  \item{... }{additional arguments passed to the methods.
  }
}
\value{
  a data.frame of class \code{collin} with one row for each parameter
  combination (parameters as in \code{sensfun}).

  Each row contains:
  
  \item{...}{for each parameter whether it is present (1) or absent (0)
    in the set,}
  \item{N}{the number of parameters in the set,}
  \item{collinearity}{the collinearity value.}

  The data.frame returned by \code{collin} has methods for the generic
  functions \code{\link{print}} and \code{\link{plot}}.
}

\details{
  The collinearity is a measure of approximate linear dependence between
  sets of parameters. The higher its value, the more the parameters are
  related.  With "related" is meant that several paraemter combinations
  may produce similar values of the output variables.
}

\note{
  It is possible to use \code{collin} for selecting parameter sets that
  can be fine-tuned based on a data set.  Thus it is a powerful
  technique to make model calibration routines more robust, because
  calibration routines often fail when parameters are strongly related.

  In general, when the collinearity index exceeds 20, the linear
  dependence is assumed to be critical (i.e. it will not be possible or
  easy to estimate all the parameters in the combination together).

  The procedure is explained in Omlin et al. (2001).

  1. First the function \code{collin} is used to test how far a dataset
    can be used for estimating certain (combinations of) parameters.
    After selection of an 'identifiable parameter set' (which has a low
    "collinearity") they are fine-tuned by calibration.
 
  2. As the sensitivity analysis is a \emph{local} analysis (i.e. its
    outcome depends on the current values of the model parameters) and
    the fitting routine is used to estimate the best values of the
    parameters, this is an iterative procedure.  This means that
    identifiable parameters are determined, fitted to the data, then a
    newly identifiable parameter set is determined, fitted, etcetera
    until convergenc is reached.

  See the paper by Omlin et al. (2001) for more information.
}
\examples{
## =======================================================================
## Test collinearity values
## =======================================================================

## linearly related set...  => Infinity
collin(cbind(1:5, 2*(1:5)))

## unrelated set            => 1
MM <- matrix(nr = 4, nc = 2, byrow = TRUE,
  data = c(-0.400, -0.374, 0.255, 0.797, 0.690, -0.472, -0.546,  0.049))

collin(MM)

## =======================================================================
## Bacterial model as in Soetaert and Herman, 2009
## =======================================================================

pars <- list(gmax = 0.5,eff = 0.5,
             ks = 0.5, rB = 0.01, dB = 0.01)

solveBact <- function(pars) {
  derivs <- function(t, state, pars) {   # returns rate of change
    with (as.list(c(state, pars)), {
      dBact <-  gmax*eff*Sub/(Sub + ks)*Bact - dB*Bact - rB*Bact
      dSub  <- -gmax    *Sub/(Sub + ks)*Bact + dB*Bact
      return(list(c(dBact, dSub)))
    })
  }
  state   <- c(Bact = 0.1, Sub = 100)
  tout    <- seq(0, 50, by = 0.5)
  ## ode solves the model by integration...
  return(as.data.frame(ode(y = state, times = tout, func = derivs,
    parms = pars)))
}

out <- solveBact(pars)

## We wish to estimate parameters gmax and eff by fitting the model to
## these data:
Data <- matrix(nc = 2, byrow = TRUE, data =
  c(  2,  0.14,  4,  0.2,    6,  0.38,  8,  0.42,
     10,  0.6,  12,  0.107, 14,  1.3,  16,  2.0,
     18,  3.0,  20,  4.5,   22,  6.15, 24,  11,
     26, 13.8,  28, 20.0,   30,  31 ,  35, 65, 40, 61)
)
colnames(Data) <- c("time","Bact")
head(Data)

Data2 <- matrix(c(2, 100, 20, 93, 30, 55, 50, 0), ncol = 2, byrow = TRUE)
colnames(Data2) <- c("time", "Sub")


## Objective function to minimise
Objective <- function (x) {                  # Model cost
 pars[] <- x
 out   <- solveBact(x)
 Cost  <- modCost(obs = Data2, model = out)  # observed data in 2 data.frames
 return(modCost(obs = Data, model = out, cost = Cost))
}

## 1. Estimate sensitivity functions - all parameters
sF <- sensFun(func = Objective, parms = pars, varscale = 1)

## 2. Estimate the collinearity
Coll <- collin(sF)

## The larger the collinearity, the less identifiable the data set
Coll

plot(Coll, log = "y")

## 20 = magical number above which there are identifiability problems
abline(h = 20, col = "red")

## select "identifiable" sets with 4 parameters
Coll [Coll[,"collinearity"] < 20 & Coll[,"N"]==4,]

## collinearity of one selected parameter set
collin(sF, c(1, 3, 5))
collin(sF, 1:5)

collin(sF, c("gmax", "eff"))
## collinearity of all combinations of 3 parameters
collin(sF, N = 3)

## The collinearity depends on the value of the parameters:
P      <- pars
P[1:2] <- 1  # was: 0.5
collin(sensFun(Objective, P, varscale = 1))

}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\references{
  Brun, R., Reichert, P. and Kunsch, H. R., 2001.
  Practical Identifiability Analysis of Large Environmental Simulation Models.
  Water Resour. Res. 37(4): 1015--1030.

  Omlin, M., Brun, R. and Reichert, P., 2001.
  Biogeochemical Model of Lake Zurich: Sensitivity, Identifiability and
  Uncertainty Analysis. Ecol. Modell. 141: 105--123.

  Soetaert, K. and Petzoldt, T. 2010.  Inverse Modelling, Sensitivity and
  Monte Carlo Analysis in R Using Package FME.  Journal of Statistical
  Software 33(3) 1--28. \doi{10.18637/jss.v033.i03}
} 
\keyword{utilities}

