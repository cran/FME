\name{modFit}
\alias{modFit}
\alias{summary.modFit}
\alias{deviance.modFit}
\alias{coef.modFit}
\alias{residuals.modFit}
\alias{df.residual.modFit}
\alias{plot.modFit}
\alias{print.summary.modFit}
\title{
  Constrained Fitting of a Model to Data
}

\description{
  Fitting a model to data, with lower and/or upper bounds
}

\usage{modFit(f, p, ..., lower = -Inf, upper = Inf,
  method = c("Marq", "Port", "Newton",
           "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
           "Pseudo"), jac = NULL,
  control = list(), hessian = TRUE)

\method{summary}{modFit}(object, cov=TRUE,\dots)

\method{deviance}{modFit}(object,  \dots)

\method{coef}{modFit}(object, \dots)

\method{residuals}{modFit}(object, \dots)

\method{df.residual}{modFit}(object, \dots)

\method{plot}{modFit}(x, ask = NULL, \dots)

\method{print}{summary.modFit}(x, digits = max(3, getOption("digits") - 3),
               \dots)
}
\arguments{
  \item{f }{a function to be minimized, with first argument the vector of
    parameters over which minimization is to take place. It should return
    either a vector of \emph{residuals} (of model versus data) or an element
    of class \emph{modCost} (as returned by a call to \code{\link{modCost}}.
  }
  \item{p }{initial values for the parameters to be optimized over.
  }
  \item{...}{additional arguments passed to function \code{f} (modFit) or passed
    to the methods.
  }
  \item{lower }{lower bounds on the parameters; if unbounded set equal to
    \code{-Inf}.
  }
  \item{upper }{upper bounds on the parameters; if unbounded set equal to
    \code{Inf}.
  }
  \item{method }{The method to be used, one of "Marq","Port","Newton",
    "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Pseudo" - see details.
  }
  \item{jac }{A function that calculates the Jacobian; it should be called as
    jac(x,...) and return the matrix with derivatives of the model residuals
    as a function of the parameters. Supplying the Jacobian can substantially
    improve performance; see last example.
  }
  \item{hessian }{\code{TRUE} if Hessian is to be estimated. Note that, if set to
    \code{FALSE}, then a summary cannot be estimated.
  }
  \item{control }{additional control arguments passed to the optimisation
    routine - see details of \code{\link[minpack.lm]{nls.lm}} ("Marq"),
    \code{\link{nlminb}} ("Port"),
    \code{\link{optim}} ("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
    \code{\link{nlm}} ("Newton") or \code{\link{pseudoOptim}}("Pseudo").
  }
  \item{object }{an object of class \code{modFit}.
  }
  \item{x }{an object of class \code{modFit}.
  }
  \item{digits }{number of significant digits in printout.
  }
  \item{cov }{when \code{TRUE} also calculates the parameter covariances.
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.
  }
}

\value{
  a list of class \emph{modFit} containing the results as returned from the
  called optimisation routines.
  
  This includes the following:
  
  \item{par }{the best set of parameters found.
  }
  \item{ssr }{the sum of squared residuals, evaluated at the best set of
    parameters.
  }
  \item{Hessian }{A symmetric matrix giving an estimate of the Hessian
    at the solution found - see note.
  }
  \item{residuals }{the result of the last \code{f} evaluation; that is, the
    residuals.
  }
  \item{ms }{the mean squared residuals, i.e. \code{ssr/length(residuals)}.
  }
  \item{var_ms }{the weighted and scaled variable mean squared residuals,
    one value per observed variable; only when \code{f}
    returns an element of class \emph{modCost} ; \code{NA} otherwise.
  }
  \item{var_ms_unscaled }{the weighted, but not scaled variable mean
    squared residuals
  }
  \item{var_ms_unweighted }{the raw variable mean squared residuals,
    unscaled and unweighted.
  }
  \item{... }{any other arguments returned by the called optimisation routine.
  }

  Note: this means that some return arguments of the original
  optimisation functions are renamed.
  
  More specifically, "objective" and "counts" from routine  \code{\link{nlminb}}
  (method=Port) are renamed; "value" and "counts"; "niter" and "minimum"
   from routine  \code{\link[minpack.lm]{nls.lm}} (method=Marq)
   are renamed; "counts" and "value"; "minimum" and "estimate" from routine
   \code{\link{nlm}} (\code{method = Newton})
   are renamed.

  The list returned by \code{modFit} has a method for the generic functions
  \code{\link{summary}},\code{\link{deviance}},\code{\link{coef}},
  \code{\link{residuals}},\code{\link{df.residual}} and \code{print.summary}
  - see note.
}

\details{

  Note that arguments after \code{...} must be matched exactly.

  The method to be used is specified by argument \code{method} which can
  be one of the methods from function \code{\link{optim}}:

  \itemize{
    \item "Nelder-Mead", the default from \code{optim},
    \item "BFGS", a quasi-Newton method,
    \item "CG", a conjugate-gradient method,
    \item "L-BFGS-B", constrained quasi-Newton method,
    \item "SANN", method of simulated annealing.
  }
  Or one of the following:
  \itemize{
    \item "Marq", the Levenberg-Marquardt algorithm
     (\code{\link[minpack.lm]{nls.lm}} from package \code{minpack}) -
     the default. Note that this method is the only least squares method.
    \item "Newton", a Newton-type algorithm (see \code{\link{nlm}}),
    \item "Port", the Port algorithm (see \code{\link{nlminb}}),
    \item "Pseudo", a pseudorandom-search algorithm
      (see \code{\link{pseudoOptim}}).
  }

  For difficult problems it may be efficient to perform some iterations with
  \code{Pseudo}, which will bring the algorithm near the vicinity of a
  (the) minimum, after which the default algorithm (\code{Marq}) is used to
  locate the minimum more precisely.

  The implementation for the routines from \code{\link{optim}} differs from
  \code{\link{constrOptim}} which implements an adaptive barrier algorithm and
  which allows a more flexible implementation of linear constraints.

  In \code{modFit}, bounds on parameters are imposed by a transformation of the
  parameters to be fitted.

  In case \emph{both lower and upper bounds} are specified, this is achieved by
  a tangens and arc tangens transformation.

  This is, parameter values, \eqn{p'}, generated by the optimisation routine, and
  which are located in the range [-Inf, Inf] are transformed, before they are
  passed to \code{f} as:

   \deqn{p = (upper + lower)/2 + (upper - lower) \cdot \arctan(p')/\pi}{%
        p = 1/2 * (upper + lower) + (upper - lower) * arctan(p')/pi}.

  which maps them into the interval [lower, upper].

  Before the optimisation routine is called, the original parameter values, as
  given by argument \code{p} are mapped from [lower,upper] to [-Inf, Inf] by:

  \deqn{p' = \tan(\pi/2 \cdot (2 p - upper - lower) / (upper - lower))}{%
       p' = tan(pi/2 * (2 * p - upper - lower) / (upper - lower))}


  In case \emph{only lower or upper bounds} are specified, this is achieved
  by a log transformation and a corresponding exponential back transformation.

  In case parameters are transformed (all methods) or in case the
  \code{method} \code{Port}, \code{Pseudo} or \code{Marq} is selected,
  the \emph{Hessian} is approximated as \eqn{2*J^T*J}, where J is the Jacobian,
  estimated by finite differences.

  This ignores the second derivative terms, but this is reasonable if the method
  has truly converged to the minimum.

  Note that finite differences are not extremely precise.

  In case the Levenberg-Marquard method (\code{Marq}) is used, and parameters
  are not transformed, 0.5 times the Hessian of the least squares problem is
  returned by \code{nls.lm}, the original marquardt algorithm. To make it
  compatible, this value is multiplied with 2 and the TRUE Hessian is thus
  returned by \code{modFit}.
}

\note{
  The \code{\link{summary}} \emph{method} is based on an estimate of the
  parameter covariance matrix.
  In computing the covariance matrix of the fitted parameters, the problem is
  treated as if it were a linear least squares problem, linearizing around
  the parameter values that minimize \eqn{Chi^2}.

  The covariance matrix is estimated as \eqn{1/(0.5*Hessian)}.

  This computation relies on several things, i.e.:

  \enumerate{
    \item the parameter values are located at the minimum (i.e. the fitting
      algorithm has converged).
    \item the observations \eqn{y_j} are subject to independent errors whose
      variances are well estimated by \eqn{1 / (n - p)} times the residual sum
      of squares (where n = number of data points, p = number of parameters).
    \item the model is not too nonlinear.
  }
  This means that the estimated covariance (correlation) matrix and the
  confidence intervals derived from it may be worthless if the assumptions
  behind the covariance computation are invalid.

  If in doubt about the validity of the summary computations, use Monte Carlo
  fitting instead, or run a \code{\link{modMCMC}}.

  Other methods included are:
  \itemize{
    \item \code{\link{deviance}}, which returns the model deviance,
    \item \code{\link{coef}}, which extracts the values of the fitted parameters,
    \item \code{\link{residuals}},which extracts the model residuals,
    \item \code{\link{df.residual}} which returns the residual degrees of freedom
    \item \code{print.summary}, producing a nice printout of the summary.
  }
  
  Specifying a function to estimate the Jacobian matrix via argument
  \code{jac} may increase speed. The Jacobian is used in the methods
  "Marq", "BFGS", "CG", "L-BFGS", "Port", and is also used at the end,
  to estimate the Hessian at the optimum.

  Specification of the \code{gradient} in routines "BFGS", "CG", "L-BFGS"
  from \code{optim} and "port" from \code{nlminb} is not allowed here.
  Within \code{modFit}, the gradient is rather estimated from the Jacobian
  \code{jac} and the function \code{f}.

}
\seealso{
  \code{\link{constrOptim}} for constrained optimization.
}

\examples{

## =======================================================================
## logistic growth model
## =======================================================================
TT    <- seq(1, 60, 5)
N0    <- 0.1
r     <- 0.5
K     <- 100

## perturbed analytical solution
Data <- data.frame(
          time = TT,
          N = K / (1+(K/N0-1) * exp(-r*TT)) * (1 + rnorm(length(TT), sd = 0.01))
        )

plot(TT, Data[,"N"], ylim=c(0, 120), pch=16, col="red",
     main = "logistic growth", xlab="time", ylab="N")


##===================================
## Fitted with analytical solution  #
##===================================

## initial "guess"
parms <- c(r = 2, K = 10, N0 = 5)

## analytical solution
model <- function(parms,time)
  with (as.list(parms), return(K/(1+(K/N0-1)*exp(-r*time))))

## run the model with initial guess and plot results
lines(TT, model(parms, TT), lwd = 2,col = "green")

## FITTING algorithm 1
ModelCost <- function(P) {
 out <- model(P,TT)
 return(Data$N-out)  # residuals
}

(Fita <- modFit(f = ModelCost, p = parms))

times <- 0:60
lines(times, model(Fita$par, times), lwd = 2, col = "blue")
summary(Fita)

##===================================
##  Fitted with numerical solution  #
##===================================

## numeric solution
logist <- function(t, x, parms) {
  with(as.list(parms), {
    dx <- r * x[1] * (1 - x[1]/K)
    list(dx)
  })
}

## model cost,
ModelCost2 <- function(P) {
 out <- ode(y = c(N = P[["N0"]]), func = logist, parms = P, times = c(0, TT))
 return(modCost(out, Data))  # object of class modCost
}

Fit <- modFit(f=ModelCost2,p=parms,lower=rep(0,3),upper=c(5,150,10))

out <- ode(y=c(N=Fit$par[["N0"]]),func=logist,parms=Fit$par,times=times)

lines(out,col="red",lty=2)
legend("right",c("data","original","fitted analytical","fitted numerical"),
       lty=c(NA,1,1,2),lwd=c(NA,2,2,1),col=c("red","green","blue","red"),
       pch=c(16,NA,NA,NA))
summary(Fit)
plot(residuals(Fit))

## =======================================================================
## the use of the Jacobian
## =======================================================================

## We use modFit to solve a set of linear equations
A <- matrix(nr=30,nc=30,runif(900))
X <- runif(30)
B <- A \%*\% X

## try to find vector "X"; the Jacobian is matrix A

## Function that returns the vector of residuals
FUN <- function(x)
  as.vector(A \%*\%x - B)

## Function that returns the Jacobian
JAC <- function(x) A

## The port algorithm
print(system.time(mf <- modFit(f=FUN,p=runif(30),method="Port")))
print(system.time(mf <- modFit(f=FUN,p=runif(30),method="Port",jac=JAC)))
max(abs(mf$par-X))  # should be very small

## BFGS
print(system.time(mf <- modFit(f=FUN,p=runif(30),method="BFGS")))
print(system.time(mf <- modFit(f=FUN,p=runif(30),method="BFGS",jac=JAC)))
max(abs(mf$par-X))

## Levenberg-Marquardt
print(system.time(mf <- modFit(f=FUN,p=runif(30),jac=JAC)))
max(abs(mf$par-X))
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,

  Thomas Petzoldt <thomas.petzoldt@tu-dresden.de>
}
\references{
  Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
  Flannery, B. P. (2007) Numerical Recipes in C. Cambridge
  University Press.

  \url{http://netlib.bell-labs.com/cm/cs/cstr/153.pdf}
}
\keyword{utilities}

