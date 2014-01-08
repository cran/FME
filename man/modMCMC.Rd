\name{modMCMC}
\alias{modMCMC}
\alias{summary.modMCMC}
\alias{plot.modMCMC}
\alias{pairs.modMCMC}
\alias{hist.modMCMC}
%%\alias(cumuplot.modMCMC)
\title{
  Constrained Markov Chain Monte Carlo}
\description{
  Performs a Markov Chain Monte Carlo simulation, using an adaptive
  Metropolis (AM) algorithm and including a delayed rejection (DR)
  procedure.
}

\usage{
modMCMC(f, p, ..., jump = NULL,  lower = -Inf, upper = +Inf, 
        prior = NULL, var0 = NULL, wvar0 = NULL, n0 = NULL, 
        niter = 1000, outputlength = niter, burninlength = 0, 
        updatecov = niter, covscale = 2.4^2/length(p),
        ntrydr = 1, drscale = NULL, verbose = TRUE)

%%cumuplot.modMCMC(x, Full=FALSE, which=1:ncol(x$pars), ...)

\method{summary}{modMCMC}(object, remove = NULL, \dots)

\method{pairs}{modMCMC}(x, Full = FALSE, which = 1:ncol(x$pars),
  remove = NULL, nsample = NULL, \dots)

\method{hist}{modMCMC}(x, Full = FALSE, which = 1:ncol(x$pars),
  remove = NULL, ask = NULL, \dots)

\method{plot}{modMCMC}(x, Full = FALSE, which = 1:ncol(x$pars),
  trace = TRUE, remove = NULL, ask = NULL, \dots)
}

\arguments{
  \item{f }{the function to be evaluated, with first argument the vector
    of parameters which should be varied. It should return either the
    model residuals, an element of class \emph{modCost} (as returned by
    a call to \code{\link{modCost}}) or -2*log(likelihood).  The latter
    is equivalent to the sum-of-squares functions when using a Gaussian
    likelihood and prior.
  }
  \item{p }{initial values for the parameters to be optimized over.
  }
  \item{...}{additional arguments passed to function \code{f} or to the
    methods.
  }
  \item{jump }{jump length, either a \emph{number}, a \emph{vector} with
    length equal to the total number of parameters, a \emph{covariance
    matrix}, or a \emph{function} that takes as input the current values
    of the parameters and produces as output the perturbed
    parameters. See details.
  }
  \item{prior }{-2*log(parameter prior probability), either a function
    that is called as \code{prior(p)} or \code{NULL}; in the latter case
    a non-informative prior is used (i.e. all parameters are equally
    likely, depending on \code{lower} and \code{upper} within min and max bounds).
  }
  \item{var0 }{initial model variance; if \code{NULL}, it is assumed
    that the model variance is \code{1}, and the return element from \code{f}
    is -2*log (likelihood). If it has a value, it is assumed that the return element
    from \code{f} contain the model residuals or a list of class
    \code{modFit}.  See details. Good options for \code{var0} are to use 
    the modelvariance (\code{modVariance}) as returned by the 
    \code{summary} method of \link{modFit}. When this option is chosen, and 
    the model has several variables, they will all be scaled similarly.
    See vignette \code{FMEdyna}. In case the model has several variables with 
    different magnitudes, then it may be better to scale each variable
    independently. In that case, one can use as \code{var0}, the mean of the 
    unweighted squared residuals from the model fit as returned from 
    \link{modFit} (\code{var_ms_unweighted}). See vignette \code{FME}.
  }
  \item{wvar0 }{"weight" for the initial model variance -- see details.
  }
  \item{n0 }{parameter used for weighing the initial model variance -
    if \code{NULL}, it is estimated as \code{n0=wvar0*n}, where n = number
    of observations. See details.
  }
  \item{lower }{lower bounds on the parameters; for unbounded parameters
    set equal to \code{-Inf}.
  }
  \item{upper }{upper bounds on the parameters; for unbounded parameters
    set equal to \code{Inf}.
  }
  \item{niter }{number of iterations for the MCMC.
  }
  \item{outputlength }{number of iterations kept in the output; should
    be smaller or equal to \code{niter}.
  }
  \item{updatecov }{number of iterations after which the parameter
    covariance matrix is (re)evaluated based on the parameters kept thus
    far, and used to update the MCMC jumps.
  }
  \item{covscale }{scale factor for the parameter covariance matrix,
    used to perform the MCMC jumps.
  }
  \item{burninlength }{number of initial iterations to be removed from
    output.
  }
  \item{ntrydr }{maximal number of tries for the delayed rejection
    procedure.  It is generally not a good idea to set this to a too
    large value.
  }
  \item{drscale }{for each try during delayed rejection, the cholesky
    decomposition of the proposal matrix is scaled with this amount; if
    \code{NULL}, it is assumed to be \code{c(0.2,0.25, 0.333, 0.333,
    ...)}  }
  \item{verbose }{if \code{TRUE}: prints extra output.
  }
  \item{object }{an object of class \code{\link{modMCMC}}.
  }
  \item{x }{an object of class \code{\link{modMCMC}}.
  }
  \item{Full }{If TRUE then not only the parameters will be plotted, but
    also the function value and (if appropriate) the model variance(s).
  }
  \item{which }{the name or the index to the parameters that should be
    plotted.  Default = all parameters. If \code{Full=TRUE}, setting
    \code{which =NULL} will plot only the function value and the model
    variance.
  }
  \item{trace }{if \code{TRUE}, adds smoothed line to the plot.
  }
  
  \item{remove }{a list with indices of the runs that should be removed (e.g.
    to remove runs during burnin.
  }
  \item{nsample }{the number of xy pairs to be plotted on the upper
    panel in the pairs plot. When \code{NULL} all xy pairs plotted. Set
    to a lower number in case the graph becomes too dense (and the
    exported picture too large).  This does not affect the histograms on
    the diagonal plot (which are estimated using all MCMC draws).
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.
  }
}
\value{
  a list of class \emph{modMCMC} containing the results as returned from the
  Markov chain.
  
  This includes the following:
  
  \item{pars }{an array with dimension (\code{outputlength},
    length(\code{p})), containing the parameters of the MCMC at each
    iteration that is kept.
  }
  \item{SS }{vector with the sum of squares function, one for each row
    in \code{pars}.}
  \item{naccepted }{the number of accepted runs.
  }
  \item{sig }{the sampled error variance \eqn{\sigma^2}{sigma^2}, a matrix with
    one row for each row in \code{pars}.
  }
  \item{bestpar }{the parameter set that gave the highest probability.
  }
  \item{bestfunp }{the function value corresponding to \code{bestpar}.
  }
  \item{prior }{the parameter prior, one value for each row in
    \code{pars}.
  }
  \item{count }{information about the MCMC chain: number of delayed
    rejection steps (\code{dr_steps}), the number of alfa steps
    \code{Alfasteps}, the number of accepted runs (\code{num_accepted})
    and the number of times the proposal covariance matrix has been
    updated (\code{num_covupdate}.)
  }
  \item{settings }{the settings for error covariance calculation,
    i.e. arguments \code{var0}, \code{n0} and \code{N} the number of
    data points.
  }
  The list returned by \code{modMCMC} has methods for the generic
  functions \code{\link{summary}}, \code{\link{plot}},
  \code{\link{pairs}} -- see note.
}

\details{
  Note that arguments after ... must be matched exactly.

  R-function \code{f} is called as \code{f(p, ...)}. It should return
  either -2 times the log likelihood of the model (one value), the
  residuals between model and data or an item of class \code{modFit} (as
  created by function \code{\link{modFit}}.

  In the latter two cases, it is assumed that the prior distribution for
  \eqn{\theta}{theta} is either non-informative or gaussian. If gaussian, it
  can be treated as a sum of squares (SS). If the measurement function
  is defined as:
  
\deqn{y=f(\theta) + \xi\\
\xi ~ N(0,\sigma^2)}{y = f(p) + xi; xi ~ N(0, sigma^2)}

  where \eqn{\xi}{xi} is the measurement error, assumed normally
  distribution, then the posterior for the parameters will be estimated
  as:
  
  \deqn{p(\theta | y,\sigma^2)\propto exp(-0.5 \cdot (\frac{SS(\theta)}{\sigma^2}
  +SS_{pri}(\theta))}{prob(p|y,sigma^2)~exp(-0.5 SS(p)/sigma^2)+SSpri(p)}

  and where \eqn{\sigma^2}{sigma^2} is the error variance, SS is the sum of squares
  function \eqn{SS(\theta)=\sum(y_i-f(\theta))^2}{SS(p)=sum(yi-f(p))^2}.
  If non-informative priors are used, then \eqn{SS_{pri}(\theta)=0}{SSpri(p)=0}.

  The error variance \eqn{\sigma^2}{sigma^2} is considered a nuisance parameter.
  A prior distribution of it should be specified and a posterior distribution
  is estimated.

  If \code{wvar0} or \code{n0} is >0, then the variances are sampled as
  conjugate priors from the inverse gamma distribution with parameters
  \code{var0} and \code{n0=wvar0*n}.  Larger values of \code{wvar0} keep
  these samples closer to \code{var0}.

  Thus, at each step, 1/ the error variance (\eqn{\sigma^{-2}}{sigma^{-2}}) is sampled
  from a gamma distribution:
  
  \deqn{p(\sigma^{-2}|y,\theta) \sim \Gamma(\frac{(n_0+n)}{2},
  \frac{(n_0 \cdot var0+SS(\theta))}{2})}{prob(sigma^2|y,p)~Gam(0.5*(n0+n),
  0.5*(n0+var0+SS(p))}

  where \code{n} is the number of data points and where
  \eqn{n0=n \cdot wvar0}{n0=n * wvar0}, and where the second argument to the gamma function
  is the shape parameter.

  The prior parameters (\code{var0} and \code{wvar0}) are the prior mean
  for \eqn{\sigma^2}{sigma^2} and the prior accuracy.

  By setting \code{wvar0} equal to 1, equal weight is given to the
  prior and the current value.

  If \code{wvar0} is 0 then the prior is ignored.

  If \code{wvar0} is \code{NULL} (the default) then the error variances are
  assumed to be fixed.

  \code{var0} estimates the variance of the measured components. In case
  independent estimates are not available, these variances can be
  obtained from the mean squares of fitted residuals. (e.g. as reported
  in \code{modFit}). See the examples. (but note that this is not truly
  independent information)

  \code{var0} is either one value, or a value for each observed
  variable, or a value for each observed data point.

  When \code{var0} is not \code{NULL}, then \code{f} is assumed to
  return the model residuals OR an instance of class \code{modCost}.

  When \code{var0=NULL}, then \code{f} should return either
  -2*log(probability of the model), or an instance of class
  \code{modCost}.

  \code{modMCMC} implements the Metropolis-Hastings method. The proposal
  distribution, which is used to generate new parameter values is the
  (multidimensional) Gaussian density distribution, with standard
  deviation given by \code{jump}.

  \code{jump} can be either one value, a vector of length = number of
  parameters or a parameter covariance matrix (nrow = ncol = number
  parameters).

  The jump parameter, \code{jump} thus determines how much the new
  parameter set will deviate from the old one.

  If \code{jump} is one value, or a vector, then the new parameter
  values are generated by sampling a normal distribution with standard
  deviation equal to \code{jump}.  A larger value will lead to larger
  jumps in the parameter space, but acceptance of new points can get
  very low. Smaller jump lengths increase the acceptance rate, but the
  algorithm may move too slowly, and too many runs may be needed to scan
  the parameter space.

  If \code{jump} is \code{NULL}, then the jump length is taken as 10\%
  of the parameter value as given in \code{p}.

  \code{jump} can also be a proposal covariance matrix. In this case,
  the new parameter values are generated by sampling a multidimensional
  normal distribution. It can be efficient to initialise \code{jump}
  using the parameter covariance as resulting from fitting the model
  (e.g. using \code{modFit}) -- see examples.

  Finally, \code{jump} can also be an R-function that takes as input the
  current values of the parameters and returns the new parameter values.

  Two methods are implemented to increase the number of accepted runs.
  \enumerate{
    \item In the \emph{"adaptive Metropolis"} method, new parameters are
      generated with a covariance matrix that is estimated from the
      parameters generated (and saved) thus far.  The idea behind this
      is that the MCMC method is more efficient if the proposal
      covariance (to generate new parameter values) is somehow tuned to
      the shape and size of the target distribution.

      Setting \code{updatecov} smaller than \code{niter} will trigger
      this functionality. In this case, every \code{updatecov}
      iterations, the jump covariance matrix will be estimated from the
      covariance matrix of the saved parameter values. The covariance
      matrix is scaled with \eqn{(2.4^2/npar)} where npar is the number
      of parameters, unless \code{covscale} has been given a different
      value.  Thus, \eqn{Jump= ( cov(\theta_1,\theta_2,....\theta_n)
      \cdot diag(np,+1e^{-16})\cdot(2.4^2/npar)}{Jump=cov(p1,p2,...pn) *
      diag(np+1e-16)*2.4^2/npar} where the small number \eqn{1e^{-16}}
      is added on the diagonal of the covariance matrix to prevent it
      from becoming singular.

      Note that a problem of adapting the proposal distribution using
      the MCMC results so far is that standard convergence results do
      not apply.  One solution is to use adaptation only for the burn-in
      period and discard the part of the chain where adaptation has been
      used.

      Thus, when using \code{updatecov} with a positive value of
      \code{burninlength}, the proposal distribution is only updated
      during burnin. If \code{burninlength} = 0 though, the updates
      occur throughout the entire simulation.

      When using the adaptive Metropolis method, it is best to start
      with a small value of the jump length.

    \item In the \emph{"delayed rejection"} method, new parameter
      values are tried upon rejection. The process of delaying rejection
      can be iterated for at most \code{ntrydr} trials. Setting
      \code{ntrydr} equal to 1 (the default) toggles off delayed
      rejection.

      During the delayed rejection procedure, new parameters are
      generated from the last accepted value by scaling the jump
      covariance matrix with a factor as specified in
      \code{drscale}. The acceptance probability of this new set depends
      on the candidates so far proposed and rejected, in such a way that
      reversibility of the Markov chain is preserved. See Haario et
      al. (2005, 2006) for more details.
    }

  Convergence of the MCMC chain can be checked via \code{plot}, which
  plots for each iteration the values of all parameters, and if
  \code{Full} is \code{TRUE}, of the function value (SS) and (if
  appropriate) the modeled variance.  If converged, there should be no
  visible drift.

  In addition, the methods from package \code{coda} become available by
  making the object returned by \code{modMCMC} of class \code{mcmc}, as
  used in the methods of \code{coda}.  For instance, if object
  \code{MCMCres} is returned by \code{modMCMC} then
  \code{as.mcmc(MCMCres$pars)} will make an instance of class
  \code{mcmc}, usable by \code{coda}.

  The \code{burninlength} is the number of initial steps that is not
  included in the output. It can be useful if the initial value of the
  parameters is far from the optimal value. Starting the MCMC with the
  best fit parameter set will alleviate the need for using
  \code{burninlength}.
}

\note{
  The following S3 methods are provided:
  \itemize{
    \item summary, produces summary statistics of the MCMC results
    \item plot, plots the MCMC results, for all parameters. Use it to
      check convergence.
    \item pairs, produces a pairs plot of the MCMC results; overrides
      the default \code{gap = 0}, \code{upper.panel = NA}, and
      \code{diag.panel}.
  }

%%\code{cumuplot.modMCMC} invokes \code{\link[coda]{cumuplot}} from the \code{coda} package.

  It is also possible to use the methods from the \code{coda} package,
  e.g.  \code{\link[coda]{densplot}}.

  To do that, first the \code{modMCMC} object has to be converted to an
  mcmc object.  See the examples for an application.
}

\seealso{
  \code{\link{modFit}} for constrained model fitting
}

\examples{

## =======================================================================
## Sampling a 3-dimensional normal distribution,
## =======================================================================
# mean = 1:3, sd = 0.1
# f returns -2*log(probability) of the parameter values

NN <- function(p) {
  mu <- c(1,2,3)
  -2*sum(log(dnorm(p, mean = mu, sd = 0.1)))   #-2*log(probability)
}

# simple Metropolis-Hastings
MCMC <- modMCMC(f = NN, p = 0:2, niter = 5000,
                outputlength = 1000, jump = 0.5)

# More accepted values by updating the jump covariance matrix...
MCMC <- modMCMC(f = NN, p = 0:2, niter = 5000, updatecov = 100,
                 outputlength = 1000, jump = 0.5)
summary(MCMC)

plot(MCMC)   # check convergence
pairs(MCMC)

## =======================================================================
## test 2: sampling a 3-D normal distribution, larger standard deviation...
## noninformative priors, lower and upper bounds imposed on parameters
## =======================================================================

NN <- function(p) {
  mu <- c(1,2,2.5)
  -2*sum(log(dnorm(p, mean = mu, sd = 0.5)))   #-2*log(probability)
}

MCMC2 <- modMCMC(f = NN, p = 1:3, niter = 2000, burninlength = 500,
  updatecov = 10, jump = 0.5, lower = c(0, 2, 1), upper = c(1, 3, 3))
plot(MCMC2)
hist(MCMC2, breaks = 20)

## Compare output of p3 with theoretical distribution
hist(MCMC2, which = "p3", breaks = 20)
lines(seq(1, 3, 0.1), dnorm(seq(1, 3, 0.1), mean = 2.5,
  sd = 0.5)/pnorm(3, 2.5, 0.5))
summary(MCMC2)

# functions from package coda...
cumuplot(as.mcmc(MCMC2$pars))
summary(as.mcmc(MCMC2$pars))
raftery.diag(MCMC2$pars)

## =======================================================================
## test 3: sampling a log-normal distribution, log mean=1:4, log sd = 1
## =======================================================================

NL <- function(p) {
  mu <- 1:4
  -2*sum(log(dlnorm(p, mean = mu, sd = 1)))      #-2*log(probability)
}
MCMCl <- modMCMC(f = NL, p = log(1:4), niter = 3000,
                 outputlength = 1000, jump = 5)
plot(MCMCl)   # bad convergence
cumuplot(as.mcmc(MCMCl$pars))

MCMCl <- modMCMC (f = NL, p = log(1:4), niter = 3000,
                  outputlength = 1000, jump = 2^(2:5))
plot(MCMCl)   # better convergence but CHECK it!
pairs(MCMCl)
colMeans(log(MCMCl$pars))
apply(log(MCMCl$pars), 2, sd)

MCMCl <- modMCMC (f = NL, p = rep(1, 4), niter = 3000, 
                  outputlength = 1000, jump = 5, updatecov = 100)
plot(MCMCl)
colMeans(log(MCMCl$pars))
apply(log(MCMCl$pars), 2, sd)

## =======================================================================
## Fitting a Monod (Michaelis-Menten) function to data
## =======================================================================

# the observations
#---------------------
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
plot(Obs, pch = 16, cex = 2, xlim = c(0, 400), ylim = c(0, 0.15),
     xlab = "mg COD/l", ylab = "1/hr", main = "Monod")

# the Monod model
#---------------------
Model <- function(p, x) data.frame(x = x, N = p[1]*x/(x+p[2]))

# Fitting the model to the data
#---------------------
# define the residual function
Residuals  <- function(p) (Obs$y - Model(p, Obs$x)$N)

# use modFit to find parameters
P      <- modFit(f = Residuals, p = c(0.1, 1))

# plot best-fit model
x      <-0:375
lines(Model(P$par, x))

# summary of fit
sP    <- summary(P)
sP[]
print(sP)

# Running an MCMC
#---------------------
# estimate parameter covariances
# (to efficiently generate new parameter values)
Covar   <- sP$cov.scaled * 2.4^2/2

# the model variance
s2prior <- sP$modVariance

# set nprior = 0 to avoid updating model variance
MCMC <- modMCMC(f = Residuals, p = P$par,jump = Covar, niter = 1000,
                var0 = s2prior, wvar0 = NULL, updatecov = 100)

plot(MCMC, Full = TRUE)
pairs(MCMC)
# function from the coda package.
raftery.diag(as.mcmc(MCMC$pars))
cor(MCMC$pars)

cov(MCMC$pars)   # covariances by MCMC
sP$cov.scaled    # covariances by Hessian of fit

x  <- 1:400
SR <- summary(sensRange(parInput = MCMC$pars, func = Model, x = x))
plot(SR, xlab="mg COD/l", ylab = "1/hr", main = "Monod")
points(Obs, pch = 16, cex = 1.5)

}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  
  Marko Laine <Marko.Laine@fmi.fi>
}
\references{
  Laine, M., 2008. Adaptive MCMC Methods With Applications in
  Environmental and Geophysical Models. Finnish Meteorological Institute
  contributions 69, ISBN 978-951-697-662-7, Finnish Meteorological
  Institute, Helsinki.

  Haario, H., Saksman, E. and Tamminen, J., 2001. An Adaptive Metropolis
  Algorithm. Bernoulli 7, pp. 223--242.

  Haario, H., Laine, M., Mira, A. and Saksman, E., 2006. DRAM: Efficient
  Adaptive MCMC. Statistics and Computing, 16(4), 339--354.
     
  Haario, H., Saksman, E. and Tamminen, J., 2005. Componentwise
  Adaptation for High Dimensional MCMC. Computational Statistics 20(2),
  265--274.

  Gelman, A. Varlin, J. B., Stern, H. S. and Rubin, D. B.,
  2004. Bayesian Data Analysis.  Second edition. Chapman and Hall, Boca
  Raton.

  Soetaert, K. and Petzoldt, T., 2010.  Inverse Modelling, Sensitivity
  and Monte Carlo Analysis in R Using Package FME.  Journal of
  Statistical Software 33(3)
  1--28. \url{http://www.jstatsoft.org/v33/i03}

}
\keyword{utilities}

