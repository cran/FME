\name{pseudoOptim}
\alias{pseudoOptim}
\title{Pseudo-random Search Optimisation Algorithm of Price (1977)}
\description{Fits a model to data, using the pseudo-random search
  algorithm of Price (1977), a random-based fitting technique.
}

\usage{pseudoOptim(f, p,..., lower, upper, control = list())}
\arguments{
  \item{f }{function to be minimised, its first argument should be the vector
    of parameters over which minimization is to take place. It should return
    a scalar result, the model cost, e.g the sum of squared residuals.}
  \item{p }{initial values of the parameters to be optimised.}
  \item{... }{arguments passed to funtion \code{f}.}
  \item{lower }{minimal values of the parameters to be optimised; these
    must be specified; they cannot be -Inf.}
  \item{upper }{maximal values of the parameters to be optimised; these
    must be specified; they cannot be +Inf.}
  \item{control}{a list of control parameters - see details.}
}
\value{
  a list containing:
  \item{par }{the optimised parameter values.}
  \item{cost }{the model cost, or function evaluation associated to the
    optimised parameter values, i.e. the minimal cost.}
  \item{iterations}{the number of iterations performed.}
    and if control\$verbose is TRUE:
  \item{poppar }{all parameter vectors remaining in the population, matrix
    of dimension (npop,length(par)).}
  \item{popcost }{model costs associated with all population parameter vectors,
    vector of length npop.}
  \item{rsstrace}{a 2-columned matrix with the iteration number and the model
    cost at each succesful iteration.}
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}
\examples{
amp    <- 6
period <- 5
phase  <- 0.5

x <- runif(20)*13 
y <- amp*sin(2*pi*x/period+phase) + rnorm(20, mean = 0, sd = 0.05)
plot(x, y, pch = 16)


cost <- function(par)
    sum((par[1] * sin(2*pi*x/par[2]+par[3])-y)^2)

p1 <- optim(par = c(amplitude = 1, phase = 1, period = 1), fn = cost)
p2 <- optim(par = c(amplitude = 1, phase = 1, period = 1), fn = cost, 
            method = "SANN")
p3 <- pseudoOptim(p = c(amplitude = 1, phase = 1, period = 1), 
            lower = c(0, 1e-8, 0), upper = c(100, 2*pi, 100), 
            f = cost, control = c(numiter = 3000, verbose = TRUE))

curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]), lty = 2, add = TRUE)
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]), lty = 3, add = TRUE)
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]), lty = 1, add = TRUE)
legend ("bottomright", lty = c(1, 2, 3),
         c("Price", "Mathematical", "Simulated annealing"))}
\details{
The \code{control} argument is a list that can supply any of the
following components:
\itemize{
  \item npop, number of elements in the population.
    Defaults to max(5*length(p),50).
  \item numiter, maximal number of iterations to be performed.
    Defaults to 10000. The algorithm either stops when \code{numiter}
    iterations has been performed or when the remaining variation is less
    than \code{varleft}.
  \item centroid, number of elements from which to estimate a new parameter
    vector, defaults to 3.
  \item varleft, relative variation remaining; if below this value the
    algorithm stops; defaults to 1e-8.
  \item verbose, if TRUE, more verbose output will contain the parameters
    in the final population, their respective population costs and the cost
    at each succesful interation. Defaults to \code{FALSE}.
}
see the book of Soetaert and Herman (2009) for a description of the algorithm
AND for a line to line explanation of the function code.
}
\references{
  Soetaert, K. and Herman, P. M. J., 2009.  A Practical Guide to Ecological
  Modelling. Using R as a Simulation Platform.  Springer, 372 pp.
 
  Price, W.L., 1977.  A Controlled Random Search Procedure for Global
  Optimisation.  The Computer Journal, 20: 367-370.
}
\keyword{optimize}

