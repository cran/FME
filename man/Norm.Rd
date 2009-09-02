\name{Norm}
\alias{Norm}

\title{
  Normal Random Distribution
}
\description{
  Generates random parameter sets that are (multi)normally distributed.
}
\usage{
Norm(parMean, parCovar, parRange = NULL, num)
}
\arguments{
  \item{parMean }{a vector, with the mean value of each parameter.
  }
  \item{parCovar }{the parameter variance-covariance matrix.
  }
  \item{parRange }{the range (min, max) of the parameters, a matrix or a
    data.frame with one row for each parameter, and two columns with the
    minimum (1st) and maximum (2nd) column.
  }
  \item{num }{the number of random parameter sets to generate.
  }
}
\value{
  a matrix with one row for each generated parameter set, and one column per
  parameter.
}

\details{
  Function \code{Norm}, draws parameter sets from a multivariate normal
  distribution, as specified through the mean value and the variance-covariance
  matrix of the parameters. In addition, it is possible to impose a minimum and
  maximum of each parameter, via \code{parRange}. This will generate a
  truncated distribution.
  Use this for instance if certain parameters cannot become negative.
}

\note{
  For function \code{Norm} to work, \code{parCovar} must be a valid
  variance-covariance matrix. (i.e. positive definite). If this is not the
  case, then the function will fail.
}
\examples{
## multinormal parameters: variance-covariance matrix and parameter mean
parCovar <- matrix(data = c(0.5, -0.2, 0.3, 0.4, -0.2, 1.0, 0.1, 0.3,
                   0.3, 0.1, 1.5, -0.7, 1.0, 0.3, -0.7, 4.5), nrow = 4)
parCovar

parMean <- 4:1

## Generated sample
Ndist <- Norm(parCovar = parCovar, parMean = parMean, num = 500)
cov(Ndist)   # check
pairs(Ndist, main = "normal")

## truncated multinormal
Ranges <- data.frame(min = rep(0, 4), max = rep(Inf, 4))

pairs(Norm(parCovar = parCovar, parMean = parMean, parRange = Ranges,
      num = 500), main = "truncated normal")

}
\seealso{
  \code{\link{Unif}} for uniformly distributed random parameter sets.

  \code{\link{Latinhyper}} to generates parameter sets using
    latin hypercube sampling.

  \code{\link{Grid}} to generate random parameter sets arranged on a regular
    grid

  \code{\link{rnorm}} the R-default for generating normally distributed random
    numbers.

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

