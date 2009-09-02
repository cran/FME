\name{Unif}
\alias{Unif}

\title{
  Uniform Random Distribution
}
\description{
  Generates uniformly distributed random parameter sets.
}

\usage{
Unif(parRange, num)
}

\arguments{
  \item{parRange }{the range (min, max) of the parameters, a matrix or a
    data.frame with one row for each parameter, and two columns with the
    minimum (1st) and maximum (2nd) value.
  }
  \item{num }{the number of random parameter sets to generate.
  }
}

\value{
  a matrix with one row for each generated parameter set, and one column
  per parameter.
}

\details{
  In the uniform sampling, each parameter is uniformly random distributed
  over its range.
}

\note{
  For small sample sizes, the latin hypercube distributed parameter sets
  (\code{\link{Latinhyper}}) may give better coverage in parameter space
  than the uniform random design.
}

\examples{
## 4 parameters
parRange <- data.frame(min = c(0, 1, 2, 3), max = c(10, 9, 8, 7))
rownames(parRange) <- c("par1", "par2", "par3", "par4")

## uniform
pairs(Unif(parRange, 100), main = "Uniformly random")
}

\seealso{
  \code{\link{Norm}} for (multi)normally distributed random parameter sets.

  \code{\link{Latinhyper}} to generates parameter sets using
    latin hypercube sampling.

  \code{\link{Grid}} to generate random parameter sets arranged on a
    regular grid

  \code{\link{runif}} the R-default for generating uniformally distributed
    random numbers.

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

