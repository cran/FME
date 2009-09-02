\name{Grid}
\alias{Grid}

\title{
  Grid Distribution
}

\description{
  Generates parameter sets arranged on a regular grid.
}

\usage{
Grid(parRange, num)
}
\arguments{
  \item{parRange }{the range (min, max) of the parameters, a matrix or a
    data.frame with one row for each parameter, and two columns with the minimum
    (1st) and maximum (2nd) column.
  }
  \item{num }{the number of random parameter sets to generate.
  }
}
\value{
  a matrix with one row for each generated parameter set, and one column per
  parameter.
}

\details{
  The grid design produces the most regular parameter distribution; there is no
  randomness involved.
  The number of parameter sets generated with \code{Grid} will be <= \code{num}.
}

\examples{

## 4 parameters
parRange <- data.frame(min = c(0, 1, 2, 3), max = c(10, 9, 8, 7))
rownames(parRange) <- c("par1", "par2", "par3", "par4")

## grid
pairs(Grid(parRange, 500), main = "Grid")

}
\seealso{
  \code{\link{Norm}} for (multi)normally distributed random parameter sets.

  \code{\link{Latinhyper}} to generates parameter sets using
    latin hypercube sampling.

  \code{\link{Unif}} for uniformly distributed random parameter sets.
  
  \code{\link{seq}} the R-default for generating regular sequences of numbers.

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

