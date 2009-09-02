\name{Latinhyper}
\alias{Latinhyper}

\title{
  Latin Hypercube Sampling
}
\description{
  Generates random parameter sets using a latin hypercube sampling algorithm.
}

\usage{
Latinhyper(parRange, num)
}

\arguments{
  \item{parRange }{the range (min, max) of the parameters, a matrix or a
    data.frame with one row for each parameter, and two columns with the
    minimum (1st) and maximum (2nd) column.
  }
  \item{num }{the number of random parameter sets to generate.
  }
}

\value{
  a matrix with one row for each generated parameter set, and one column
  per parameter.
}

\details{
  In the latin hypercube sampling, the space for each parameter is subdivided
  into \code{num} equally-sized segments and one parameter value in each of
  the segments drawn randomly.
}

\note{
  The latin hypercube distributed parameter sets give better coverage in
  parameter space than the uniform random design (\code{\link{Unif}}).
  It is a reasonable choice in case the number of parameter sets is limited.
}

\examples{
## 4 parameters
parRange <- data.frame(min = c(0, 1, 2, 3), max = c(10, 9, 8, 7))
rownames(parRange) <- c("par1", "par2", "par3", "par4")

## Latin hypercube
pairs(Latinhyper(parRange, 100), main = "Latin hypercube")
}
\references{
  Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
  Flannery, B. P. (2007) Numerical Recipes in C. Cambridge
  University Press.
}
\seealso{
  \code{\link{Norm}} for (multi)normally distributed random parameter sets.

  \code{\link{Unif}} for uniformly distributed random parameter sets.

  \code{\link{Grid}} to generate random parameter sets arranged on a regular
    grid
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

