\name{cross2long}
\alias{cross2long}

\title{
  Convert a dataset in wide (crosstab) format to long (database) format
}

\description{
  Rearranges a data frame in cross tab format by putting all relevant columns below
  each other, replicating the independent variable and, if necessary, other specified
  columns. Optionally, an err column is added.
}

\usage{
cross2long( data, x, select = NULL, replicate = NULL, 
            error = FALSE,  na.rm = FALSE)
}

\arguments{

  \item{data }{a data frame (or matrix) with crosstab layout
  }
  \item{x }{name of the independent variable to be replicated
  }
  \item{select }{a vector of column names to be included (see details).
     All columns are included if not specified.
  }
  \item{replicate }{a vector of names of variables (apart from the
     independent variable that have to be replicated for every included
     column (e.g. experimental treatment specification).
  }
  \item{error }{boolean indicating whether the final dataset in long format
     should contain an extra column for error values (cf. \link{modCost});
     here filled with 1's.
  }
  \item{na.rm }{whether or not to remove the \code{NA}s.
  }
}

\value{
  A data frame with the following columns:
  
  \item{name }{Column containing the column names of the original crosstab data frame, \code{data}
  }
  \item{x }{A replication of the independent variable
  }
  \item{y }{The actual data stacked upon each other in one column
  }
  \item{err }{Optional column, filled with NA values (necessary for some other functions)
  }
  \item{... }{all other columns from the original dataset that had to be replicated
     (indicated by the parameter \code{replicate})
  }
}

\details{
  The original data frame is converted from a wide (crosstab) layout (one
  variable per column) to a long (database) layout (all variable value in
  one column).

  As an example of both formats consider the data, called \code{Dat} consisting
  of two observed variables, called "Obs1" and "Obs2", both containing two
  observations, at time 1 and 2:
  \tabular{llll}{
  name    \tab time   \tab   val \tab err \cr
  Obs1    \tab 1      \tab   50  \tab 5  \cr
  Obs1    \tab 2      \tab  150  \tab 15  \cr
  Obs2    \tab 1      \tab  1    \tab 0.1  \cr
  Obs2    \tab 2      \tab  2    \tab 0.2  \cr
  }
  for the long format and

  \tabular{lll}{
   time   \tab   Obs1 \tab Obs2 \cr
   1      \tab  50    \tab 1 \cr
   2      \tab  150   \tab 2 \cr
  }
  for the crosstab format.

  The parameters \code{x}, \code{select}, and  \code{replicate} should
  be disjoint. Although the independent variable always has to be replicated
  it should not be given by the \code{replicate} parameter.

}

\examples{

## =======================================================================
## Suppose we have measured sediment oxygen concentration profiles
## =======================================================================

depth  <- 0:7
O2mud  <- c( 6,   1,   0.5, 0.1, 0.05,0,   0,   0)
O2silt <- c( 6,   5,   3,   2,   1.5, 1,   0.5, 0)
O2sand <- c( 6,   6,   5,   4,   3,   2,   1,   0)
zones  <- c("a", "b", "b", "c", "c", "d", "d", "e")
oxygen <- data.frame(depth = depth,
                     zone  = zones,
                     mud   = O2mud,
                     silt  = O2silt,
                     sand  = O2sand
          )

 cross2long(data = oxygen, x = depth, 
            select = c(silt, mud), replicate = zone)

 cross2long(data = oxygen, x = depth, 
            select = c(mud, -silt), replicate = zone)

# twice the same column name: replicates
 colnames(oxygen)[4] <- "mud"    

 cross2long(data=oxygen, x = depth, select = mud)
}
\author{
  Tom Van Engeland <tom.vanengeland@nioz.nl>
}

\references{
  Soetaert, K. and Petzoldt, T., 2010.  Inverse Modelling, Sensitivity and
  Monte Carlo Analysis in R Using Package FME.  Journal of Statistical
  Software 33(3) 1--28.
  \url{http://www.jstatsoft.org/v33/i03}
}
\keyword{utilities}

