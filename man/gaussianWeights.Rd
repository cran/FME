\name{gaussianWeights}
\alias{gaussianWeights}

\title{
  A kernel average smoother function to weigh residuals according to a Gaussian
  density function
  
  This function is still experimental...  use with care
}

\description{
  A calibration dataset in database format (cf. modCost for the database
  format) is extended in order to fit model output using a weighted least
  squares approach. To this end, the observations are replicated for a
  certain number of times, and weights are assigned to the replicates
  according to a Gaussian density function. This density has the relevant
  observation as mean value. The standard deviation, provided as a parameter,
  determines the number of inserted replicate observations (see Detail).

  This weighted regression approach may be interesting when discontinuities
  exist in the observational data. Under these circumstances small changes 
  in the timing (or more general the position along the axis of the independent
  variable) of the model output may have a disproportional impact on the
  overall goodness-of-fit (e.g. timing of nutrient depletion). Additionally,
  this approach may be used to model uncertainty in the independent variable
  (e.g. slices of sediment profiles, or the timing of a sampling).

}
\usage{
gaussianWeights (obs, x = x, y = y, xmodel, spread, weight = "none",
                 aggregation = x ,ordering)
}

\arguments{

  \item{obs }{dataset in long (database) format as is typically used by
     modCost
  }
  \item{x }{name of the independent variable (typically x, cf. modCost) in
     \code{obs}. Defaults to x (not given as character string; cf. subset)
  }
  \item{y }{name of the dependent variable in \code{obs}. Defaults to y.
  }
  \item{xmodel }{an ordered vector of unique times at which model output
     is produced. If not given, the independent variable of the observational
     dataset is used.
  }
  \item{spread }{standard deviation used to calculate the weights from a
     normal density function. This value also determines the number of
     points from the model output that are compared to a specific observa-
     tion in \code{obs} (2 * 3 * spread + 1; containing 99.7\% of the
     Gaussian distribution, centered around the observation of interest).
  }
  \item{weight }{scaling factor of the modCost function ("sd", "mean", or
     "none"). The Gaussian weights are multiplied by this factor to account
     for differences in units.
  }
  \item{aggregation }{vector of column names from the dataset that are used
    to aggregate observations while calculating the scaling factor. Defaults
    to the variable name, "name".
  }
  \item{ordering }{Optional extra grouping and ordering of observations. Given
    as a vector of variable names. If none given, ordering will be done by
    variable name and independent variable. If both aggregation and ordering
    variables are given, ordering will be done as follows:
    x within ordering (in reverse order) within aggregation (in reverse order).
    Aggregation and ordering should be disjoint sets of variable names.
  }
}
\value{
  
  A modified version of \code{obs} is returned with the following extensions:

  1. Each observation obs[i] is replicated n times were n represents the number
     of \code{modelx} values within the interval [obs_i - (3 * spread), obs_i +
     3 * spread)].

  2. These replicate observations get the same \code{x} values as their
     modeled counterparts (\code{xmodel}).

  3. Weights are given in column, called "err"

  The returned data frame has the following columns:
  \itemize{
    \item{"name" or another name specified by the first element of 
    \code{aggregation}. Usually this column contains the names of the
    observed variables.
    }
    \item{"x" or another name specified by \code{x}
    }
    \item{"y" or another name specified by \code{y}
    }
    \item{"err" containing the calculated weights
    }
    \item{The rest of the columns of the data frame given by \code{obs} in
    that order.
    }
  }
}

\details{

  Suppose: spread = 1/24 (days; = 1 hour)
           x      = time in days, 1 per hour

  Then:
   obs_i is replicated 7 times (spread = observational periodicity = 1 hour):

   =>  obs_i-3 = ... = obs_i-1 = obs_i = obs_i+1 = ... = obs_i+3

  The weights (W_i+j, for j = -3 ...3) are calculated as follows:
       W'_i+j = 1/(spread * sqrt(2pi)) * exp(-1/2 * ((obs_i+j - obs_i)/spread)^2

       W_i+j  = W'_i+j/sum(W_i-3,...,W_i+3)
         (such that their sum equals 1)

}

\examples{
## =======================================================================
## A Sediment example
## =======================================================================

## Sediment oxygen concentration is measured every
## centimeter in 3 sediment types
depth <- 0:7
observations <- data.frame(
                    profile = rep(c("mud","silt","sand"), each=8),
                    depth   = depth,
                    O2      = c(c(6,1,0.5,0.1,0.05,0,0,0),
                                c(6,5,3,2,1.5,1,0.5,0),
                                c(6,6,5,4,3,2,1,0)
                              )
                )

## A model generates profiles with a depth resolution of 1 millimeter
modeldepths <- seq(0, 9, by = 0.05)

## All these model outputs are compared with  weighed observations.
gaussianWeights(obs = observations, x = depth, y = O2,
                xmodel = modeldepths,
                spread = 0.1, weight = "none", 
                aggregation = profile)



# Weights of one observation in silt at depth 2:
Sub <- subset(observations, subset = (profile == "silt" & depth == 2))
plot(Sub[,-1])
SubWW <- gaussianWeights(obs = Sub, x = depth, y = O2, 
                xmodel = modeldepths, spread = 0.5, 
                weight="none", aggregation = profile)
SubWW[,-1]
}

\author{
  Tom Van Engeland <tom.vanengeland@nioz.nl>
}

\keyword{utilities}

