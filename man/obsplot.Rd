\name{obsplot}
\alias{obsplot}
\title{
 Plot Method for observed data
}
\description{
Plot all observed variables in matrix formalt
}
\usage{
obsplot(x, ..., which = NULL, xyswap = FALSE, ask = NULL)
}

\arguments{
  \item{x }{a \code{matrix} or \code{data.frame}, containing the observed
    data to be plotted. The 'x'-values (first axis) should be the first
    column.      
    
    Several other matrices or data.frames can be passed in the \code{...},
    after \code{x} (unnamed) - see second example.
    
    If the first column of \code{x} consists of factors, or characters (strings),
    then it is assumed that the data are presented in long (database) format,
    where the first three columns contain (name, x, y). See last example.
    }

  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables, except the first column. }

  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.}

  \item{xyswap }{ if \code{TRUE}, then x-and y-values are swapped and the 
    y-axis is from top to bottom. Useful for drawing vertical profiles.}

  \item{\dots}{additional arguments.
  
    The graphical arguments are passed to \code{\link{plot.default}}
    and \code{points}.

    The dots may contain other matrices and data.frames with observed data
    to be plotted on the same graphs as \code{x}  - see second example.
        
    The arguments after \ldots must be matched exactly.     
  } 
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow = c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.
  Set \code{mfrow} equal to \code{NULL} to avoid the plotting function to 
  change user-defined \code{mfrow} or \code{mfcol} settings.

  Other graphical parameters can be passed as well. Parameters
  are vectorized, either according to the number of plots 
  (\code{xlab, ylab}, \code{main, sub}, \code{xlim, ylim}, \code{log},
   \code{asp, ann, axes, frame.plot},\code{panel.first,panel.last},
   \code{cex.lab,cex.axis,cex.main}) or 
  according to the number of lines within one plot (other parameters 
  e.g. \code{col}, \code{lty}, \code{lwd} etc.) so it is possible to
  assign specific axis labels to individual plots, resp. different plotting 
  style. Plotting parameter \code{ylim}, or \code{xlim} can also be a list 
  to assign different axis limits to individual plots.
}

\seealso{
 \link[deSolve]{print.deSolve}, \link[deSolve]{ode}, \link[deSolve]{deSolve}  
}
\examples{

## 'observed' data
AIRquality <- cbind(DAY = 1:153, airquality[, 1:4])
head(AIRquality)
obsplot(AIRquality, type="l", xlab="Day since May")

## second set of observed data
AIR2 <- cbind( 1:100, Solar.R = 250 * runif(100), Temp = 90-30*cos(2*pi*1:100/365) )

obsplot(AIRquality, AIR2, type = "l", xlab = "Day since May" , lwd = 1:2)

obsplot(AIRquality, AIR2, type = "l", xlab = "Day since May" , 
         lwd = 1 : 2, which =c("Solar.R", "Temp"), 
         xlim = list(c(0, 150), c(0, 100)))

obsplot(AIRquality, AIR2, type = "l", xlab = "Day since May" , 
         lwd = 1 : 2, which =c("Solar.R", "Temp"), log = c("y", ""))

obsplot(AIRquality, AIR2, which = 1:3, xyswap = c(TRUE,FALSE,TRUE))

## ' a data.frame, with 'treatments', presented in long database format
Data <- ToothGrowth[,c(2,3,1)]
head   (Data)
obsplot(Data, ylab = "len", xlab = "dose")

# same, plotted as two observed data sets
obsplot(subset(ToothGrowth, supp == "VC", select = c(dose, len)),
        subset(ToothGrowth, supp == "OJ", select = c(dose, len)))

}

\keyword{ hplot }

