\name{plot.LISA}
\alias{plot.LISA}
\title{Plot of LISA functions}
\description{Plot method for LISA functions of a point pattern}
\usage{plot.LISA(pl, ...)}
\arguments{
\item{pl}{An object of class "list", containing four variables attached to LISA functions (usually calculated through \code{gLISA} function).}
\item{\dots}{Arguments to be passed to methods, such as graphical parameters (see par)}
}
\details{If necessary, more details than the description above}
\value{
\item{Invisible}{object of class "NULL"  giving the meaning of the different line types and colours.}
}
\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Baddeley, A., Rubak, E. and Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. In: Chapman & Hall
Interdisciplinary Statistics Series, CRC Press, Boca Raton, Florida.
}

\author{Jonatan A. Gonzaleza <jmonsalv@uji.es> and Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}}
\examples{
## Not run:
#################

# Realisations of the homogeneous Poisson processes
pattern <- humberside
pl <- gLISA(pattern, ks = "epanech", correction = "isotropic")
plot.LISA(pl)

## End(Not run)
#################
}
