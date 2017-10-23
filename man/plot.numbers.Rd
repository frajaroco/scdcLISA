\name{plot.numbers}
\alias{plot.numbers}
\title{Plot of label points}
\description{Draws a text label at a given x and y coordinate into a spatial point pattern}
\usage{\method{plot}{numbers}(x,...)}

\arguments{
\item{x}{The spatial point pattern to be plotted. An object of class "ppp"}
\item{...}{Arguments to be passed to methods, such as graphical parameters (see par)}
}
\details{If necessary, more details than the description above}
\value{
\item{Invisible}{object of class "symbolmap" giving the correspondence between labeled points and plotting characters.}
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
plot.numbers(pattern, main = "Child Leukaemia")

## End(Not run)
#################
}