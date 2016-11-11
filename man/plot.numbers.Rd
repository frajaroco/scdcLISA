\name{plot.numbers}
\alias{plot.numbers}
\title{Plot of label points}
\usage{
plot.numbers(xy, ...)
}
\arguments{
\item{xy}{The spatial point pattern to be plotted. An object of class "ppp".}

\item{...}{Arguments to be passed to methods, such as graphical parameters (see par). }
\value{
(Invisible) object of class "symbolmap" giving the correspondence between labeled points and plotting characters.
}
\description{
Draws a text label at a given x and y coordinate into a spatial point pattern.
}
\examples{
## Not run:
#################

# Realisations of the homogeneous Poisson processes
pattern <- humberside
plot.numbers(pattern, main = "Child Leukaemia")

## End(Not run)
}
\author{
Jonatan A Gonzalez
}

