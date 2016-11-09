\name{plot.numbers}
\alias{plot.numbers}
\title{Plot of label points}
\usage{plot.numbers(pattern, ...)}
\arguments{
  \item{pattern}{The spatial point pattern to be plotted. An object of class "ppp".}
  \item{\dots}{extra arguments that will be passed to the plotting functions}
}
\details{This is the plot method for point pattern datasets (of class "ppp").}
\value{(Invisible) object of class "symbolmap" giving the correspondence between labeled points and plotting characters.}
\author{Jonatan A. Gonzalez}
\examples{
## Not run:
#################

# The dataset humberside is an object of the same format, representing the same point pattern data
pattern <- humberside
plot.numbers(pattern, main = "Child Leukaemia")

## End(Not run)
}