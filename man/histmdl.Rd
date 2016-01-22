\name{histmdl}
\alias{histmdl}
\title{
  Most Informative Histograms
}
\description{
  Compute an estimate of the maximum likelihood parameter of a
  histogram-like model.  If \code{plot=TRUE}, the resulting object of
  \link{class} \code{"histogram"} is plotted according to
  \code{\link{plot.histogram}}, before it is returned.
}
\usage{
histmdl(x, model = "Witteveen", gain = 0, precision = 0, support = 4,
        plot = TRUE, main = paste("Histogram of", xname),
        xlab = xname, ylab = "Density", \ldots)
}
\arguments{
  \item{x}{a vector of values for which the histogram is desired.}
  \item{model}{a character string naming the desired histogram-like
    model.  Currently, only \code{"Witteveen"} is implemented.}
  \item{gain}{minimum required complexity reduction before an additional
    interval is accepted.}
  \item{precision}{a value giving the minimum resolution of the data.
    When computing complexities, boundary values are blurred by this
    amount to mitigate the effects of, for example, rounding.}
  \item{support}{minimum number of data points desired per interval.}
  \item{plot}{logical.  If \code{TRUE} (default), a histogram is
    plotted, otherwise a list of breaks and densities is returned.}
  \item{main, xlab, ylab}{these arguments to \code{title} have useful
    defaults here.}
  \item{\dots}{further arguments and \link{graphical parameters} passed
    to \code{\link{plot.histogram}} and thence to \code{\link{title}}
    and \code{\link{axis}} (if \code{plot=TRUE}).}
}
%\details{
%}
\value{
  an object of class \code{"histogram"} which is a list of components:
  \item{breaks}{the boundaries of intervals.  Note that consecutive
    values are not the same as the boundaries of the intervals that,
    recursively, define the model instance.}
  \item{density}{densities of the data inside the model-interval that a
    section is part of.}
  \item{xname}{a character string with the actual \code{x} argument name.}
}
%\references{ TBA }
\author{
  Jouke Witteveen and Richard Gill
}
\seealso{
  \code{\link{hist}}
}
\examples{
  set.seed (28011988)

  x <- c (rnorm (1000, -6), rnorm (1000, 6))
  histmdl (x, gain=2, col="peru")
  hist (x, freq=FALSE, add=TRUE, col="#80808080")

  x <- c (runif (50), runif (50, max=3))
  histmdl (x, col="peru", ylim=0:1)
  hist (x, freq=FALSE, add=TRUE, col="#80808080")
}
\keyword{histogram}
\keyword{dplot}
\keyword{hplot}
\keyword{distribution}
