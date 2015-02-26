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
histmdl(x, model = "Witteveen", gain = 0, support = 4,
        plot = TRUE, add = FALSE,
        density = NULL, angle = 45, col = NULL,
        border = par("fg"), lty = NULL,
        main = paste("Histogram of", xname), sub = NULL,
        xlab = xname, ylab = "Density",
        xlim = range(r$breaks), ylim = range(0, r$density),
        axes = TRUE, ann = TRUE, \ldots)
}
\arguments{
  \item{x}{a vector of values for which the histogram is desired.}
  \item{model}{a character string naming the desired histogram-like
    model.}
  \item{gain}{minimum required complexity reduction before an additional
    interval is accepted.}
  \item{support}{minimum number of data points desired per interval.}
  \item{plot}{logical.  If \code{TRUE} (default), a histogram is
    plotted, otherwise a list of breaks and counts is returned.}
  \item{add}{logical.  If \code{FALSE} (default), a new plot is created,
    otherwise the histogram is added to the active plot.}
  \item{density}{the density of shading lines, in lines per inch.
    The default value of \code{NULL} means that no shading lines
    are drawn. Non-positive values of \code{density} also inhibit the
    drawing of shading lines.}
  \item{angle}{the slope of shading lines, given as an angle in
    degrees (counter-clockwise).}
  \item{col}{a colour to be used to fill the bars.
    The default of \code{NULL} yields unfilled bars.}
  \item{border}{the color of the border around the bars.  The default
    is to use the standard foreground color.}
  \item{lty}{the line type used for the bars, see also
    \code{\link{lines}}.}
  \item{main, sub, xlab, ylab}{these arguments to \code{title} have
    useful defaults here.}
  \item{xlim, ylim}{the range of x and y values with sensible defaults.
    Note that \code{xlim} is \emph{not} used to define the histogram
    (breaks), but only for plotting (when \code{plot = TRUE}).}
  \item{axes}{logical.  If \code{TRUE} (default), axes are draw if the
    plot is drawn.}
  \item{ann}{logical.  Should annotations (titles and axis titles) be
    plotted?}
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
