# Binary search, equivalent to max (which (x <= value)), but faster
startIndex <- function (x, value) .Call (maxLE, x, value)


# If possible, return an interval covering more than `support` points
bestInterval <- function (x, support=3) {
	m <- length (x)
	M <- x[m] - x[1]
	x <- rle (x)
	# The forest of Kruskals algorithm as a vector of interval left sides
	K <- x$values
	best <- sapply (1:length (K), function (i)
		list (left=K[i], right=K[i], length=x$lengths[i], complexity=Inf))
	for (left in x$values[order (diff (x$values))]) {
		i <- startIndex (K, left)
		# Construct the candidate interval details
		ch <- list (left=best[, i]$left, right=best[, i + 1]$right,
		            length=as.integer (NA), complexity=as.numeric (NA))
		MIn <- ch$right - ch$left
		mIn <- sum (x$lengths[startIndex (x$values, ch$left) :
		                      startIndex (x$values, ch$right)])
		ch$length <- mIn
		ch$complexity <- mIn * log2 (MIn / mIn) + if (mIn != m)
			(m - mIn) * log2 ((M - MIn) / (m - mIn)) else 0
		# Unite the intervals: this is a speed bottleneck
		K <- K[-(i + 1)]
		if (best[, i]$complexity < best[, i + 1]$complexity ||
		    best[, i + 1]$length < support)
			best <- best[, -(i + 1), drop=FALSE]
		else
			best <- best[, -i, drop=FALSE]
		# Consider the candidate interval
		if (ch$complexity < best[, i]$complexity || best[, i]$length < support)
			best[, i] <- ch
	}
	drop (best)
}


# The returned densities are not probability densities
# x: a sorted non-empty double-precision vector
# gain: minimum complexity gain (units depend on the base of the logarithms)
# support: the minimum number of data points per bin
recursiveIntervals <- function (x, gain=0, support=3) {
	m <- length (x)
	breaks <- c (x[1], x[m])
	density <- m / (breaks[2] - breaks[1])
	if (m >= 2 * support &&
	    (ch <- bestInterval (x, support))$complexity +
	    m * log2 (density) + 2 * log2 (m) - 1 < -gain) {
		width <- ch$right - ch$left
		iL <- startIndex (x, ch$left)
		while (identical (x[iL], x[iL - 1])) iL <- iL - 1
		iR <- startIndex (x, ch$right)
		rIn <- recursiveIntervals (x[iL : iR], gain, support)
		# It is crucial that the equality is in the non-shifted part because
		# ch$right - width may not evaluate equal to ch$left (rounding)
		rOut <- recursiveIntervals (c (head (x, iL),
		                               tail (x, -iR) - width), gain, support)
		breaks <- c (rOut$breaks[rOut$breaks < ch$left],
		             rIn$breaks,
		             rOut$breaks[rOut$breaks > ch$left] + width)
		density <- c (rOut$density[head (rOut$breaks, -1) < ch$left],
		              rIn$density,
		              rOut$density[tail (rOut$breaks, -1) > ch$left])
		if (length (breaks) != length (density) + 1)
			stop ("Malformed output. Internal error.")
	}
	list (breaks=breaks, density=density)
}


histmdl <- function (x, model="Witteveen", gain=0, support=4, plot=TRUE,
                     add=FALSE, density=NULL, angle=45, col=NULL,
                     border=par ("fg"), lty=NULL,
                     main=paste ("Histogram of" , xname), sub=NULL,
                     xlab=xname, ylab="Density", xlim=range (r$breaks),
                     ylim=range (0, r$density), axes=TRUE, ann=TRUE, ...) {
	xname <- paste (deparse (substitute (x), 500), collapse="\n")
	x <- sort (as.numeric (x))
	if (length(x) == 0 || anyNA (x))
		stop ("'x' must be non-empty and solely numeric")
	if (model != "Witteveen")
		stop (paste ("Unsupported model:", model))
	r <- structure (c (recursiveIntervals (x, gain, support),
	                   list (equidist=FALSE, xname=xname)),
	                class="histogram")
	# By contracting intervals, some points may get counted multiple times
	# r$density <- r$density / length (x)
	r$density <- r$density / sum (r$density * diff (r$breaks))
	if (plot) {
		dev.hold (); on.exit (dev.flush ())
		if (!add) {
			plot.new ()
			plot.window (xlim, ylim, "")
			if (ann)
				title (main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
			if (axes) {
				axis (1, ...)
				axis (2, ...)
			}
		}
		rect (head (r$breaks, -1), 0, tail (r$breaks, -1), r$density,
		      col=col, border=border, angle=angle, density=density, lty=lty)
		invisible (r)
	} else {
		r
	}
}

