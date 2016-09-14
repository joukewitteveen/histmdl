# Binary search, equivalent to max (which (x <= value)), but faster
startIndex <- function (x, value) .Call (maxLE, x, value)


bestInterval <- function (x, precision=0, support=3) {
	m <- length (x)
	M <- x[m] - x[1]
	x <- rle (x)
	# The forest of Kruskals algorithm as a vector of interval left sides
	K <- x$values
	best <- mapply (function (value, count) {
		if (precision > 0 && count >= support)
			complexity <- count * log2 (precision / count) +
		                  (m - count) * log2 (M / (m - count))
		else
			# Prevent infinite densities
			complexity <- Inf
		list (left=value, right=value, count=count, complexity=complexity)
	}, x$values, x$lengths)
	for (left in x$values[order (diff (x$values))]) {
		i <- startIndex (K, left)
		# Construct the candidate interval details
		ch <- list (left=best[, i]$left, right=best[, i + 1]$right,
		            count=as.integer (NA), complexity=as.numeric (NA))
		MIn <- ch$right - ch$left
		mIn <- sum (x$lengths[startIndex (x$values, ch$left) :
		                      startIndex (x$values, ch$right)])
		ch$count <- mIn
		ch$complexity <- mIn * log2 ((MIn + precision) / mIn) + if (mIn != m)
			(m - mIn) * log2 ((M - MIn) / (m - mIn)) else 0
		# Unite the intervals: this is a speed bottleneck
		K <- K[-(i + 1)]
		if (best[, i]$complexity < best[, i + 1]$complexity ||
		    best[, i + 1]$count < support)
			best <- best[, -(i + 1), drop=FALSE]
		else
			best <- best[, -i, drop=FALSE]
		# Consider the candidate interval
		if (ch$complexity < best[, i]$complexity || best[, i]$count < support)
			best[, i] <- ch
	}
	drop (best)
}


# The returned densities are not probability densities
# x: a sorted non-empty double-precision vector
# gain: minimum complexity gain (units depend on the base of the logarithms)
# precision: the minimum resolution at which `x` is measured
# support: the minimum number of data points per bin
recursiveIntervals <- function (x, gain=0, precision=0, support=3) {
	m <- length (x)
	breaks <- c (x[1], x[m])
	if (breaks[1] < breaks[2])
		density <- m / (breaks[2] - breaks[1])
	else
		density <- m / precision
	if (m >= 2 * support && breaks[1] < breaks[2] &&
	    (ch <- bestInterval (x, precision, support))$complexity +
	    m * log2 (density) + 2 * log2 (m) - 1 < -gain) {
		width <- ch$right - ch$left
		iL <- startIndex (x, ch$left)
		while (identical (x[iL], x[iL - 1])) iL <- iL - 1
		iR <- startIndex (x, ch$right)
		rIn <- recursiveIntervals (x[iL : iR], gain, precision, support)
		# It is crucial that the left boundary point is kept because
		# ch$right - width may not evaluate equal to ch$left (rounding)
		rOut <- recursiveIntervals (c (head (x, iL), tail (x, -iR) - width),
		                            gain, precision, support)
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


histmdl <- function (x, model="Witteveen", gain=0, precision=0, support=4,
                     plot=TRUE, main=paste ("Histogram of", xname),
                     xlab=xname, ylab="Density", ...) {
	xname <- paste (deparse (substitute (x), 500), collapse="\n")
	if (!is.numeric (x))
		stop ("'x' must be numeric")
	# Sorting removes NAs
	x <- sort (as.double (x))
	if (length (x) == 0)
		stop ("'x' must be non-empty")
	if (model != "Witteveen")
		stop (paste ("Unsupported model:", model))
	if (!is.numeric (precision) || precision < 0)
		stop ("'precision' must be non-negative")
	if (!is.numeric (support) || support < 2)
		stop ("'support' must be at least 2")
	if (precision == 0)
		for (w in with (rle (x),
		                paste ("infinite density ignored at value", values,
		                       "occurring", lengths, "times",
		                       "while 'precision' is 0")[lengths >= support]))
			warning (w)
	r <- structure (c (recursiveIntervals (x, gain, precision, support),
	                   list (equidist=FALSE, xname=xname)),
	                class="histogram")
	# By contracting intervals, some points may get counted multiple times
	# r$density <- r$density / length (x)
	r$density <- r$density / sum (r$density * diff (r$breaks))
	if (plot) {
		plot (r, freq=FALSE, labels=FALSE,
		      main=main, xlab=xlab, ylab=ylab, ...)
		invisible (r)
	} else {
		r
	}
}

