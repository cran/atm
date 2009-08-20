oospplot.atm = function (m, which=NA, ...)
{	r = oospresiduals (m)
	if (is.na (which) )
	{	dims = matrix (c (1,1,1,  2,1,2,  3,1,3,  4,2,2,  6,2,3,  8,2,4, 9,3,3,  12,3,4,
		15,3,5,  18,3,6,  20,4,5,  24,4,6,  30,5,6), ncol=3, byrow=TRUE)
		dims = rbind (dims [dims [,1] >= m$nt,])
		if (nrow (dims) > 0) dims = dims [1, 2:3]
		else stop ("currently oospplot.atm limited to 30 terms")
		parst = par (mfrow=dims, oma=c (0, 0, 0, 0), mar=c (2.5, 2.5, 3.0, 1.0) )
		for (t in m$ts) oospplot (t, oospfitted (t) + r)
		par (parst)
	}
	else
	{	t = m$ts [[which]]
		oospplot (t, oospfitted (t) + r)
	}
}

oospplot.term = function (t, y=NULL, s, newplot=TRUE,
	main=term.label (t), lwd=2, col=rgb (0, 0.6, 0.1), ...)
{	if (missing (s) )
	{	s = if (is.null (y) ) paste ("fh(", t$s, ")", sep="")
		else deparse (substitute (y) )
	}
	if (t$nc > 1) stop ("oospplot.term not applicable where nc > 1")
	fx = seq (t)
	fy = evaluate (t, fx)
	if (newplot)
	{	rx = range (c (fx, t$x), na.rm=TRUE)
		ry = range (c (fy, y), na.rm=TRUE)
		plot.new ()
		plot.window (xlim=rx, ylim=ry)
		box ();
		axis (1);
		axis (2);
		title (main=main, xlab=t$s, ylab=s)
	}
	if (!is.null (y) ) points (t$x, y)
	lines (fx, fy, lwd=lwd, col=col)
}

oosplines.term = function (t, ...) oospplot (t, ..., newplot=FALSE)

oospplot.categorical = function (t, y=NULL, s,
	main=term.label (t), lwd=2, col=rgb (0, 0.6, 0.1), ...)
{	if (missing (s) )
	{	s = if (is.null (y) ) paste ("fh(", t$s, ")", sep="")
		else deparse (substitute (y) )
	}
	u = 1:t$d$np
	fy = t$e$th
	plot.new ()
	plot.window (xlim=c (0.5, t$d$np + 0.5), ylim=range (c (fy, y), na.rm=TRUE))
	box ();
	axis (1, 1:t$d$np, t$e$labs);
	axis (2);
	title (main=main, xlab=t$s, ylab=s)
	if (!is.null (y) )
	{	xst = as.integer (t$x) + runif (t$nr, -0.325, 0.325)
		points (xst, y)
	}
	lines (u, fy, lwd=lwd, col=col)
	points (u, fy, pch=16, cex=1.5, col=col)
	points (u, fy, cex=1.5)
}

oosplines.categorical = function (t, ...)
	stop ("oosplines not applicable for categorical term")

oospplot.interaction = function (t, ...)
{	plot.new ()
	plot.window (xlim=0:1, ylim=0:1)
	title (main=term.label (t) )
}

