atm = function (y, ..., fit=TRUE, rs=TRUE, maxfits=25, w=NULL)
{	ts = list (...)
	if (inherits (ts [[1]], "list") ) ts = ts [[1]]
	m = extend (termenv (), "atm")
	m$th = 0
	m$ts = ts
	m$nt = length (ts)
	m$y = y
	m$nr = length (y)
	m$converge = FALSE
	m$nfits = 0
	if (!all (is.finite (y) ) ) stop ("all y must be valid")
	if (!is.null (w) ) for (t in ts) term.reweight (t, w)
	if (fit) fit (m, rs, maxfits)
	else m
}

fit.atm = function (m, rs=TRUE, maxfits=25, intermediate.diagnostics=FALSE, ...)
{	ts = m$ts
	if (rs) for (t in ts) reset (t)
	mar1 = mar2 = marb = NA
	m$th= mean (oospresiduals (m) )
	while (!m$converge && m$nfits <= maxfits)
	{	es = list ()
		for (i in 1:m$nt)
			es [[i]] = fit (m$ts [[i]], oospresiduals (m, partial=i), FALSE)
		for (i in 1:m$nt)
		{	#does this destroy the environment?
			m$ts [[i]]$e = NULL
			m$ts [[i]]$e = es [[i]]
		}
		m$th = m$th + mean (oospresiduals (m) )
		r = oospresiduals (m)
		mar0 = mar1
		mar1 = mean (abs (r) )
		if (is.na (mar0) ) marb = mar1
		#back the front (should be 1 - ...), however doesn't matter
		changeb = (marb - mar1) / marb
		changei = (mar0 - mar1) / mar0
		#convergence criterion satisfied if
		#ith fit 1000 times better than first fit,
		#ith fit has no great improvement over previous fit
		if (!is.na (changei) )
			if (changeb > 0.999 || changei < 0.001) m$converge = TRUE
		if (intermediate.diagnostics)
		{	oospplot (m)
			cat ("changeb: ", changeb, ", changei: ", changei, "\n", sep="") 
			locator (1)
		}
		m$nfits = m$nfits + 1
	}
	if (!m$converge) warning ("atm didn't converge")
	m
}

atmx = function (m)
{	x = list ()
	for (i in 1:m$nt) x [[i]] = m$ts [[i]]$x
	x
}

evaluate.atm = function (m, us, ...)
{	eh = m$th
	for (i in 1:m$nt) eh = eh + evaluate (m$ts [[i]], us [[i]])
	eh
}

gf.atm = function (m, ...)
{	r = oospresiduals (m)
	yb = mean (m$y)
	yh = oospfitted (m)
	simple.gf (yb, yh, r, ...)
}

oospsummary.atm = function (m, which=NA, ...)
{	if (is.na (which) )
	{	s = structure (list (), class="summary.atm")
		r = oospresiduals (m)
		s$th = m$th
		s$ts = list ()
		for (t in m$ts) s$ts = c (s$ts, oospsummary (t, r) )
		s$gf = data.frame (gf (m) )
		s$converge = m$converge
		s
	}
	else oospsummary (m$ts [[which]], oospresiduals (m, which) )
}

oospprint.summary.atm = function (s, ...) print.default (s)

oospfitted.atm = function (m, ...)
{	#too slow????
	#atm.eval (m, atmx (m) )
	yh = m$th
	for (t in m$ts) yh = yh + oospfitted (t)
	yh
}

#atms should always force residuals
#r != y - yh, r = (1 - n) * y - th + r1 + r2 + ..., if partial one term is excluded
oospresiduals.atm = function (m, y=m$y, partial=NULL, ...)
{	r = (1 - m$nt) * y - m$th
	for (t in m$ts) r = r + force.residuals (t, y)
	if (!is.null (partial) ) r = fitted (m$ts [[partial]]) + r
	r
}


