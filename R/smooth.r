smooth = function (x, y=NULL, ns=50, refit=2, deg=2, E=atm.hs,
	wf=stdwf, smoothness=NULL, s, w=NULL)
{	if (missing (s) ) s = deparse (substitute (x) )
	if (deg == 0) stop ("degree must be at least one")
	d = smooth.design (ns, refit, deg, wf, smoothness, x)
	e = smooth.estimate (ns)
	t = term (d, E, e, s, x, w=w)
	term.pump (t, y, "smooth")
}

smooth.design = function (ns, refit, deg, wf, smoothness, x)
{	d = extend (term.design (), "smooth.design")
	xr = range (x, na.rm=TRUE)
	d$ns = ns
	d$sx = seq (xr [1], xr [2], length=ns)
	d$refit = refit
	d$deg = deg
	d$wf = wf
	d$smoothness = if (is.null (smoothness) ) 0.85 * diff (xr) else smoothness
	d
}

smooth.estimate = function (ns, sy=rep (0, ns) )
{	e = extend (term.estimate (), "smooth.estimate")
	e$sy = sy
	e
}

reset.smooth = function (t, ...) t$e$sy [] = 0

fit.smooth = function (t, y, pump, ...)
{	sy = t$E (t$z, clean.response (t, y), t$d$ns, t$d$sx, t$d$deg, 
		t$d$wf, t$d$smoothness, clean.weights (t) )
	if (t$d$refit > 1) for (i in 2:t$d$refit)
		sy = t$E (t$d$sx, sy, t$d$ns, t$d$sx, t$d$deg, t$d$wf, t$d$smoothness, NULL)
	if (pump) t$e$sy = sy
	else smooth.estimate (t$d$ns, sy)
}

evaluate.smooth = function (t, u, ...) spline (t$d$sx, t$e$sy, xout=u)$y

oospsummary.smooth = function (t, y, ...)
{	s = structure (list (), class="summary.smooth")
	s$label = term.label (t)
	s$ns = t$d$ns
	s$deg = t$deg
	s$sx = t$d$sx
	s$sy = t$e$sy
	s$gf = data.frame (gf (t, y) )
	s
}

oospprint.summary.smooth = function (s, ...)
	cat (s$label, "\nns: ", s$ns, "\ndeg: ",
		s$deg, "\nsx: ", s$sx, "\nsy", s$sy, "\n", s$gf, "\n")

oospprint.smooth.estimate = function (e, ...) print (t$e$sy)

atm.ss = function (x, y, ns, sx, deg, wf, smoothness, w)
{	z = NULL
	for (i in 0:deg) z = cbind (z, x ^ i)
	sy = numeric (ns)
	hs = smoothness / 2
	yb = mean (y)
	for (i in 1:ns)
	{	u = sx [i]
		v = x >= u - hs & x <= u + hs
		sy [i] = if (sum (v) > deg)
		{	xsub = x [v]
			ysub = y [v]
			zsub = z [v,]
			wst = wf (u, smoothness, xsub)
			if (!is.null (w) ) wst = w * wst
			sum (lm.wfit (zsub, ysub, wst)$coefficients * u ^ (0:deg) )
		}
		else yb
	}
	sy
}

atm.hs = function (x, y, ns, sx, deg, wf, smoothness, w)
{	a = min (x)
	b = max (x)
	if (b - a > smoothness) atm.ss (x, y, ns, sx, deg, wf, smoothness, w)
	z = NULL
	for (i in 0:deg) z = cbind (z, x ^ i)
	sy = numeric (ns)
	hs = smoothness / 2
	a = a + hs
	b = b - hs
	yb = mean (y)
	for (i in 1:ns)
	{	u = sx [i]
		v = x >= u - hs & x <= u + hs
		sy [i] = if (sum (v) > deg)
		{	xsub = x [v]
			ysub = y [v]
			zsub = z [v,]
			wst = hardwf (a, b, wf, u, smoothness, xsub)
			if (!is.null (w) ) wst = w * wst
			sum (lm.wfit (zsub, ysub, wst)$coefficients * u ^ (0:deg) )
		}
		else yb
	}
	sy
}

hardwf = function (a, b, wf, m, w, x)
{	if (m < a)
	{	i1 = x < m
		i3 = x > a
		x1 = x [i1]
		n2 = sum (!(i1 | i3) )
		x3 = x [i3]
		y1 = wf (m, w, x1)
		y2 = rep (1, n2)
		y3 = wf (a, w, x3)
		c (y1, y2, y3)
	}
	else if (m > b)
	{	i1 = x < b
		i3 = x > m
		x1 = x [i1]
		n2 = sum (!(i1 | i3) )
		x3 = x [i3]
		y1 = wf (b, w, x1)
		y2 = rep (1, n2)
		y3 = wf (m, w, x3)
		c (y1, y2, y3)
	}
	else wf (m, w, x)
}

stdwf = function (m, w, x)
{	if (m != 0) x = x - m
	z = 2 * x / w
	1 - z * z
}

