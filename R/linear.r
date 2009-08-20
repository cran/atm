linear = function (x, y=NULL, key="x", E=atm.ls, s, w=NULL)
{	if (missing (s) ) s = deparse (substitute (x) )
	d = linear.design (key)
	e = linear.estimate (key)
	t = term (d, E, e, s, x, 0, w=w)
	z = NULL
	for (k in d$ks) z = cbind (z, linear.expand (t$nr, k (x) ) )
	if (!t$clean) z = cbind (z [t$valid,])
	if (!all (is.finite (z) ) )
		stop ("key functions must produce valid z given valid x")
	t$z = z
	term.pump (t, y, "linear")
}

categorical = function (x, y=NULL, E=atm.ls, s, w=NULL)
{	if (missing (s) ) s = deparse (substitute (x) )
	if (!inherits (x, "factor") ) x = factor (x)
	kstr = NULL
	levs = levels (x)
	for (lev in levs)
		kstr = c (kstr, paste ("as.integer (x=='", lev, "')", sep="") )
	t = linear (x, NULL, kstr, E, s, w)
	t$e$labs = levs
	term.pump (t, y, "categorical")
}

polynomial = function (x, y=NULL, deg=1, intercept=(!is.null (y) ), E=atm.ls, s, w=NULL)
{	if (missing (s) ) s = deparse (substitute (x) )
	kstr = character (deg)
	for (i in 1:deg) kstr [i] = paste ("x^", i, sep = "")
	if (intercept) kstr = c ("1", kstr)
	term.pump (linear (x, NULL, kstr, E, s, w), y, "polynomial")
}

linear.design = function (key)
{	d = extend (term.design (), "linear.design")
	d$np = length (key)
	d$ks = list ()
	for (i in 1:d$np) d$ks [[i]] = mutate (function (x) NULL, key [i])
	d
}

linear.estimate = function (labs, th=rep (0, length (labs) ) )
{	e = extend (term.estimate (), "linear.estimate")
	e$labs = labs
	e$th = th
	e
}

reset.linear = function (t, ...) t$e$th [] = 0

fit.linear = function (t, y, pump, ...)
{	th = t$E (t$z, clean.response (t, y), clean.weights (t) )
	if (pump) t$e$th = th
	else linear.estimate (t$e$labs, th)
}

atm.ls = function (x, y, w)
	as.vector ( (if (is.null (w) ) lm.fit (x, y) else lm.wfit (x, y, w) )$coefficients)

evaluate.linear = function (t, u, ...)
{	v = 0
	for (i in 1:t$d$np) v = v + t$e$th [i] * t$d$ks [[i]] (u)
	v
}

#not correct, ignores w...
#fails for one parameter term
linear.ftest = function (t, y, ...)
{	r1 = clean.residuals (t, y)
	r2 = clean.response (t, y) - mean (y)
	ssres1 = sum (r1 ^ 2)
	ssres2 = sum (r2 ^ 2)
	dfa = t$d$np - 1
	dfb = length (r1) - t$d$np
	fstat = ( (ssres2 - ssres1) / dfa) / (ssres1 / dfb)
	pvalue = 1 - pf (fstat, dfa, dfb)
	list ("fstat" = fstat, "pvalue" = pvalue)
}

#not correct, ignores w...
oospsummary.linear = function (t, y, ...)
{	s = structure (list (), class="summary.linear")
	s$label = term.label (t)
	s$weighted = (!is.null (t$w) )
	df = t$nv - t$d$np
	resvar = sum (clean.residuals (t, y) ^ 2) / df
	e = list ()
	e$labs = t$e$labs
	e$th = t$e$th
	e$se = sqrt (diag (resvar * solve (base::t (t$z) %*% t$z) ) )
	e$tv = e$th / e$se
	e$pv = 2 * (1 - pt (abs (e$tv), df) )
	e$stars = pstars (e$pv)
	s$e = data.frame (e)
	s$ft = linear.ftest (t, y)
	s$gf = data.frame (gf (t, y) )
	s
}

oospprint.summary.linear = function (s, ...)
{	cat (s$label, "\n")
	if (s$weighted) cat ("WARNING: WEIGHTING IGNORED IN SUMMARY\n")
	cat (s$e, "\n", s$ft, "\n", s$gf, "\n")
}

oospfitted.linear = function (t, ...)
{	if (t$clean) t$z %*% t$e$th
	else evaluate (t, t$x)
}
oospprint.linear.design = function (d, ...) print (d$ks)
oospprint.linear.estimate = function (e, ...) print (data.frame (labs=e$labs, th=e$th) )

linear.expand = function (n, x)
{	if (length (x) == 1 && n > 1) rep (x, n)
	else (x)
}


