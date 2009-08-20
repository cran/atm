#abstract superclass
term = function (d, E, e, s, x, z=NULL, validate=TRUE,
	clean=TRUE, nr=length (x), nc=1, nv=nr, valid=NULL, w=NULL)
{	t = extend (termenv (), "term")
	t$d = d	#term design, depends on subclass
	t$E = E	#term estimator (function)
	t$e = e	#term estimate, depends on subclass
	t$s = s	#symbol, a character representing the name of x
	t$x = x	#realisations of a real world variable (anything theoretically measurable)
	t$z = z	#realisations of a model world variable (e.g. a vector or matrix)
	t$clean = clean	#all realisations valid (nr equals nv)
	t$nr = nr	#number of realisations of x (for matrix, number of rows)
	t$nc = nc	#number of component variables of x (for matrix, number of columns)
	t$nv = nv	#number of valid realisations of x
	t$valid = valid	#valid index
	t$w = w	#weights, unless null, nr should equal length (w)
	if (validate)
	{	term.validate (t)
		if (is.null (t$z) ) t$z = clean.explanatory (t)
	}
	t
}

term.label = function (t) paste (class (t) [1], t$s, sep = ":")

term.validate = function (t, vz)
{	x = t$x
	valid = NULL
	if (inherits (x, "matrix") || inherits (x, "data.frame") )
	{	t$nr = nrow (x)
		t$nc = ncol (x)
		valid = logical (t$nr)
		for (i in 1:t$nr) valid [i] = all (is.finite (x [i,]) )
	}
	else if (is.vector (x) || is.factor (x) ) valid = is.finite (x)
	else stop ("term.validate can not be applied")
	clean = all (valid)
	if (! clean)
	{	t$clean = FALSE
		t$nv = sum (valid)
		t$valid = valid
	}
	if (t$nv < 2) stop ("terms must have at least 2 valid realisations")
}

#maybe rename these
#"clean" is with respect to x
clean.explanatory = function (t)
	if (t$clean) t$x else if (t$nc > 1) t$x [t$valid,] else t$x [t$valid]
clean.response = function (t, y) if (t$clean) y else y [t$valid]
clean.weights = function (t) if (t$clean || is.null (t$w) ) t$w else t$w [t$valid]

term.pump = function (t, y, sc=NULL, ...)
{	t = extend (t, sc)
	if (!is.null (y) ) fit (t, y, TRUE, ...)
	t
}

termenv = function () extend (compenv (), "termenv")
term.design = function () extend (termenv (), "term.design")
term.estimate = function () extend (termenv (), "term.estimate")

gf.term = function (t, y, method=c ("rsq", "mar"), ...)
{	yb = yh = NULL
	if (any (method == "rsq") )
	{	yb = mean (clean.response (t, y) )
		yh = clean.fitted (t)
	}
	r = clean.residuals (t, y)
	simple.gf (yb, yh, r, method)
}

simple.gf = function (yb, yh, r, method=c ("rsq", "mar") )
{	gf = list ()
	if (any (method == "mar") ) gf$mar = mean (abs (r) )
	if (any (method == "msr") ) gf$msr = mean (r ^ 2)
	if (any (method == "ssr") ) gf$ssr = sum (r ^ 2)
	if (any (method == "rsq") ) gf$rsq = sum ( (yh - yb) ^ 2) / sum (r ^ 2)
	gf
}

#these don't look right????
term.confband = function (t, y, u=seq (t), type="standard", cl=0.95, ...)
{	r = clean.residuals (t, y)
	n = length (r)
	v = evaluate (t, u)
	e = if (type == "standard") qt (0.5 + cl / 2, n - 1) * sd (r) / n
	else if (type == "prediction") qt (0.5 + cl / 2, length (r) - 1) * sd (r)
	else stop ("confband type can only be standard or prediction")
	cbind (v - e, v + e)
}

term.predband = function (t, y, ...) term.confband (t, y, type="prediction", ...)
term.reweight = function (t, w) t$w = w
oospprint.term = function (t, ...) cat (term.label (t), "\n")

oospsummary.term = function (t, ...)
{	s = structure (list (), class="summary.term")
	s$label = term.label (t)
	s
}

oospprint.summary.term = function (s, ...) cat (s$label, "\n")

oosppredict.term = function (t, y, u=t$x, type=c("argument", "value"), cl=0.95, ...)
{	p = list ()
	if (any (type == "argument") ) p$u = u
	if (any (type == "value") ) p$v = evaluate (t, u)
	if (any (type == "standard") )
	{	cb = term.confband (t, y, u, cl=cl)
		p$confl = cb [,1]
		p$confu = cb [,2]
	}
	if (any (type == "prediction") )
	{	cb = term.predband (t, y, u, cl=cl)
		p$predl = cb [,1]
		p$predu = cb [,2]
	}
	data.frame (p)
}

#valid x ****MUST**** produce valid fitted values
oospfitted.term = function (t, ...) evaluate (t, t$x)
oospresiduals.term = function (t, y, ...) y - oospfitted (t)
oospweights.term = function (t, ...) t$w

oospprint.term.estimate = function (e, ...) print.default (e)

seq.term = function (t, n=200, ...)
{	if (t$nc > 1) stop ("seq.term reqiures nc = 1")
	if (inherits (t$x, "factor") ) 1:length (levels (t$x) )
	else
	{	x = range (t$x, na.rm=TRUE)
		seq (x [1], x [2], length=n)
	}
}

clean.fitted = function (t)
{	yh = oospfitted (t)
	if (!t$clean) yh = yh [t$valid]
	yh
}

#assuming that all y valid
#assuming that invalid x is the ****ONLY**** possible cause of invalid r
clean.residuals = function (t, y)
{	r = oospresiduals (t, y)
	if (!t$clean) r = r [t$valid]
	r
}

#same comments
force.residuals = function (t, y)
{	r = oospresiduals (t, y)
	if (!t$clean) r [!t$valid] = mean (r [t$valid])
	r
}

pstars = function (ps)
{	v = character ()
	for (p in ps)
	{	if (p > 0.05 && p <= 0.1) v = c (v, "   .")
		else
		{	n = 0
			if (p <= 0.000001) n = 4
			else if (p <= 0.001) n = 3
			else if (p <= 0.01) n = 2
			else if (p <= 0.05) n = 1
			str1 = paste (rep (" ", 4 - n), collapse = "")
			str2 = paste (rep ("*", n), collapse = "")
			v = c (v, paste (str1, str2, sep = "") )
		}
	}
	v
}

