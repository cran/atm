interaction = function (x, y=NULL, E=atm.ls, s, w=NULL)
{	pd = interaction.predesign (x)
	if (missing (s) ) s = pd$s
	t = linear (x, y, pd$kstr, E, s, w)
	t$e$labs = pd$labs
	t$ncat = pd$ncat
	t$ncont = pd$ncont
	term.pump (t, y, "interaction")
}

#characters?
interaction.predesign = function (x)
{	d = list ()
	if (!inherits (x, "data.frame") ) stop ("interaction x must be data.frame")
	d$s = paste (names (x), collapse = ":")
	d$ncat = d$ncont = 0
	catk = catl = contk = contl = list ()
	for (i in 1:length (x) )
	{	v = x [[i]]
		if (inherits (v, "factor") )
		{	levs = levels (v) [table (v) > 0]
			nlevs = length (levs)
			if (nlevs == 0) stop ("interaction x contains corrupt variable")
			d$ncat = d$ncat + 1
			k = NULL
			for (lev in levs)
				k = c (k, paste ("as.integer (x[,", i, "], =='", lev, "')", sep="") )
			catk [[d$ncat]] = k
			catl [[d$ncat]] = levs		
		}
		else if (mode (v) == "numeric")
		{	d$ncont = d$ncont + 1
			contk [[d$ncont]] = paste ("x[,", i, "]", sep="")
			contl [[d$ncont]] = names (x) [i]
		}
		else stop ("unsupported variable class in interaction x")
	}
	d$kstr = interaction.combine (c (list (interaction.combine (catk, "*") ), contk), "*")
	d$labs = interaction.combine (c (list (interaction.combine (catl, ":") ), contl), ":")
	d
}

#x is a list of vectors
interaction.combine = function (x, op)
{	n = length (x)
	if (n > 1)
	{	y = x [[1]]
		for (i in 2:n)
		{	y = outer (y, x [[i]], paste, sep=op)
			y = as.vector (t (y) )
		}
		y
	}
	else x [[1]]
}

